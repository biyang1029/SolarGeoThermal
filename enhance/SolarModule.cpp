// SolarModule.cpp
#include "SolarModule.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>

static inline double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

static inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

static inline bool ends_with_icase(const std::string& s, const std::string& suf) {
    if (s.size() < suf.size()) return false;
    for (size_t i = 0; i < suf.size(); ++i) {
        char a = s[s.size() - suf.size() + i];
        char b = suf[i];
        if (std::tolower((unsigned char)a) != std::tolower((unsigned char)b)) return false;
    }
    return true;
}

bool SolarModule::loadCSVColumn_(const std::string& path,
                                const std::string& column,
                                std::vector<double>& out)
{
    out.clear();
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;

    // EPW not supported here (no header); use Weather_gansu.csv-like format.
    if (ends_with_icase(path, ".epw")) return false;

    std::string line;
    if (!std::getline(ifs, line)) return false; // header

    std::vector<std::string> headers;
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) headers.push_back(trim(cell));
    }

    int col = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == column) { col = (int)i; break; }
    }
    if (col < 0) return false;

    struct Rec { int m = 0, d = 0, h = 0; double v = 0.0; };
    std::vector<Rec> recs;
    recs.reserve(9000);

    bool hasMDH = false;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::stringstream ls(line);
        std::string cell;
        std::vector<std::string> cells;
        while (std::getline(ls, cell, ',')) cells.push_back(trim(cell));
        if (cells.empty()) continue;

        // Prefer year,month,day,hour format (like Weather_gansu.csv)
        if (cells.size() >= 4) {
            try {
                Rec r;
                r.m = std::stoi(cells[1]);
                r.d = std::stoi(cells[2]);
                r.h = std::stoi(cells[3]);
                if (col < (int)cells.size()) r.v = std::stod(cells[col]);
                recs.push_back(r);
                hasMDH = true;
                continue;
            } catch (...) {
                // fall through to generic handling
            }
        }

        if (col < (int)cells.size()) {
            try { out.push_back(std::stod(cells[col])); } catch (...) {}
        }
    }

    if (hasMDH && !recs.empty()) {
        // Keep consistent with LoadModel: rotate to start at first 10/15 if found.
        size_t startIdx = 0;
        for (size_t i = 0; i < recs.size(); ++i) {
            if (recs[i].m == 10 && recs[i].d == 15) { startIdx = i; break; }
        }
        out.reserve(recs.size());
        for (size_t k = 0; k < recs.size(); ++k) {
            const Rec& r = recs[(startIdx + k) % recs.size()];
            out.push_back(r.v);
        }
    }

    return !out.empty();
}

bool SolarModule::loadWeather() {
    loaded_ = false;
    irr_Wm2_.clear();
    if (!cfg_.solar.enable) return false;
    loaded_ = loadCSVColumn_(cfg_.solar.weather_csv, cfg_.solar.irr_column, irr_Wm2_);
    return loaded_;
}

SolarModule::Result SolarModule::step(int hourIndex, double T_in_C) const {
    Result r;
    r.T_in_C = T_in_C;
    r.T_out_C = T_in_C;

    if (!cfg_.solar.enable) return r;
    if (!loaded_ || irr_Wm2_.empty()) return r;
    if (hourIndex < 0) return r;

    double irr = irr_Wm2_[(size_t)(hourIndex % (int)irr_Wm2_.size())];
    if (!std::isfinite(irr) || irr < 0.0) irr = 0.0;
    r.irr_Wm2 = irr;

    const double A = std::max(0.0, cfg_.solar.area_m2);
    const double eta = std::max(0.0, cfg_.solar.eta);
    const double irr_on = std::max(0.0, cfg_.solar.irr_on_Wm2);
    if (!(A > 0.0 && eta > 0.0)) return r;
    if (irr < irr_on) return r; // pump off (night/low irradiance)

    const double cp = std::max(1.0, cfg_.fluid.cp);
    const double Qraw_W = irr * A * eta;
    r.Q_raw_kW = Qraw_W / 1000.0;

    double mdot = std::max(1e-9, cfg_.solar.mdot_nominal_kgps);
    const double mdot_max = std::max(mdot, cfg_.solar.mdot_max_kgps);

    const double Tmax = cfg_.solar.T_out_max_C;
    const double dTmax = cfg_.solar.dT_max_C;

    // If inlet already at/above max outlet, do not run.
    if (std::isfinite(Tmax) && T_in_C >= Tmax - 1e-9) return r;

    // Flow modulation to avoid excessive outlet temperature and/or too-large pass temperature rise.
    if (Qraw_W > 0.0) {
        double mdot_req = mdot;
        if (std::isfinite(Tmax)) {
            double denom = cp * std::max(1e-9, (Tmax - T_in_C));
            mdot_req = std::max(mdot_req, Qraw_W / denom);
        }
        if (std::isfinite(dTmax) && dTmax > 1e-9) {
            double denom = cp * dTmax;
            mdot_req = std::max(mdot_req, Qraw_W / denom);
        }
        mdot = std::min(mdot_req, mdot_max);
    }

    // Compute outlet with chosen flow.
    double Quse_W = Qraw_W;
    double dT = Quse_W / (mdot * cp);

    // If still too hot because flow hit mdot_max, limit useful heat (defocus / spill).
    if (std::isfinite(dTmax) && dTmax > 1e-9 && dT > dTmax) {
        dT = dTmax;
        Quse_W = mdot * cp * dT;
    }
    if (std::isfinite(Tmax) && (T_in_C + dT) > Tmax) {
        dT = std::max(0.0, Tmax - T_in_C);
        Quse_W = mdot * cp * dT;
    }

    r.mdot_kgps = mdot;
    r.T_out_C = T_in_C + dT;
    r.Q_kW = Quse_W / 1000.0;
    r.on = (r.mdot_kgps > 0.0 && r.Q_kW > 1e-9);
    return r;
}

