// ========================= SimulationController.cpp =========================
#define _CRT_SECURE_NO_WARNINGS 1
#include "SimulationController.h"
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <fstream>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <windows.h>
#endif

// local helpers (avoid C++17-specific std::clamp/std::filesystem)
static inline double clampd(double v, double lo, double hi){ return v<lo?lo:(v>hi?hi:v); }
static std::string dirname_of(const std::string& p){ size_t pos=p.find_last_of("/\\"); 
return (pos==std::string::npos)?std::string():p.substr(0,pos); }

// Progressbar
static void printProgress(long long done_min, long long total_min, int barWidth = 40) {
    if (done_min < 0) done_min = 0;
    if (done_min > total_min) done_min = total_min;
    double frac = (total_min > 0) ? (double)done_min / (double)total_min : 0.0;
    int filled = (int)std::lround(frac * barWidth);
    double done_h  = done_min  / 60.0;
    double total_h = total_min / 60.0;

    char tail[128];
    std::snprintf(tail, sizeof(tail), " %.1f%% (%.1f/%.1f h)", frac * 100.0, done_h, total_h);

    std::string bar; bar.reserve((size_t)barWidth + 3 + sizeof(tail));
    bar.push_back('[');
    for (int i = 0; i < barWidth; ++i) bar.push_back(i < filled ? '=' : ' ');
    bar.push_back(']');
    bar.append(tail);

    static size_t lastLen = 0;
    size_t pad = (lastLen > bar.size() ? lastLen - bar.size() : 0);
    // Append padding spaces to fully overwrite previous characters
    if (pad) bar.append(pad, ' ');
#ifdef _WIN32
    // Try Win32 console cursor reposition to avoid line stacking in some hosts
    static bool conInit = false; static HANDLE hOut = INVALID_HANDLE_VALUE; static SHORT baseRow = -1;
    if (!conInit) {
        hOut = GetStdHandle(STD_OUTPUT_HANDLE);
        if (hOut != INVALID_HANDLE_VALUE) {
            CONSOLE_SCREEN_BUFFER_INFO csbi; if (GetConsoleScreenBufferInfo(hOut, &csbi)) baseRow = csbi.dwCursorPosition.Y;
        }
        conInit = true;
    }
    if (hOut != INVALID_HANDLE_VALUE && baseRow >= 0) {
        COORD pos; pos.X = 0; pos.Y = baseRow;
        SetConsoleCursorPosition(hOut, pos);
        DWORD written = 0; WriteConsoleA(hOut, bar.c_str(), (DWORD)bar.size(), &written, nullptr);
        // no std::cout here to avoid mixing buffered/unbuffered outputs
    } else {
        std::cout << "\r" << bar; std::cout.flush();
    }
#else
    std::cout << "\r" << bar; std::cout.flush();
#endif
    lastLen = bar.size();
}

static void printProgressThrottled(long long done_min, long long total_min) {
    // 每小时刷新一次；如需更频繁可设 PROGRESS_STEP_PCT=1（按百分比）
    static long long lastHour = -1;
    static bool usePct = false;
    static int stepPct = 1;
    static bool inited = false;
    if (!inited) {
        if (const char* v = std::getenv("PROGRESS_STEP_PCT")) {
            try { int t = std::stoi(v); if (t > 0) { stepPct = t; usePct = true; } } catch(...) {}
        }
        inited = true;
    }
    if (usePct && stepPct > 1) {
        double frac = (total_min > 0) ? (double)done_min / (double)total_min : 0.0;
        static int lastPct = -1;
        int pct = (int)std::floor(frac * 100.0 + 1e-9);
        if ((pct / stepPct) != (lastPct / stepPct) || done_min == total_min) {
            printProgress((int)done_min, (int)total_min);
            lastPct = pct;
        }
        return;
    }
    long long hour = (done_min >= 0 ? done_min / 60 : 0);
    if (hour != lastHour || done_min == total_min) {
        printProgress((int)done_min, (int)total_min);
        lastHour = hour;
    }
}

// date helpers (non-leap year)
static int days_in_month(int m){ static const int dm[12]={31,28,31,30,31,30,31,31,30,31,30,31}; if(m<1||m>12) return 30; return dm[m-1]; }
static void add_minutes(int& y,int& m,int& d,int& hh,int& mm,int add_min){
    mm += add_min; hh += mm/60; mm%=60; d += hh/24; hh%=24;
    while(true){ int md=days_in_month(m); if(d<=md) break; d-=md; m++; if(m>12){ m=1; y++; } }
}
static int day_of_year(int y,int m,int d){ int sum=0; for(int i=1;i<m;i++) sum+=days_in_month(i); return sum+d; }
//Determine if in heating season, if no 
static bool in_heating_season(const DataConfig& cfg, int y, int m, int d) {
    if (!cfg.season.enable) return true;
    int ds = day_of_year(y, cfg.season.start_month, cfg.season.start_day);
    int de = day_of_year(y, cfg.season.end_month, cfg.season.end_day);
    int dc = day_of_year(y, m, d);
    if (de >= ds) return dc >= ds && dc <= de; // same-year window
    // wrap across year-end: Oct..Apr
    return (dc >= ds) || (dc <= de);
    // if in heating season run heat pump
}

bool SimulationController::run(const std::string& csv_path) {
    // 
    ground_.initialize(cfg_);
    hp_.initialize(cfg_);
    tank_.initialize(cfg_);
    solar_.initialize(cfg_);
    solar_.loadWeather();

    bool summaryOnly = false; if (const char* v = std::getenv("SUMMARY_ONLY")) summaryOnly = (*v!='0' && *v!='n' && *v!='N') ;
    if (!summaryOnly && !logger_.open(csv_path)) {
        std::cerr << "[Error] Failed to open log file: " << csv_path << "\n";
        return false;
    }
    // open solar results file alongside results.csv
    std::ofstream solar_ofs;
    bool solar_log = false;
    std::string solar_path;
    if (!summaryOnly) {
        auto baseName = [](const std::string& p)->std::string {
            size_t pos = p.find_last_of("/\\");
            return (pos == std::string::npos) ? p : p.substr(pos + 1);
        };
        auto replacePrefix = [](std::string s, const std::string& from, const std::string& to)->std::string {
            if (s.rfind(from, 0) == 0) s.replace(0, from.size(), to);
            return s;
        };
        std::string dir = dirname_of(csv_path);
        std::string csvBase = baseName(csv_path);
        std::string solarBase = replacePrefix(csvBase, "results", "results_solar");
        if (solarBase == csvBase) solarBase = "results_solar_" + csvBase;
        solar_path = dir.empty() ? solarBase : (dir + std::string("/") + solarBase);
        solar_ofs.open(solar_path.c_str(), std::ios::out | std::ios::trunc);
        if (solar_ofs) {
            solar_log = true;
            solar_ofs << "hour,date,time,solar_on,solar_mode,solar_irr_Wm2,solar_mdot_kgps,solar_T_in_C,solar_T_out_C,solar_Q_kW,solar_Q_raw_kW,"
                      << "T_load_supply_C,T_load_return_C,T_ground_in_C,T_ground_out_C,flow_src_kgps,flow_load_kgps,v_src_ann_mps,v_src_in_mps,"
                      << "HP_on,COP,Q_out_kW,P_el_kW,Q_geo_kW,P_pump_kW,Q_space_req_kW,Q_dhw_req_kW,Q_space_served_kW,Q_dhw_served_kW,Q_unmet_kW,T_tank_C,T_outdoor_C\n";
        }
    }
    // open debug file
    logger_.openDebug("debug.csv");

    // 
    double T_return = cfg_.fluid.inlet_T_C;
    double T_solar_return = cfg_.fluid.inlet_T_C; // used in non-heating season: ground outlet -> solar inlet
    tank_.reset(cfg_.tank.setpoint_C);
    const int    maxIter = 600;
    const double tolT    = 0.01; 

    // Variable flow configuration (min/max kg/s), default around base mass flow
    static bool flow_cfg_inited = false;
    static double flow_src_min_kgps = 0.0, flow_src_max_kgps = 0.0;
    static double flow_load_min_kgps = 0.0, flow_load_max_kgps = 0.0;
if (!flow_cfg_inited) {
    double base_src = cfg_.fluid.massFlow_kgps;
    double base_load = base_src;
    //if (const char* v = std::getenv("FLOW_SRC_KGPS")) { try { base_src = std::stod(v); } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_LOAD_KGPS")) { try { base_load = std::stod(v); } catch (...) {} }
    flow_src_min_kgps = std::max(0.1, 0.9 * base_src);
    flow_src_max_kgps = std::max(flow_src_min_kgps, 1.1 * base_src);
    flow_load_min_kgps = std::max(0.1, 1.0 * base_load);
    flow_load_max_kgps = std::max(flow_load_min_kgps, 2.0 * base_load);
    //if (const char* v = std::getenv("FLOW_SRC_MIN_KGPS")) { try { flow_src_min_kgps = std::stod(v); } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_SRC_MAX_KGPS")) { try { flow_src_max_kgps = std::stod(v); } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_LOAD_MIN_KGPS")) { try { flow_load_min_kgps = std::stod(v); } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_LOAD_MAX_KGPS")) { try { flow_load_max_kgps = std::stod(v); } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_MIN_KGPS")) { try { double t = std::stod(v); flow_src_min_kgps = t; } catch (...) {} }
    //if (const char* v = std::getenv("FLOW_MAX_KGPS")) { try { double t = std::stod(v); flow_src_max_kgps = t; } catch (...) {} }
    flow_cfg_inited = true;
}

    // Physical temperature guardrails (avoid runaway states)
    const double T_min_phys = cfg_.hp.min_source_return_C; // 硬限（默认≥0），避免被夹到负值
    const double T_max_phys = cfg_.T_surface_C + cfg_.geograd_C_per_m * cfg_.well.depth_m + 60.0;

    int curY = cfg_.sim_start_year, curM = cfg_.sim_start_month, curD = cfg_.sim_start_day, curH=0, curMin=0;
    const int stepMin = (int)std::round(cfg_.time.timeStep_s/60.0);
    bool logHourly = false; if (const char* v = std::getenv("LOG_HOURLY")) logHourly = (*v!='0' && *v!='n' && *v!='N');
    int stepsPerHour = (int)std::max(1.0, std::round(3600.0/cfg_.time.timeStep_s));

    // --- override enhanced segment depth by env var ---

        // Variable-step schedule and per-hour logging
    long long target_minutes = (long long)std::llround(cfg_.time.totalSteps * (std::max(1.0, cfg_.time.timeStep_s) / 60.0));
    long long elapsed_min = 0;
    int stepIdx = 0;
    double last_dt_s = 600.0; // default 10 min
    while (elapsed_min < target_minutes) {
        // schedule: first 2000 steps = 10 min; then +6 s per step to 1h; then keep 2h
        double dt_s = 600.0;

        if (stepIdx < 12000) dt_s = 600.0; else if (last_dt_s < 3600.0) dt_s = std::min(600.0, last_dt_s + 3.0); else dt_s = 600.0;
        last_dt_s = dt_s;
        int stepMin = (int)std::max(1.0, std::round(dt_s/60.0));

        double Q_space_req_kW = 0.0;
        bool heating_on = in_heating_season(cfg_,curY,curM,curD);
        // fail-safe: off-season forces HP off
        if ( !heating_on) hp_on_ = false; 
        if (cfg_.load.enable_weather && heating_on) {
            static bool loaded = false;
            if (!loaded) { load_.loadCSV(cfg_.load.weather_csv, cfg_.load.column_name); loaded = true; }
            int hourIndex = (int)(elapsed_min / 60);
            Q_space_req_kW = load_.demandAt(hourIndex, cfg_.load);
        }
        int hourIndex = (int)(elapsed_min / 60);
        double T_outdoor_C = load_.toutAt(hourIndex, cfg_.load.indoor_T_C);
        double Q_dhw_req_kW = (cfg_.dhw.enable ? load_.dhwAt(hourIndex, cfg_.dhw) : 0.0);

        // 负荷缩放：可按需调整
        const double load_scale = CFG.load.loadscale;
        Q_space_req_kW *= load_scale;
        Q_dhw_req_kW   *= load_scale;
        double present_req_kW = std::max(0.0, Q_space_req_kW + Q_dhw_req_kW);
        bool any_demand = (present_req_kW > 1e-6);

        const double Ttank = tank_.temperature_C();

        // Solar concentrator: heating season -> tank; non-heating season -> ground recharge loop.
        int solar_mode = heating_on ? 1 : 2; // 1=tank, 2=ground
        SolarModule::Result solar = solar_.step(hourIndex, heating_on ? Ttank : T_solar_return);
        double Q_solar_to_tank_kW = (heating_on ? solar.Q_kW : 0.0);
        // 水箱温控：<=45 开机，>=55 关机，其余保持原状态
        const double on_th  = 40.0;
        const double off_th = 45.0;
        if (!hp_on_ && Ttank <= on_th && heating_on && any_demand) hp_on_ = true;
        if (hp_on_  && (Ttank >= off_th || !heating_on || !any_demand)) hp_on_ = false;
        double present_req = std::max(0.0, Q_space_req_kW + Q_dhw_req_kW);


        double cmd = 0.0;

        if (hp_on_) hp_.setDemand(cfg_.hp.max_Q_out_kW);
        else        hp_.setDemand(0.0);

        double frac = 0.0;
        if (hp_on_) {
            double qmax = std::max(1e-6, cfg_.hp.max_Q_out_kW);
            frac = clampd(cmd / qmax, 0.0, 1.0);
        }
        // Fixed flow: in heating season use HP loop; in off-season use solar-driven ground recharge loop
        double flow_src_now = 0.0;
        double flow_load_now = 0.0;
        if (heating_on && any_demand && hp_on_) {
            flow_src_now = cfg_.pump.mdot_ground_kgps;
            flow_load_now = cfg_.pump.mdot_load_kgps;
        } else if ((!heating_on) && solar.on) {
            flow_src_now = solar.mdot_kgps;
            flow_load_now = 0.0;
        }
        ground_.setMassFlow(flow_src_now);
        hp_.setSourceMassFlow(flow_src_now);

        double COP = 0, Q_out_kW = 0, P_el_kW = 0, Q_geo_kW = 0;
        // 迭代主变量改为热泵进水温度（地源出水）
        double T_source_in_iter = clampd(heating_on ? T_return : T_solar_return, T_min_phys, T_max_phys);
        double T_geo_in_now = heating_on ? T_return : solar.T_out_C;   // 热泵出水 / 地源进水
        double T_source_out = T_source_in_iter;
        double T_return_next = T_return; // next HP-side source return temp (kept in heating season)
        double T_return_log = T_return;  // what gets written to results.csv for this step

        static bool warned_src_low = false;

        if (heating_on && any_demand && hp_on_) {
            for (int it = 0; it < maxIter - 1; ++it) {
                // 动态松弛：前期偏大加速收敛，后期偏小防震荡
                double relax_eff = 0.2;
                if (it > 30)      relax_eff = 0.08;
                else if (it > 15) relax_eff = 0.10;
                else if (it > 5)  relax_eff = 0.15;

                // 1) 热泵步：基于当前地源出水温度，得到热泵出水与 Q_geo
                T_geo_in_now = hp_.step(T_source_in_iter, tank_.temperature_C(), flow_load_now, dt_s, COP, Q_out_kW, P_el_kW);
                if (!std::isfinite(T_geo_in_now)) T_geo_in_now = T_source_in_iter;
                T_geo_in_now = clampd(T_geo_in_now, T_min_phys, T_max_phys);
                Q_geo_kW = Q_out_kW - P_el_kW;
                // 保护：源侧进水过低时提示一次
                if (!warned_src_low && hp_on_ && T_source_in_iter <= cfg_.hp.min_source_return_C + 1e-6) {
                    std::cerr << "[Warning] Source inlet temperature clamped at/near minimum (" << T_source_in_iter
                              << " C <= limit " << cfg_.hp.min_source_return_C << " C)."
                              << " Check: HP_MIN_SRC_RETURN_C/HP_FREEZE_GUARD_C, brine type (antifreeze), flow_src, bore length.\n";
                    warned_src_low = true;
                }

                // 2) 地源步：基于热泵出水温度与当前 Q_geo，预测新的地源出水温度
                double T_next_source_in = ground_.step(T_geo_in_now, Q_geo_kW, dt_s, /*advanceSoil=*/false);
                if (!std::isfinite(T_next_source_in)) T_next_source_in = T_source_in_iter;
                T_next_source_in = clampd(T_next_source_in, T_min_phys, T_max_phys);

                // 限制迭代跳变并检查收敛
                const double max_dT_iter = 5.0;
                double dT = clampd(T_next_source_in - T_source_in_iter, -max_dT_iter, max_dT_iter);
                T_next_source_in = T_source_in_iter + dT;
                if (std::fabs(T_next_source_in - T_source_in_iter) < tolT) {
                    T_source_in_iter = T_next_source_in;
                    break;
                }

                // 松弛更新
                T_source_in_iter = relax_eff * T_next_source_in + (1.0 - relax_eff) * T_source_in_iter;
                T_source_in_iter = clampd(T_source_in_iter, T_min_phys, T_max_phys);
            }

            // 用收敛后的热泵出水温度推进土壤（advanceSoil=true）
            T_source_out = ground_.step(T_geo_in_now, Q_geo_kW, dt_s, /*advanceSoil=*/true);
            if (!std::isfinite(T_source_out)) T_source_out = T_source_in_iter;
            T_source_out = clampd(T_source_out, T_min_phys, T_max_phys);
            // 最终热泵步，得到下一步回水温度
            T_return_next = hp_.step(T_source_out, tank_.temperature_C(), flow_load_now, dt_s, COP, Q_out_kW, P_el_kW);
            // 确保用于日志/后处理的 Q_geo 与最终一步一致
            Q_geo_kW = Q_out_kW - P_el_kW;
            if (!std::isfinite(T_return_next)) T_return_next = T_source_in_iter;
            bool hit_min = (T_return_next <= T_min_phys + 1e-9);
            if (hit_min) {
                // 触发源侧冻结保护：这一小时视为 HP 输出为 0（或降到很小）
                Q_out_kW = 0; P_el_kW = 0; COP = 0;
                // 不写回 clamp 值，保持上一小时回水作为状态（避免自激）
                T_return_next = T_return;
            }
            T_return_log = T_return_next;
        } else {
            // No heating-season HP operation: still advance soil; in non-heating season, solar can drive a simple recharge loop.
            double Q_ground_kW = 0.0;
            if (!heating_on) {
                if (solar.on) {
                    // Solar outlet -> ground inlet; ground outlet -> solar inlet (return) for next step.
                    T_geo_in_now = clampd(solar.T_out_C, T_min_phys, T_max_phys);
                    T_source_out = ground_.step(T_geo_in_now, Q_ground_kW, dt_s, /*advanceSoil=*/true);
                    if (!std::isfinite(T_source_out)) T_source_out = T_solar_return;
                    T_source_out = clampd(T_source_out, T_min_phys, T_max_phys);
                    T_solar_return = T_source_out;
                    // Convention: Q_geo_kW >0 means ground extraction; recharge is negative.
                    Q_geo_kW = -Q_ground_kW;
                    T_return_log = T_geo_in_now; // log ground inlet temperature in off-season
                } else {
                    // Pump off: evolve soil naturally.
                    T_source_out = ground_.step(T_solar_return, Q_ground_kW, dt_s, /*advanceSoil=*/true);
                    if (!std::isfinite(T_source_out)) T_source_out = T_solar_return;
                    T_source_out = clampd(T_source_out, T_min_phys, T_max_phys);
                    Q_geo_kW = -Q_ground_kW;
                    T_return_log = T_solar_return;
                }
            } else {
                // Heating season but HP off: evolve soil naturally.
                T_geo_in_now = clampd(T_return, T_min_phys, T_max_phys);
                T_source_out = ground_.step(T_geo_in_now, Q_ground_kW, dt_s, /*advanceSoil=*/true);
                if (!std::isfinite(T_source_out)) T_source_out = T_source_in_iter;
                T_source_out = clampd(T_source_out, T_min_phys, T_max_phys);
                Q_geo_kW = 0.0;
                T_return_log = T_return;
            }
            COP = 0.0; Q_out_kW = 0.0; P_el_kW = 0.0;
            T_return_next = T_return; // keep HP-side state unchanged outside heating season
        }

        double Q_space_served_kW = 0.0, Q_dhw_served_kW = 0.0, Q_unmet_kW = 0.0;
        tank_.applyHour(dt_s, Q_out_kW + Q_solar_to_tank_kW, Q_space_req_kW, Q_dhw_req_kW,
                        Q_space_served_kW, Q_dhw_served_kW, Q_unmet_kW);

        // 估算负荷侧供回水温度（基于水箱温度与负荷流量）
        double Q_served_kW = Q_space_served_kW + Q_dhw_served_kW;
        // 负荷侧回水即水箱温度，供水按当前输出功率和流量计算升温，热泵出水=供水
        double T_load_return = tank_.temperature_C();
        double T_load_supply = T_load_return;
        double T_hp_load_out = T_load_return;
        if (flow_load_now > 1e-9) {
            double dT_load = (Q_out_kW * 1000.0) / (flow_load_now * cfg_.fluid.cp);
            T_load_supply = T_load_return + dT_load;
            T_hp_load_out = T_load_supply;
        }

        auto safe = [](double v, double lo){ return (v>lo? v: lo); };
        const double rho = cfg_.fluid.rho;
        const double mu  = cfg_.fluid.mu;
        const double L   = std::max(0.0, cfg_.well.depth_m);
        const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
        const double D_o_out = std::max(1e-6, cfg_.well.D_outer_m);
        const double D_o_inn = std::max(1e-6, cfg_.well.D_inner_m);
        const double D_o_inner_wall = std::max(1e-6, D_o_out - 2.0 * t_p);
        const double D_i_inner = std::max(1e-6, D_o_inn - 2.0 * t_p);
        const double Dh_ann = std::max(1e-6, D_o_inn > D_o_inner_wall ? (D_o_inn - D_o_inner_wall) : (D_o_inner_wall - D_o_inn));
        const double A_ann = std::max(1e-12, 0.25 * PI * (D_o_inner_wall*D_o_inner_wall - D_o_inn*D_o_inn));
        const double A_in  = std::max(1e-12, 0.25 * PI * (D_i_inner*D_i_inner));
        const double m_dot = std::max(1e-9, flow_src_now);
        const double v_ann = m_dot / (rho * A_ann);
        const double v_in  = m_dot / (rho * A_in);
        double v_src_ann_mps = 0.0;
        double v_src_in_mps = 0.0;
        if (flow_src_now > 1e-9 && rho > 1e-9) {
            v_src_ann_mps = flow_src_now / (rho * A_ann);
            v_src_in_mps = flow_src_now / (rho * A_in);
        }
        const double Re_ann = rho * v_ann * Dh_ann / safe(mu,1e-9);
        const double Re_in  = rho * v_in  * std::max(1e-6, D_i_inner) / safe(mu,1e-9);
        auto fric = [&](double Re, double rel_rough)->double{
            if (Re < 1e-12) return 0.0;
            if (Re < 2300.0) return 64.0 / safe(Re,1e-9);
            double term = rel_rough/3.7 + 5.74/std::pow(Re,0.9);
            double f = 0.25/std::pow(std::log10(term),2.0);
            if (!std::isfinite(f) || f<=0.0) f = 0.02;
            return f;
        };
        const double f_ann = fric(Re_ann, cfg_.pump.rel_roughness);
        const double f_in  = fric(Re_in,  cfg_.pump.rel_roughness);
        const double q_ann = 0.5 * rho * v_ann * v_ann;
        const double q_in  = 0.5 * rho * v_in  * v_in;
        double f_ann_eff = f_ann;
        if (cfg_.enh.enable) {
            double frac_enh = 1.0;
            if (cfg_.enh.z_end_m > cfg_.enh.z_start_m && L > 0.0) {
                frac_enh = clampd((cfg_.enh.z_end_m - cfg_.enh.z_start_m) / L, 0.0, 1.0);
            }
            f_ann_eff = f_ann * (1.0 + frac_enh * (std::max(1.0, cfg_.enh.f_mult) - 1.0));
        }
        const double dP_ann = f_ann_eff * (L / std::max(1e-6, Dh_ann)) * q_ann + cfg_.pump.K_minor_annulus * q_ann;
        const double dP_in  = f_in  * (L / std::max(1e-6, D_i_inner)) * q_in  + cfg_.pump.K_minor_inner   * q_in;
        const double dP_total_Pa = dP_ann + dP_in;
        const double eta_p = std::max(0.05, std::min(0.95, cfg_.pump.eff));
        const double P_pump_W = (m_dot / rho) * dP_total_Pa / eta_p;
        const double P_pump_kW = P_pump_W / 1000.0;
        const double dP_kPa = dP_total_Pa / 1000.0;

        if (!summaryOnly) {
            int hours_before = (int)(elapsed_min / 60);
            int hours_after  = (int)((elapsed_min + stepMin) / 60);
            auto dbg = hp_.lastDebug();
            std::string model = dbg.used_coolprop ? "CP" : "FB";
            for (int hh = hours_before; hh < hours_after; ++hh) {
                int yL=curY, mL=curM, dL=curD, hL=curH, minL=curMin;
                int deltaMin = ((hh+1)*60) - (int)elapsed_min;
                add_minutes(yL,mL,dL,hL,minL, deltaMin);
                char bufDate[16]; snprintf(bufDate,sizeof(bufDate),"%04d-%02d-%02d",yL,mL,dL);
                char bufTime[8];  snprintf(bufTime,sizeof(bufTime), "%02d:%02d",hL,0);
                logger_.writeHour(hh, T_source_out, T_return_log, COP, Q_out_kW, P_el_kW, Q_geo_kW,
                                  T_outdoor_C, T_load_supply, T_load_return, T_hp_load_out, model,
                                  std::string(bufDate), std::string(bufTime),
                                  tank_.temperature_C(), hp_on_?1:0,
                                  Q_space_req_kW, Q_dhw_req_kW,
                                  Q_space_served_kW, Q_dhw_served_kW, Q_unmet_kW,
                                  flow_src_now, flow_load_now, dP_kPa, P_pump_kW,
                                  solar.on?1:0, solar_mode, solar.irr_Wm2, solar.mdot_kgps,
                                  solar.T_in_C, solar.T_out_C, solar.Q_kW, solar.Q_raw_kW);
                logger_.writeDebugHour(hh, dbg.fluid, dbg.used_coolprop,
                                       dbg.T_evap_sat_K, dbg.T_cond_sat_K,
                                       dbg.P_evap_kPa, dbg.P_cond_kPa,
                                       dbg.h1, dbg.h2s, dbg.h2, dbg.h3,
                                       T_source_in_iter, T_geo_in_now, Q_geo_kW,
                                       flow_src_now, flow_load_now, T_outdoor_C);
                if (solar_log) {
                    solar_ofs << hh << ',' << bufDate << ',' << bufTime << ','
                              << (solar.on?1:0) << ',' << solar_mode << ','
                              << std::fixed << std::setprecision(4)
                              << solar.irr_Wm2 << ',' << solar.mdot_kgps << ','
                              << solar.T_in_C << ',' << solar.T_out_C << ','
                              << solar.Q_kW << ',' << solar.Q_raw_kW << ','
                              << T_load_supply << ',' << T_load_return << ','
                              << T_return_log << ',' << T_source_out << ','
                              << flow_src_now << ',' << flow_load_now << ','
                              << v_src_ann_mps << ',' << v_src_in_mps << ','
                              << (hp_on_?1:0) << ',' << COP << ','
                              << Q_out_kW << ',' << P_el_kW << ',' << Q_geo_kW << ',' << P_pump_kW << ','
                              << Q_space_req_kW << ',' << Q_dhw_req_kW << ','
                              << Q_space_served_kW << ',' << Q_dhw_served_kW << ',' << Q_unmet_kW << ','
                              << tank_.temperature_C() << ',' << T_outdoor_C << "\n";
                }
            }
        }

        if (heating_on) T_return = T_return_next;
        long long done = std::min(target_minutes, elapsed_min + stepMin);
        printProgressThrottled(done, target_minutes);
        add_minutes(curY,curM,curD,curH,curMin, stepMin);
        elapsed_min += stepMin;
        stepIdx++;
    }
    // finalize progress line
    printProgress(target_minutes, target_minutes);
    std::cout << std::endl;
    logger_.close();
    if (solar_ofs.is_open()) solar_ofs.close();

    // Economic summary (optional scaling to annual)
    // Re-parse this run's results.csv
    try {
        // Determine directory from csv_path (without <filesystem>)
        std::string dir = dirname_of(csv_path);
        std::ifstream ifs(csv_path.c_str());
        std::string line; if (!ifs) throw 1;
        // header
        std::getline(ifs, line);
        double sumPel_kWh = 0.0, sumPump_kWh = 0.0, sumQspace_kWh = 0.0, sumQdhw_kWh = 0.0;
        double peakPel_kW = 0.0, peakPump_kW = 0.0, peakQout_kW = 0.0;
        int onSteps=0, starts=0, stops=0; int prevOn=0;
        // results.csv is logged per-hour, so each row represents 1 hour.
        const double dt_h = 1.0;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            // crude CSV parse
            // columns (fixed order): step,date,time,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,...,P_pump_kW,(solar...)
            // indexes used below: Q_out_kW=6, P_el_kW=7, Q_space_served_kW=18, Q_dhw_served_kW=19, HP_on=15, P_pump_kW=24
            std::vector<std::string> cells; cells.reserve(20);
            std::string cur; for (char c: line) { if (c==',') { cells.push_back(cur); cur.clear(); } else { cur.push_back(c); } } cells.push_back(cur);
            auto rd = [&](int idx){ if (idx >=0 && idx < (int)cells.size()) { try { return std::stod(cells[idx]); } catch(...) { return 0.0; } } return 0.0; };
            double Pel = rd(7), Ppump = rd(24), Qspace = rd(18), Qdhw = rd(19);
            sumPel_kWh += Pel * dt_h;
            sumPump_kWh += Ppump * dt_h;
            sumQspace_kWh += Qspace * dt_h;
            sumQdhw_kWh += Qdhw * dt_h;
            peakPel_kW = std::max(peakPel_kW, Pel);
            peakPump_kW = std::max(peakPump_kW, Ppump);
            peakQout_kW = std::max(peakQout_kW, rd(6));
            int on = int(rd(15)); if (on==1){ onSteps++; }
            if (prevOn==0 && on==1) starts++; if (prevOn==1 && on==0) stops++; prevOn = on;
        }
        double elec_kWh = sumPel_kWh + sumPump_kWh;
        double heat_kWh = sumQspace_kWh + sumQdhw_kWh;
        double cost_elec = elec_kWh * std::max(0.0, cfg_.econ.elec_price_per_kWh);
        double hoursSim = (double)cfg_.time.totalSteps * cfg_.time.timeStep_s / 3600.0;
        double scale = (hoursSim > 0.0) ? (8760.0 / hoursSim) : 0.0;
        double elec_kWh_annual = elec_kWh * scale;
        double heat_kWh_annual = heat_kWh * scale;
        double r = std::max(0.0, cfg_.econ.discount_rate);
        double N = std::max(1.0, cfg_.econ.lifetime_years);
        double CRF = (r > 0.0) ? (r * std::pow(1.0 + r, N)) / (std::pow(1.0 + r, N) - 1.0) : (1.0/N);
        // CAPEX breakdown if not provided
        double capex = cfg_.econ.capex;
        if (capex <= 0.0){
            double L = cfg_.well.depth_m;
            // EHEP total length
            double L_enh = 0.0; if (!cfg_.enh_segments.empty()){
                for (auto& s: cfg_.enh_segments){ L_enh += std::max(0.0, s.z_end_m - s.z_start_m); }
            } else if (cfg_.enh.enable){ L_enh = std::max(0.0, cfg_.enh.z_end_m - cfg_.enh.z_start_m); }
            // Borehole drilling cost by formula (CNY): 50000 + 0.206*D^2 - 50*D (D in m)
            const double Dm = L;
            double cost_drill = 50000.0 + 0.206*Dm*Dm - 50.0*Dm;
            if (!(std::isfinite(cost_drill) && cost_drill>0.0)) cost_drill = cfg_.econ.drill_cost_per_m * L;
            // Pipe/casing
            const double cost_pipe = cfg_.econ.casing_cost_per_m * L;
            // EHEP premium: explicit per-m if provided else 5% over casing per m on enhanced length
            double cost_ehep = (cfg_.econ.ehep_cost_per_m > 0.0) ? (cfg_.econ.ehep_cost_per_m * L_enh)
                                                                 : (0.05 * cfg_.econ.casing_cost_per_m * L_enh);
            // HP & pump equipment
            double cost_hp = (cfg_.econ.hp_unit_cost_per_kW > 0.0) ? (cfg_.econ.hp_unit_cost_per_kW * peakQout_kW)
                                                                   : cfg_.econ.hp_cost_fixed;
            double cost_pump = (cfg_.econ.pump_unit_cost_per_kW > 0.0) ? (cfg_.econ.pump_unit_cost_per_kW * peakPump_kW)
                                                                       : cfg_.econ.pump_cost_fixed;
            const double cost_tank = cfg_.econ.tank_cost_fixed;
            const double cost_misc = cfg_.econ.misc_fixed;
            capex = cost_drill + cost_pipe + cost_ehep + cost_hp + cost_pump + cost_tank + cost_misc;
        }
        const double om_annual = cfg_.econ.om_factor * capex;
        double LCOH_annual = (heat_kWh_annual > 1e-9) ? ((capex * CRF) + om_annual + elec_kWh_annual * cfg_.econ.elec_price_per_kWh) / heat_kWh_annual : 0.0;
        //std::string sumPath = dir.empty()? std::string("summary.csv") : (dir + std::string("/summary.csv"));
        // derive summary file name from csv_path: results_XXX.csv -> summary_XXX.csv
        auto baseName = [](const std::string& p)->std::string {
            size_t pos = p.find_last_of("/\\");
            return (pos == std::string::npos) ? p : p.substr(pos + 1);
            };
        auto replacePrefix = [](std::string s, const std::string& from, const std::string& to)->std::string {
            if (s.rfind(from, 0) == 0) s.replace(0, from.size(), to);
            return s;
            };
        std::string csvBase = baseName(csv_path);
        std::string sumBase = replacePrefix(csvBase, "results", "summary");
        if (sumBase == csvBase) sumBase = "summary_" + csvBase; // fallback

        std::string sumPath = dir.empty() ? sumBase : (dir + std::string("/") + sumBase);


        std::ofstream ofs(sumPath.c_str(), std::ios::out | std::ios::trunc);
        if (ofs) {
            ofs << "elec_kWh,comp_kWh,pump_kWh,heat_kWh,space_kWh,dhw_kWh,cost_elec,capex,lifetime_y,discount,CRF,LCOH_annual,SCOP_sys,peak_Pel_kW,peak_Ppump_kW,peak_Qout_kW,on_hours,starts,stops,om_annual\n";
            ofs << (elec_kWh) << ',' << (sumPel_kWh) << ',' << (sumPump_kWh) << ',' << (heat_kWh) << ','
                << (sumQspace_kWh) << ',' << (sumQdhw_kWh) << ',' << (cost_elec) << ','
                << (capex) << ',' << (cfg_.econ.lifetime_years) << ',' << (cfg_.econ.discount_rate) << ',' << (CRF) << ',' << (LCOH_annual) << ','
                << ((elec_kWh>1e-9)? (heat_kWh/elec_kWh):0.0) << ',' << peakPel_kW << ',' << peakPump_kW << ',' << peakQout_kW << ','
                << (onSteps * dt_h) << ',' << starts << ',' << stops << ',' << om_annual << "\n";
        }
    } catch(...) {
        // ignore errors
    }
    std::cout << std::endl;
    std::cout << "Simulation finished. Results at: " << csv_path << "\n";
    return true;
}
