// SolarModule.h
#pragma once

#include "DataConfig.h"
#include <string>
#include <vector>

class SolarModule {
public:
    struct Result {
        bool   on = false;
        double irr_Wm2 = 0.0;
        double mdot_kgps = 0.0;
        double T_in_C = 0.0;
        double T_out_C = 0.0;
        double Q_kW = 0.0;      // useful heat to fluid (kW)
        double Q_raw_kW = 0.0;  // irradiance*area*eta (kW), before temp limits
    };

    void initialize(const DataConfig& cfg) { cfg_ = cfg; }
    bool loadWeather();
    Result step(int hourIndex, double T_in_C) const;
    bool isLoaded() const { return loaded_; }

private:
    static bool loadCSVColumn_(const std::string& path,
                              const std::string& column,
                              std::vector<double>& out);

private:
    DataConfig cfg_{};
    std::vector<double> irr_Wm2_;
    bool loaded_ = false;
};

