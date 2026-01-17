// DataConfig.h
#pragma once
#include <functional>
#include <string>
#include <vector>
#include <cmath>

constexpr double PI = 3.14159265358979323846;

struct FluidProperty {
    double rho = 997.0;   // kg/m3
    double cp  = 4186.0;  // J/kg-K
    double k   = 0.6;     // W/m-K
    double mu  = 0.001;   // Pa·s
    // inlet
    double inlet_T_C      = 12.0; // C
    double massFlow_kgps  = 20.0; // kg/s
};

struct MaterialProperty {
    double k   = 4.0;     // W/m-K
    double rho = 2600.0;  // kg/m3
    double cp  = 3600.0;  // J/kg-K
};

struct WellGeometry {
    double depth_m               = 2800.0;
    int    segments              = 500;
    double D_outer_m             = 0.20;   // m
    double D_inner_m             = 0.10;   // m
    double borehole_D_m          = 0.22;   // m
    double pipe_wall_thickness_m = 0.008;  // m
};

struct SoilMesh {
    double rMax_m = 100.0;
    double dr0_m  = 0.01;
    double growth = 1.2;
};

struct TimeControl {
    double year = 10;
    double    timeStep_s  = 600;      // s
    double    totalSteps  = 6 * 24 * 10;
    int    maxIter     = 500;
    double residualTol = 1e-5;
};

// Weather-driven load configuration
struct LoadConfig {
    bool        enable_weather = true;
    std::string weather_csv    = "Weather_gansu.csv";
    std::string column_name    = "T_outdoor_C";
    double      indoor_T_C     = 17.0;
    double      UA_kW_per_K = 20;
    double      base_kW        = 5.0;
    bool        enable_cutoff  = false;  // if true, apply heat_cutoff_C rule
    double      heat_cutoff_C  = 18.0;   // no space-heating load above this Tout
    double   loadscale = 1;
};

// Domestic Hot Water (DHW) load configuration (simple schedule-based)
struct DHWConfig {
    bool   enable           = true;
    double base_kW          = 0.3;
    int    morning_start_h  = 6;
    int    morning_hours    = 2;
    double morning_kW       = 3.0;
    int    evening_start_h  = 19;
    int    evening_hours    = 2;
    double evening_kW       = 3.0;
};

// Solar thermal concentrator (very simple, lossless collector model)
struct SolarConfig {
    bool        enable = true;
    std::string weather_csv = "Weather_gansu.csv";
    std::string irr_column = "DNI";      // e.g., DNI for concentrator
    double      irr_on_Wm2 = 80.0;       // below this -> pump off
    double      area_m2 = 500.0;
    double      eta = 0.65;             // optical/thermal efficiency (lumped)
    double      mdot_nominal_kgps = 1.0;
    double      mdot_max_kgps = 50.0;   // used when limiting outlet temperature
    double      T_out_max_C = 100.0;    // do not exceed boiling risk
    double      dT_max_C = 50.0;        // limit temperature rise per pass
};

// Simple single-node buffer tank
struct TankConfig {
    double volume_m3   = 30.0;
    double setpoint_C  = 40.0;
    double deadband_K  = 30.0;
    double min_T_C     = 0.0;
    double ambient_T_C = 20.0;
    double UA_kW_per_K = 0.0002;
};

// Pump hydraulics configuration
struct PumpConfig {
    double eff            = 0.70;   // pump-wire efficiency
    double rel_roughness  = 1e-5;   // epsilon/D
    double K_minor_annulus = 5.0;   // sum of minor-loss K on annulus
    double K_minor_inner   = 5.0;   // sum of minor-loss K on inner pipe
    double mdot_ground_kgps = 10.0;
    double mdot_load_kgps = 23.0;
};

// Enhanced Heat Exchange Pipe (EHEP) configuration (single-segment quick config)
struct EnhanceConfig {
    bool   enable    = true;   // enable enhancement multipliers
    double nu_mult   = 1.25;   // multiplier on convective Nu (thus on h)
    double f_mult    = 1.41;   // multiplier on Darcy friction factor (annulus)
    double z_start_m = 0.0;    // start depth (m); if z_end_m<=z_start_m, apply to whole length
    double z_end_m   = 0.0;    // end depth (m)
};

// Multi-segment enhancement
struct EHEPSegment {
    double z_start_m = 0.0;
    double z_end_m   = 0.0;
    double nu_mult   = 1.25;
    double f_mult    = 1.41;
};

struct HeatPumpConfig {
    double Q_nom_kW = 350.0;              // 
    double design_T_source_out_C = 12.0;  // 名义工况：源侧出水（蒸发侧）输入
    double design_T_load_in_C = 45.0;   // 名义工况：负荷侧回水（冷凝器入口）
    double design_m_dot_load_kgps = 20.0; // 名义工况：负荷侧水流量（用于推 60℃时的量级）

    double defaultCOP = 3.0;
    double Q_demand_kW = 350.0;
    double target_T_load_out_C = 53.0;
    double freeze_guard_C = 1.0;
    double evap_approach_K = 5.0;
    double cond_approach_K = 5.0;
    double cond_T_max_C = 65.0;
    double min_temp_lift_K = 5.0;
    double eff_carnot = 0.55;      // Carnot效率系数，用于COP估算
    double evap_T_min_C = -60.0;   // 蒸发器最低允许温度（防冻保护参考）
    double mod_range_C = 1.0;      // 源侧温度调制区间（定频模式可忽略）
    // basic performance switches
    bool   use_coolprop = true;     // allow disabling CoolProp runtime model
    double max_Q_out_kW = 275.0;    // nameplate max heating output (kW)
    double min_source_return_C = -100.0; // minimum allowed source return temp (C)
    // For environment overrides (optional; these may be set elsewhere in project)
    std::string fluid = "R134a";
    double eta_isentropic = 0.7;
    double superheat_K = 5.0;
    double subcool_K = 3.0;
    double max_source_dT_per_hour = 0.0; // 0 = disabled cap
    double eta_motor = 0.88;
    double m_ref_nom_kgps = 0.73;
    double rho_suc_nom_kgpm3 = 17.61;
    double pressure_ratio = 4.0;


};

// Global configuration bundle
struct DataConfig {
    // geometry/mesh/time
    WellGeometry well;
    SoilMesh     mesh;
    TimeControl  time;

    // ground initial condition
    double T_surface_C      = 8.0;     // C at z=0
    double geograd_C_per_m  = 0.038;   // C per meter

    // radial mesh coarse controls (duplicated for historical reasons)
    double rMax_m    = 200.0;
    double dr0_m     = 0.005;
    double dr_growth = 1.05;

    // materials
    FluidProperty    fluid;
    MaterialProperty soil;
    MaterialProperty grout;
    MaterialProperty pipe;
    // Optional explicit inner/outer wall conductivities; if >0, override pipe.k
    double pipe_k_inner = 0.002; // W/m-K (for inner wall resistance R_w)
    double pipe_k_outer = 40.0; // W/m-K (for outer wall resistance R_p)

    // plant side
    HeatPumpConfig   hp;
    TankConfig       tank;
    DHWConfig        dhw;
    SolarConfig      solar;
    PumpConfig       pump;
    EnhanceConfig    enh;
    std::vector<EHEPSegment> enh_segments; // if empty, fall back to `enh`
    // Inner-pipe top insulation zone: apply lower equivalent conductivity within 0..top_len_m
    struct InnerInsulConfig { bool enable = false; double top_len_m = 100.0; double k_inner = 0.01; double thickness_m = 0.002; } insul;

    // Legacy simplified coupling coefficient (kept for compatibility; not used in new model)
    double h_io_Wm2K = 10.0;

    // Simulation calendar and heating season
    int sim_start_year  = 2025;
    int sim_start_month = 10;
    int sim_start_day   = 1;
    struct SeasonConfig {
        bool enable = true;
        int  start_month = 10;
        int  start_day   = 1;
        int  end_month   = 3;
        int  end_day     = 1;
    } season;

    struct EconConfig {
        double elec_price_per_kWh = 0.0; // currency/kWh
        double capex = 0.0;              // if >0, overrides breakdown
        double lifetime_years = 20.0;
        double discount_rate  = 0.08;
        double om_factor      = 0.01;    // annual O&M fraction of CAPEX
        // Breakdown (used if capex==0)
        double drill_cost_per_m = 1000.0;
        double casing_cost_per_m = 150.0;
        double ehep_cost_per_m   = 0.0;  // if 0 => use heuristic premium
        // Equipment costs
        double hp_cost_fixed     = 30000.0;
        double tank_cost_fixed   = 5000.0;
        double pump_cost_fixed   = 2000.0;
        double misc_fixed        = 5000.0;
        // Unit costs
        double hp_unit_cost_per_kW   = 200.0; // CNY/kW
        double pump_unit_cost_per_kW = 80.0;  // CNY/kW
    } econ;

    // Correlation: Nu = f(Re, Pr, z)
    std::function<double(double, double, double)> NuFunc;

    inline double calcRe(double m_dot_kgps, double D_m) const {
        // Re = 4 m_dot / (pi D mu)
        return (4.0 * m_dot_kgps) / (PI * D_m * fluid.mu);
    }
    inline double calcPr() const {
        // Pr = cp * mu / k
        return (fluid.cp * fluid.mu) / fluid.k;
    }

    // Weather/load configuration
    LoadConfig load;
};

extern DataConfig CFG;
