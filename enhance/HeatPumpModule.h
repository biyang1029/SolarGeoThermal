#pragma once
#include "DataConfig.h"


class HeatPumpModule {
public:
    HeatPumpModule() = default;
    void initialize(const DataConfig& cfg) { cfg_ = cfg; }


    // 输入：地源侧出水温（热泵蒸发器入口）、负荷侧回水温、负荷侧流量
    // 输出：回到地热回路的回水温；通过引用返回 COP、供热量Q_out(kW)、电功率P_el(kW)
    double step(double T_source_out_C, double T_load_in_C, double m_dot_load_kgps,
                double dt_s, double& COP, double& Q_out_kW, double& P_el_kW);

    struct HPDebug {
        bool used_coolprop = false;
        std::string fluid;
        double P_evap_kPa = 0.0, P_cond_kPa = 0.0;
        double T_evap_sat_K = 0.0, T_cond_sat_K = 0.0;
        double h1 = 0.0, h2s = 0.0, h2 = 0.0, h3 = 0.0;
    };
    void setDemand(double q_kW) { cfg_.hp.Q_demand_kW = q_kW; }
    void setMassFlow(double kgps) { cfg_.fluid.massFlow_kgps = kgps; }
    const HPDebug& lastDebug() const { return last_debug_; }
    void setSourceMassFlow(double m_dot_src_kgps) { m_dot_src_kgps_ = m_dot_src_kgps; }
private:
    DataConfig cfg_;
    HPDebug    last_debug_;
    bool   design_ready_ = false;
    double q_cond_design_Jperkg_ = 0.0;     // (h2-h3)_design
    double rho_suc_nom_kgpm3_ = 0.0;     // 吸气密度名义值
    double m_dot_src_kgps_ = 0.0;
};
