#pragma once
#include "DataConfig.h"

class TankModule {
public:
    TankModule() = default;

    void initialize(const DataConfig& cfg) {
        cfg_ = cfg;
        const double T0 = cfg_.tank.setpoint_C;
        reset(T0);
    }

    void reset(double T_init_C) { T_C_ = T_init_C; }
    double temperature_C() const { return T_C_; }

    // Compute tank heat loss at current temperature (kW)
    double loss_kW() const {
        double dT = (T_C_ - cfg_.tank.ambient_T_C);
        return (dT > 0.0) ? (cfg_.tank.UA_kW_per_K * dT) : 0.0;
    }

    // Apply one time step. Q_hp_in_kW is heat added by HP to the tank (kW).
    // Q_space_req_kW & Q_dhw_req_kW are requested draws (kW) for this step.
    // Outputs the actually served draws and unmet (kW average over the step).
    void applyHour(double dt_s,
                   double Q_hp_in_kW,
                   double Q_space_req_kW,
                   double Q_dhw_req_kW,
                   double& Q_space_served_kW,
                   double& Q_dhw_served_kW,
                   double& Q_unmet_kW)
    {
        const double rho = cfg_.fluid.rho;       // kg/m3
        const double cp  = cfg_.fluid.cp;        // J/kg-K
        const double m   = rho * cfg_.tank.volume_m3; // kg
        const double Tmin = cfg_.tank.min_T_C;
        const double Tmax = cfg_.tank.setpoint_C + cfg_.tank.deadband_K; // do not store above off-threshold

        // Energy accounting in Joules
        const double E_hp   = std::max(0.0, Q_hp_in_kW) * 1000.0 * dt_s;
        const double E_loss = loss_kW() * 1000.0 * dt_s; // use current T for loss estimate
        const double E_space_req = std::max(0.0, Q_space_req_kW) * 1000.0 * dt_s;
        const double E_dhw_req   = std::max(0.0, Q_dhw_req_kW)   * 1000.0 * dt_s;

        const double E_above_min_init = std::max(0.0, (T_C_ - Tmin) * m * cp);
        double E_available = std::max(0.0, E_above_min_init + E_hp - E_loss);

        // DHW priority (typical). Serve DHW first, then space heat.

        double E_dhw_served = std::min(E_dhw_req, E_available);
        E_available -= E_dhw_served;
        double E_space_served = std::min(E_space_req, E_available);
        E_available -= E_space_served;

        Q_dhw_served_kW = E_dhw_served / dt_s / 1000.0;
        Q_space_served_kW = E_space_served / dt_s / 1000.0;

        const double E_unmet = (E_space_req + E_dhw_req) - (E_space_served + E_dhw_served);
        Q_unmet_kW = (E_unmet > 0.0) ? (E_unmet / dt_s / 1000.0) : 0.0;

        // Remaining energy above Tmin is E_available
        double T_new = Tmin + E_available / (m * cp);
        if (T_new > Tmax) T_new = Tmax;
        if (T_new < Tmin) T_new = Tmin;
        T_C_ = T_new;

    }

private:
    DataConfig cfg_;
    double T_C_ = 45.0;
};

