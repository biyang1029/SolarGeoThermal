// GroundModule.h
#pragma once

#include "DataConfig.h"
#include <vector>

class GroundModule {
public:
    void initialize(const DataConfig& cfg);
    // Returns outlet temperature of inner leg (top), fills extracted kW
    double step(double T_return_C, double& Q_extracted_kW, double dt_s, bool advanceSoil);
    void setMassFlow(double m_dot_kgps);

private:
    void update_outer_downstream_(double T_inlet_C, const std::vector<double>& T_soil_local);
    void update_inner_upstream_(double T_bottom_in_C, const std::vector<double>& T_outer_local);
    double integrate_Q_extracted_kW_(const std::vector<double>& q_wall_Wm2) const;

    void build_soil_grid_();
    void init_soil_field_();
    void soil_step_ADI_(double dt_s);
    void sweep_r_(double dt_s);
    void sweep_z_(double dt_s);
    void thomas_solve_(const std::vector<double>& a,
                       const std::vector<double>& b,
                       const std::vector<double>& c,
                       std::vector<double>& d,
                       std::vector<double>& x);

private:
    DataConfig cfg_{};
    int N_  = 0;
    int Nz_ = 0;
    int Nr_ = 0;
    double dz_ = 0.0;
    double r_bore_ = 0.0;
    double m_dot_kgps_ = 0.0;
    std::vector<double> zc_;
    std::vector<double> z_;
    std::vector<double> r_;
    std::vector<double> dr_;

    std::vector<double> T_outer_;
    std::vector<double> T_inner_;
    std::vector<std::vector<double>> Tsoil_;
    std::vector<std::vector<double>> Ttmp_;
    std::vector<double> q_line_wall_;

    static constexpr double PI = 3.14159265358979323846;
};
