// GroundModule.cpp
#include "GroundModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>

// helper clamp
static inline double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

void GroundModule::setMassFlow(double m_dot_kgps) {
    m_dot_kgps_ = m_dot_kgps;
}

void GroundModule::initialize(const DataConfig& cfg) {
    cfg_ = cfg;
    N_ = cfg_.well.segments;
    dz_ = cfg_.well.depth_m / static_cast<double>(N_);

    zc_.resize(N_);
    T_outer_.assign(N_, cfg_.fluid.inlet_T_C);
    T_inner_.assign(N_, cfg_.fluid.inlet_T_C);
    for (int i = 0; i < N_; ++i) zc_[i] = (i + 0.5) * dz_;

    double D_bore = (cfg_.well.borehole_D_m > 0.0 ? cfg_.well.borehole_D_m : cfg_.well.D_outer_m);
    r_bore_ = std::max(1e-6, 0.5 * D_bore);

    Nz_ = N_;
    z_.resize(Nz_);
    for (int i = 0; i < Nz_; ++i) z_[i] = zc_[i];

    build_soil_grid_();
    init_soil_field_();

    double Twall0 = Tsoil_.size() ? Tsoil_[0][0] : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[0]);
    std::fill(T_outer_.begin(), T_outer_.end(), Twall0);
    std::fill(T_inner_.begin(), T_inner_.end(), Twall0);

    q_line_wall_.assign(Nz_, 0.0);
}

double GroundModule::step(double T_return_C, double& Q_extracted_kW, double dt_s, bool advanceSoil) {
    // 近壁土壤温度（用于对流边界）
    std::vector<double> T_soil_local(N_);
    for (int i = 0; i < N_; ++i) {
        T_soil_local[i] = (Nr_ > 0) ? Tsoil_[0][i] : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[i]);
    }
    // ---- ADD: pump-off short-circuit ----
    if (m_dot_kgps_ <= 0.0) {
        Q_extracted_kW = 0.0;

        // 可选：如果你希望“热泵停机时土壤也按自然边界继续演化”，保留这一句
        if (advanceSoil) soil_step_ADI_(dt_s);

        // 不改变管内温度（保持上一步状态）
        return T_inner_.empty() ? T_return_C : T_inner_.front();
    }
    // ---- END ADD ----

    // ---- REPLACE: keep numeric safety but not fake-flow ----
    const double m_dot = m_dot_kgps_;   // 这里不再 max(1e-9,...)
    // 如果你仍担心除零，在后面用 std::max(1e-12, m_dot) 只做分母保护即可

    const double cp    = cfg_.fluid.cp;
    const double rho   = cfg_.fluid.rho;
    const double mu    = std::max(1e-9, cfg_.fluid.mu);
    const double kf    = std::max(1e-9, cfg_.fluid.k);

    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double D_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double D_inn = std::max(1e-6, cfg_.well.D_inner_m);
    const double D_inner_wall = std::max(1e-6, D_out - 2.0 * t_p);
    const double D_i_inner = std::max(1e-6, D_inn - 2.0 * t_p);
    const double Dh_outer = std::max(1e-6, D_out - D_inn);
    const double Dh_ann   = std::max(1e-6, std::abs(D_inner_wall - D_inn));
    const double P_outer = PI * D_out;
    const double P_inner = PI * D_i_inner;
    const double r_po = 0.5 * D_out;
    const double r_pi = std::max(1e-6, r_po - t_p);
    const double r_bo = std::max(r_bore_, r_po + 1e-6);

    auto Nu_turb = [](double Re, double Pr){ return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, 0.4); };
    auto Nu_lam  = [](double, double){ return 4.36; };
    auto blendNu = [&](double Re, double Pr){
        if (Re < 2000.0) return Nu_lam(Re,Pr);
        if (Re > 3000.0) return Nu_turb(Re,Pr);
        double a = (Re - 2000.0) / 1000.0;
        return (1.0 - a) * Nu_lam(Re,Pr) + a * Nu_turb(Re,Pr);
    };
    const double Pr = cfg_.calcPr();
    const double h_min = 5.0;
    const double h_boost = 1.0; // 适度增强对流换热系数

    std::vector<double> T_outer_new(N_, 0.0), T_inner_new(N_, 0.0), q_line(N_, 0.0);
    const int maxIter = 6;
    const double tolT = 0.01;

    for (int it = 0; it < maxIter; ++it) {
        // 外管：从上到下
        double T_prev = T_return_C;
        double Re_outer = cfg_.calcRe(m_dot, Dh_outer);
        for (int i = 0; i < N_; ++i) {
            double Nu = cfg_.NuFunc ? cfg_.NuFunc(Re_outer, Pr, zc_[i]) : blendNu(Re_outer, Pr);
            double h_o  = h_boost * std::max(h_min, Nu * kf / std::max(1e-6, Dh_outer));
      
            // 仅计算管壁和灌浆的导热热阻 R_cond (即原 Rp + Rg)
            double k_outer = (cfg_.pipe_k_outer > 0.0 ? cfg_.pipe_k_outer : cfg_.pipe.k);
            double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, k_outer)); // 外管壁导热热阻
            double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k)); // 灌浆导热热阻
            double R_cond = Rp + Rg; // R_cond：从外管内壁到近壁土壤的总导热热阻

            double T_wall_soil = T_soil_local[i]; // 近壁土壤温度 T_soil

            // (1) 定义 q_wall_conv: 外管流体与外管内壁的对流换热热流
            // (2) 定义 q_wall_cond: 外管内壁通过导热传给土壤的热流
            // 在稳定状态下，q_wall_conv = q_wall_cond = q_wall_total

            // T_wall_inner 是未知的“外管内壁温度”
            // q_wall_conv = h_o * P_outer * (T_prev - T_wall_inner)
            // q_wall_cond = (T_wall_inner - T_wall_soil) / R_cond
            // 联立解 T_wall_inner:
            // T_wall_inner * (h_o * P_outer + 1.0 / R_cond) = h_o * P_outer * T_prev + T_wall_soil / R_cond

            double hA_conv = h_o * P_outer;
            double G_cond = 1.0 / std::max(1e-12, R_cond); // 导热热导
            double T_wall_inner = (hA_conv * T_prev + G_cond * T_wall_soil) / std::max(1e-12, (hA_conv + G_cond));

            // (3) 用已求得的 T_wall_inner 计算实际的壁面换热热流 q_wall
            double q_wall = hA_conv * (T_prev - T_wall_inner); // W/m
            // --- END OF MODIFICATION ---

            // 约 111 行：原代码中的 q_wall 现在已经由新的方法计算得到
            // double q_wall = (T_prev - T_wall) / std::max(1e-12, Rb); // DELETE THIS LINE

            // 与内管耦合（用上一轮的内管温度）
            // 与内管耦合（用上一轮的内管温度）
            double T_inner_old = T_inner_.empty() ? T_prev : T_inner_[i];
            double q_inner = 0.0;
            {
                // ---- 1) 速度/雷诺数：仍沿用你原来的写法 ----
                double A_ann = std::max(1e-12, 0.25 * PI * (D_inner_wall * D_inner_wall - D_inn * D_inn));
                double A_in = std::max(1e-12, 0.25 * PI * (D_i_inner * D_i_inner));
                double v_in = m_dot / (rho * A_in);
                double v_ann = m_dot / (rho * A_ann);
                double Re_in = rho * v_in * std::max(1e-6, D_i_inner) / mu;
                double Re_ann = rho * v_ann * std::max(1e-6, Dh_ann) / mu;

                double Nu_in = blendNu(Re_in, Pr);
                double Nu_ann = blendNu(Re_ann, Pr);

                // 内管内对流换热系数（inner fluid -> inner pipe inner wall）
                double h_in = std::max(h_min, Nu_in * kf / std::max(1e-6, D_i_inner));

                // 环空侧对流换热系数（annulus fluid -> inner pipe outer wall）
                // 注意：这里特征长度你原来用 Dh_ann，我保留
                double h_ann = std::max(h_min, Nu_ann * kf / std::max(1e-6, Dh_ann));

                // ---- 2) 单位长度周长：一定要用“对应的传热表面” ----
                // 内管内表面周长
                const double P_in = PI * std::max(1e-6, D_i_inner);
                // 内管外表面周长（环空侧换热就是对着内管外壁）
                const double P_ann = PI * std::max(1e-6, D_inn);

                // ---- 3) 单位长度热阻网络 R' (K·m/W) ----
                // (a) 内管内对流热阻
                const double R_conv_in = 1.0 / std::max(1e-12, h_in * P_in);

                // (b) 内管管壁导热热阻（从内管内半径 -> 内管外半径）
                // 半径：r_i = D_i_inner/2, r_o = D_inn/2
                double k_pipe_inner = (cfg_.pipe_k_inner > 0.0 ? cfg_.pipe_k_inner : cfg_.pipe.k);
                const double r_i = 0.5 * std::max(1e-6, D_i_inner);
                const double r_o = 0.5 * std::max(1e-6, D_inn);
                const double R_cond_pipe = std::log(std::max(1e-12, r_o / r_i)) / (2.0 * PI * std::max(1e-12, k_pipe_inner));

                double R_cond_ins = 0.0;
                // 只在保温层启用 + 深度范围内才加
                if (cfg_.insul.enable && zc_[i] <= std::max(0.0, cfg_.insul.top_len_m)) {
                    const double t_ins = cfg_.insul.thickness_m;
                    const double k_ins = cfg_.insul.k_inner;

                    if (t_ins > 1e-9 && k_ins > 1e-9) {
                        const double r_ins_i = r_o;
                        const double r_ins_o = r_o + t_ins;
                        R_cond_ins =
                            std::log(std::max(1e-12, r_ins_o / r_ins_i)) /
                            (2.0 * PI * k_ins);
                    }
                }


                // (d) 环空侧对流热阻（annulus fluid -> inner pipe outer wall）
                const double R_conv_ann = 1.0 / std::max(1e-12, h_ann * P_ann);

                // 总热阻（单位长度）
                const double R_tot = R_conv_in + R_cond_pipe + R_cond_ins + R_conv_ann;

                // ---- 4) 单位长度换热 q' (W/m) ----
                // 保持你原来的符号：T_prev 是环空温度，T_inner_old 是内管流体温度
                q_inner = (T_prev - T_inner_old) / std::max(1e-12, R_tot);
            }


            double q_prime = q_wall + q_inner; // 总换热
            double dT = - (q_prime * dz_) / (m_dot * cp);
            dT = clampd(dT, -55.0, 55.0);
            T_outer_new[i] = T_prev + dT;
            q_line[i] = q_wall; // 壁面传给土壤的热流，供后续土壤更新
            T_prev = T_outer_new[i];
        }

        // 内管：从下到上，耦合外管温度
        double A_ann  = std::max(1e-12, 0.25 * PI * (D_inner_wall*D_inner_wall - D_inn*D_inn));
        double A_in   = std::max(1e-12, 0.25 * PI * (D_i_inner*D_i_inner));
        double v_in  = m_dot / (rho * A_in);
        double v_ann = m_dot / (rho * A_ann);
        double Re_in  = rho * v_in  * std::max(1e-6, D_i_inner)    / mu;
        double Re_ann = rho * v_ann * std::max(1e-6, Dh_ann) / mu;
        double Nu_in  = blendNu(Re_in, Pr);
        double Nu_ann = blendNu(Re_ann,Pr);
        double h_i = h_boost * std::max(h_min, Nu_in  * kf / std::max(1e-6, D_i_inner));
        double h_o = h_boost * std::max(h_min, Nu_ann * kf / std::max(1e-6, Dh_ann));

        T_prev = T_outer_new.back();
        for (int i = N_ - 1; i >= 0; --i) {
            // 用与下行段一致的热阻网络来算内外管换热（单位长度）
            double A_ann2 = std::max(1e-12, 0.25 * PI * (D_inner_wall * D_inner_wall - D_inn * D_inn));
            double A_in2 = std::max(1e-12, 0.25 * PI * (D_i_inner * D_i_inner));
            double v_in2 = m_dot / (rho * A_in2);
            double v_ann2 = m_dot / (rho * A_ann2);
            double Re_in2 = rho * v_in2 * std::max(1e-6, D_i_inner) / mu;
            double Re_ann2 = rho * v_ann2 * std::max(1e-6, Dh_ann) / mu;

            double Nu_in2 = blendNu(Re_in2, Pr);
            double Nu_ann2 = blendNu(Re_ann2, Pr);

            double h_in2 = h_boost * std::max(h_min, Nu_in2 * kf / std::max(1e-6, D_i_inner));
            double h_ann2 = h_boost * std::max(h_min, Nu_ann2 * kf / std::max(1e-6, Dh_ann));

            const double P_in2 = PI * std::max(1e-6, D_i_inner);
            const double P_ann2 = PI * std::max(1e-6, D_inn);

            const double R_conv_in2 = 1.0 / std::max(1e-12, h_in2 * P_in2);
            const double R_conv_ann2 = 1.0 / std::max(1e-12, h_ann2 * P_ann2);

            double k_pipe_inner2 = (cfg_.pipe_k_inner > 0.0 ? cfg_.pipe_k_inner : cfg_.pipe.k);
            const double r_i2 = 0.5 * std::max(1e-6, D_i_inner);
            const double r_o2 = 0.5 * std::max(1e-6, D_inn);
            const double R_cond_pipe2 =
                std::log(std::max(1e-12, r_o2 / r_i2)) / (2.0 * PI * std::max(1e-12, k_pipe_inner2));

            // 绝热层（与下行段保持一致：先用 1A 的写法）
            double R_cond_ins2 = 0.0;

            if (cfg_.insul.enable && zc_[i] <= std::max(0.0, cfg_.insul.top_len_m)) {
                const double t_ins2 = cfg_.insul.thickness_m;
                const double k_ins2 = cfg_.insul.k_inner;

                if (t_ins2 > 1e-9 && k_ins2 > 1e-9) {
                    const double r_ins_i2 = r_o2;
                    const double r_ins_o2 = r_o2 + t_ins2;
                    R_cond_ins2 =
                        std::log(std::max(1e-12, r_ins_o2 / r_ins_i2)) /
                        (2.0 * PI * k_ins2);
                }
            }


            const double R_tot2 = R_conv_in2 + R_cond_pipe2 + R_cond_ins2 + R_conv_ann2;

            // 上行段：T_prev 是内管流体温度，T_outer_new[i] 是环空温度
            double q_prime = (T_prev - T_outer_new[i]) / std::max(1e-12, R_tot2); // W/m
            double dT = -(q_prime * dz_) / (m_dot * cp);

            dT = clampd(dT, -15.0, 15.0);
            T_inner_new[i] = T_prev + dT;
            T_prev = T_inner_new[i];
        }

        // 收敛检查
        double maxDiff = 0.0;
        for (int i = 0; i < N_; ++i) {
            maxDiff = std::max(maxDiff, std::abs(T_outer_new[i] - T_outer_[i]));
            maxDiff = std::max(maxDiff, std::abs(T_inner_new[i] - T_inner_[i]));
        }
        T_outer_ = T_outer_new;
        T_inner_ = T_inner_new;
        q_line_wall_ = q_line;
        if (maxDiff < tolT) break;
    }

 
    if (advanceSoil) soil_step_ADI_(dt_s);

    double Qtot_W = 0.0;
    for (int i = 0; i < N_; ++i) Qtot_W += q_line_wall_[i] * dz_;
    Q_extracted_kW = Qtot_W / 1000.0;
    return T_inner_.front();
}

void GroundModule::update_outer_downstream_(double T_inlet_C, const std::vector<double>& T_soil_local) {
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp = cfg_.fluid.cp;

    const double Dh_outer = std::max(1e-6, cfg_.well.D_outer_m - cfg_.well.D_inner_m);
    const double Re_outer = cfg_.calcRe(m_dot, Dh_outer);
    const double Pr = cfg_.calcPr();

    const double D_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double r_po = 0.5 * D_out;
    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double r_pi = std::max(1e-6, r_po - t_p);
    const double r_bo = std::max(r_bore_, r_po + 1e-6);
    const double P_f = PI * D_out;

    double T_prev = T_inlet_C;
    for (int i = 0; i < N_; ++i) {
        const double z = zc_[i];
        const double Nu = cfg_.NuFunc ? cfg_.NuFunc(Re_outer, Pr, z) : (0.023 * std::pow(Re_outer, 0.8) * std::pow(Pr, 0.3));
        const double h = Nu * cfg_.fluid.k / D_out;

        const double Rf = 1.0 / std::max(1e-9, h * P_f);
        const double k_outer = (cfg_.pipe_k_outer > 0.0 ? cfg_.pipe_k_outer : cfg_.pipe.k);
        const double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, k_outer));
        const double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k));
        const double Rtot = Rf + Rp + Rg;

        const double T_wall = T_soil_local[i];
        const double q_prime = (T_prev - T_wall) / std::max(1e-12, Rtot);
        const double dT = - (q_prime * dz_) / (m_dot * cp);
        T_outer_[i] = T_prev + dT;
        T_prev = T_outer_[i];
    }
}

void GroundModule::update_inner_upstream_(double T_bottom_in_C, const std::vector<double>& T_outer_local) {
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp    = cfg_.fluid.cp;
    const double rho   = cfg_.fluid.rho;
    const double mu    = std::max(1e-9, cfg_.fluid.mu);
    const double kf    = std::max(1e-9, cfg_.fluid.k);

    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double D_o_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double D_o_inn = std::max(1e-6, cfg_.well.D_inner_m);
    const double D_o_inner_wall = std::max(1e-6, D_o_out - 2.0 * t_p);
    const double D_i_inner = std::max(1e-6, D_o_inn - 2.0 * t_p);
    const double r_i_inner = 0.5 * D_i_inner;
    const double r_i_outer = 0.5 * D_o_inn;

    const double Dh_ann = std::max(1e-6, std::abs(D_o_inner_wall - D_o_inn));
    const double A_ann  = std::max(1e-12, 0.25 * PI * (D_o_inner_wall*D_o_inner_wall - D_o_inn*D_o_inn));
    const double A_in   = std::max(1e-12, 0.25 * PI * (D_i_inner*D_i_inner));

    auto Nu_turb = [](double Re, double Pr){ return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, 0.4); };
    auto Nu_lam  = [](double, double){ return 4.36; };
    auto blendNu = [&](double Re, double Pr){
        if (Re < 2000.0) return Nu_lam(Re,Pr);
        if (Re > 3000.0) return Nu_turb(Re,Pr);
        double a = (Re - 2000.0) / 1000.0; return (1.0 - a) * Nu_lam(Re,Pr) + a * Nu_turb(Re,Pr);
    };

    const double Pr = cfg_.calcPr();
    const double h_min = 5.0;

    double T_prev = T_bottom_in_C;
    for (int i = N_ - 1; i >= 0; --i) {
        const double T_outer = T_outer_local[i];
        const double v_in  = m_dot / (rho * A_in);
        const double v_ann = m_dot / (rho * A_ann);
        const double Re_in  = rho * v_in  * std::max(1e-6, D_i_inner) / mu;
        const double Re_ann = rho * v_ann * std::max(1e-6, Dh_ann)    / mu;
        const double Nu_in  = blendNu(std::max(1e-6, Re_in),  std::max(1e-6, Pr));
        const double Nu_ann = blendNu(std::max(1e-6, Re_ann), std::max(1e-6, Pr));
        const double h_i = std::max(h_min, Nu_in  * kf / std::max(1e-6, D_i_inner));
        const double h_o = std::max(h_min, Nu_ann * kf / std::max(1e-6, Dh_ann));

        const double R_i = 1.0 / std::max(1e-9, h_i * 2.0 * PI * std::max(1e-9, r_i_inner));

        // 管壁导热永远用管材
        double k_pipe = (cfg_.pipe_k_inner > 0.0 ? cfg_.pipe_k_inner : cfg_.pipe.k);
        const double R_w =
            std::log(std::max(1e-12, r_i_outer / r_i_inner)) /
            (2.0 * PI * std::max(1e-12, k_pipe));

        // 保温层作为额外热阻（深度范围内启用）
        double R_ins = 0.0;
        double r_conv_outer = r_i_outer; // 默认对流发生在内管外表面

        if (cfg_.insul.enable && zc_[i] <= std::max(0.0, cfg_.insul.top_len_m)) {
            const double t_ins = cfg_.insul.thickness_m;
            const double k_ins = cfg_.insul.k_inner;

            if (t_ins > 1e-9 && k_ins > 1e-9) {
                const double r_ins_i = r_i_outer;
                const double r_ins_o = r_i_outer + t_ins;

                R_ins =
                    std::log(std::max(1e-12, r_ins_o / r_ins_i)) /
                    (2.0 * PI * k_ins);

                r_conv_outer = r_ins_o; // 对流边界在保温层外表面
            }
        }

        // 环空侧对流（用保温层外半径）
        const double R_o =
            1.0 / std::max(1e-12, h_o * 2.0 * PI * std::max(1e-12, r_conv_outer));

        const double U_per_len =
            1.0 / std::max(1e-12, R_i + R_w + R_ins + R_o);

        const double q_prime_Wpm = U_per_len * (T_prev - T_outer); // W/m
        const double dT_raw = - (q_prime_Wpm * dz_) / (m_dot * cp);
        double dT = dT_raw; if (dT > 10.0) dT = 10.0; else if (dT < -10.0) dT = -10.0;
        T_inner_[i] = T_prev + dT; T_prev = T_inner_[i];
    }
}

double GroundModule::integrate_Q_extracted_kW_(const std::vector<double>& q_wall_Wm2) const {
    const double D_ref = std::max(1e-6, cfg_.well.D_outer_m);
    const double P_outer = PI * D_ref;
    double Q_total_W = 0.0;
    for (int i = 0; i < N_; ++i) Q_total_W += q_wall_Wm2[i] * (P_outer * dz_);
    return Q_total_W / 1000.0;
}

void GroundModule::build_soil_grid_() {
    r_.clear(); dr_.clear();
    const double rMax = cfg_.rMax_m; double r_acc = r_bore_; double dr = cfg_.dr0_m;
    while (r_acc < rMax) {
        double r_cell_center = r_acc + 0.5 * dr; 
        r_.push_back(r_cell_center); 
        dr_.push_back(dr); 
        r_acc += dr; 
        dr *= cfg_.dr_growth; }
    Nr_ = static_cast<int>(r_.size());
    if (Nr_ < 3) { while (Nr_ < 3) { double back = dr_.empty() ? cfg_.dr0_m : dr_.back(); r_.push_back((r_acc + 0.5 * back)); dr_.push_back(back); r_acc += back; ++Nr_; } }
}

void GroundModule::init_soil_field_() {
    Tsoil_.assign(Nr_, std::vector<double>(Nz_, 0.0));
    Ttmp_.assign(Nr_, std::vector<double>(Nz_, 0.0));
    for (int iz = 0; iz < Nz_; ++iz) { double Tz = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[iz]; for (int ir = 0; ir < Nr_; ++ir) Tsoil_[ir][iz] = Tz; }
}

void GroundModule::soil_step_ADI_(double dt_s) {
    if (Nr_ <= 0 || Nz_ <= 0) return; double half = std::max(1e-9, 0.5 * dt_s);
    sweep_r_(half); sweep_z_(half);
}

void GroundModule::sweep_r_(double dt_s) {
    const double rho = cfg_.soil.rho; const double cp = cfg_.soil.cp; const double k = cfg_.soil.k;
    #pragma omp parallel for schedule(static)
    for (int iz = 0; iz < Nz_; ++iz) {
        std::vector<double> a(Nr_, 0.0), b(Nr_, 0.0), c(Nr_, 0.0), d(Nr_, 0.0), x(Nr_, 0.0);
        for (int ir = 0; ir < Nr_; ++ir) Ttmp_[ir][iz] = Tsoil_[ir][iz];
        for (int ir = 0; ir < Nr_; ++ir) {
            double r = r_[ir]; double dr = dr_[ir]; double h_z = (iz == 0 || iz == Nz_ - 1) ? ((Nz_ > 1) ? (z_[std::min(iz + 1, Nz_ - 1)] - z_[std::max(iz - 1, 0)]) / 2.0 : 10.0) : (z_[iz + 1] - z_[iz - 1]) * 0.5;
            double area_r = 2.0 * PI * r * h_z; double vol = area_r * dr; double cap = rho * cp * vol / dt_s;
            double km = k, kp = k; double drm = (ir == 0) ? dr : 0.5 * (dr_[ir - 1] + dr_[ir]); double drp = (ir == Nr_ - 1) ? dr : 0.5 * (dr_[ir] + dr_[ir + 1]);
            double aW = (ir > 0 ? km * area_r / drm : 0.0); double aE = (ir < Nr_ - 1 ? kp * area_r / drp : 0.0);
            a[ir] = -aW; c[ir] = -aE; b[ir] = cap + aW + aE; d[ir] = cap * Tsoil_[ir][iz];
        }
        double q_line = q_line_wall_[iz]; double q_face = q_line / (2.0 * PI * std::max(1e-6, r_bore_));
        double h_z2 = (iz == 0 || iz == Nz_ - 1) ? ((Nz_ > 1) ? (z_[std::min(iz + 1, Nz_ - 1)] - z_[std::max(iz - 1, 0)]) / 2.0 : 10.0) : (z_[iz + 1] - z_[iz - 1]) * 0.5;
        double area_face = 2.0 * PI * std::max(1e-6, r_bore_) * h_z2; d[0] += q_face * area_face;
        double T_far = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[iz]; a[Nr_ - 1] = 0.0; c[Nr_ - 1] = 0.0; b[Nr_ - 1] = 1.0; d[Nr_ - 1] = T_far;
        thomas_solve_(a, b, c, d, x); for (int ir = 0; ir < Nr_; ++ir) Tsoil_[ir][iz] = x[ir];
    }
}

void GroundModule::sweep_z_(double dt_s) {
    const double rho = cfg_.soil.rho; const double cp = cfg_.soil.cp; const double k = cfg_.soil.k; const double dz = dz_;
    #pragma omp parallel for schedule(static)
    for (int ir = 0; ir < Nr_; ++ir) {
        double r = r_[ir]; double dr = dr_[ir]; double area_z = 2.0 * PI * r * dr; double vol = area_z * dz;
        std::vector<double> a(Nz_, 0.0), b(Nz_, 0.0), c(Nz_, 0.0), d(Nz_, 0.0), x(Nz_, 0.0);
        for (int iz = 0; iz < Nz_; ++iz) { double cap = rho * cp * vol / dt_s; double aS = (iz > 0 ? k * area_z / dz : 0.0); double aN = (iz < Nz_ - 1 ? k * area_z / dz : 0.0); a[iz] = -aS; c[iz] = -aN; b[iz] = cap + aS + aN; d[iz] = cap * Tsoil_[ir][iz]; }
        if (Nz_ >= 1) { double T_top = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[0]; a[0] = 0.0; c[0] = 0.0; b[0] = 1.0; d[0] = T_top; }
        if (Nz_ >= 2) { int izb = Nz_ - 1; double T_bot = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[izb]; a[izb] = 0.0; c[izb] = 0.0; b[izb] = 1.0; d[izb] = T_bot; }
        thomas_solve_(a, b, c, d, x); for (int iz = 0; iz < Nz_; ++iz) Tsoil_[ir][iz] = x[iz];
    }
}

void GroundModule::thomas_solve_(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, std::vector<double>& d, std::vector<double>& x) {
    const int n = (int)b.size(); std::vector<double> c2(n, 0.0), d2(n, 0.0); c2[0] = c[0] / std::max(1e-12, b[0]); d2[0] = d[0] / std::max(1e-12, b[0]);
    for (int i = 1; i < n; ++i) { double denom = std::max(1e-12, b[i] - a[i] * c2[i - 1]); c2[i] = (i < n - 1) ? c[i] / denom : 0.0; d2[i] = (d[i] - a[i] * d2[i - 1]) / denom; }
    x[n - 1] = d2[n - 1]; for (int i = n - 2; i >= 0; --i) x[i] = d2[i] - c2[i] * x[i + 1];
}
