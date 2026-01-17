#define _CRT_SECURE_NO_WARNINGS 1
#include "HeatPumpModule.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <windows.h>
#endif




static inline double clampd(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

// ... CoolPropDyn 结构体保持不变 ...

#ifdef _WIN32
struct CoolPropDyn {
    HMODULE h = nullptr;
    typedef double(__cdecl* PropsSI_t)(const char*, const char*, double, const char*, double, const char*);
    PropsSI_t PropsSI = nullptr;
    bool load_once = false;
    bool try_load(const std::string& hintDir = std::string()) {
        if (load_once) return h && PropsSI;
        load_once = true;
        auto try_path = [&](const std::string& path) {
            h = LoadLibraryA(path.c_str());
            if (h) { PropsSI = (PropsSI_t)GetProcAddress(h, "PropsSI"); }
            return h && PropsSI;
            };
        if (!hintDir.empty()) {
            std::string p = hintDir;
            char back = p.empty() ? '\0' : p.back();
            if (back != '\\' && back != '/') p += "\\";
            if (try_path(p + "CoolProp.dll")) return true;
        }
        char* env = nullptr; size_t len = 0;
        if (_dupenv_s(&env, &len, "COOLPROP_DIR") == 0 && env) {
            std::string base(env); free(env);
            if (try_path(base + "\\bin\\CoolProp.dll")) return true;
            if (try_path(base + "\\CoolProp.dll")) return true;
        }
        if (try_path("CoolProp.dll")) return true;
        if (_dupenv_s(&env, &len, "USERPROFILE") == 0 && env) {
            free(env);
            std::string guess = "D:/yangb/soft/anaconda/Lib/site-packages/CoolProp/CoolProp.dll";
            try_path(guess);
        }
        return h && PropsSI;
    }
};
static CoolPropDyn g_coolprop;
#endif
#ifdef _WIN32
static void calibrate_design_point(
    const DataConfig& cfg,
    CoolPropDyn& cp,
    double& q_cond_design_Jperkg,
    double& rho_suc_nom,
    bool& design_ready)
{
    if (design_ready) return;

    const std::string& fluid = cfg.hp.fluid;

    // 1) 用“名义工况”构造蒸发/冷凝饱和温度（与 step 里一致的建模逻辑）
    double Tevap_sat_K = (cfg.hp.design_T_source_out_C - cfg.hp.evap_approach_K) + 273.15;
    Tevap_sat_K = std::max(cfg.hp.evap_T_min_C + 273.15, Tevap_sat_K);

    double Tcond_sat_K = (cfg.hp.target_T_load_out_C + cfg.hp.cond_approach_K) + 273.15;

    // 2) 安全夹紧（跟你 step 里同样的 trip/crit 逻辑）
    double Tcrit = 1e9, Ttrip = 50.0;
    try { Tcrit = cp.PropsSI("Tcrit", "", 0.0, "", 0.0, fluid.c_str()); }
    catch (...) {}
    try { Ttrip = cp.PropsSI("T_triple", "", 0.0, "", 0.0, fluid.c_str()); }
    catch (...) {}
    const double Kmin = 2.0;

    auto clampd = [](double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); };

    Tcond_sat_K = clampd(Tcond_sat_K, Ttrip + Kmin, Tcrit - 5.0);

    if (Tcond_sat_K - Tevap_sat_K < cfg.hp.min_temp_lift_K) {
        Tcond_sat_K = Tevap_sat_K + cfg.hp.min_temp_lift_K;
        Tcond_sat_K = clampd(Tcond_sat_K, Ttrip + Kmin, Tcrit - 5.0);
    }

    Tevap_sat_K = std::min(Tevap_sat_K, Tcond_sat_K - cfg.hp.min_temp_lift_K);
    Tevap_sat_K = std::max(Tevap_sat_K, Ttrip + Kmin);

    // 3) 求蒸发/冷凝压力
    double P_cond_Pa = cp.PropsSI("P", "T", Tcond_sat_K, "Q", 0.0, fluid.c_str());
    double P_evap_Pa = cp.PropsSI("P", "T", Tevap_sat_K, "Q", 1.0, fluid.c_str());

    // 4) 状态点 1：过热
    double T_superheat = Tevap_sat_K + std::max(Kmin, cfg.hp.superheat_K);
    double h1 = cp.PropsSI("H", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());
    double s1 = cp.PropsSI("S", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());

    // 5) 状态点 2：压缩后
    double h2s = cp.PropsSI("H", "P", P_cond_Pa, "S", s1, fluid.c_str());
    double h2 = h1 + (h2s - h1) / std::max(0.1, cfg.hp.eta_isentropic);

    // 6) 状态点 3：过冷（按你希望“比水箱温度高3℃”）
    // 名义水箱回水 = design_T_load_in_C
    double T3_K_base = (cfg.hp.design_T_load_in_C + 3.0) + 273.15; // +3℃
    double T_subcool_min_K = std::max(3.0, cfg.hp.subcool_K);
    double T3_K_max = Tcond_sat_K - T_subcool_min_K;
    double T3_K = clampd(T3_K_base, Ttrip + Kmin, T3_K_max);

    double h3 = cp.PropsSI("H", "T", T3_K, "P", P_cond_Pa, fluid.c_str());

    // 7) 计算名义单位供热量 & 名义吸气密度
    q_cond_design_Jperkg = (h2 - h3); // J/kg
    rho_suc_nom = cp.PropsSI("D", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());

    // 合理性保护
    if (!(std::isfinite(q_cond_design_Jperkg) && q_cond_design_Jperkg > 1e3)) q_cond_design_Jperkg = 0.0;
    if (!(std::isfinite(rho_suc_nom) && rho_suc_nom > 1e-6)) rho_suc_nom = 0.0;

    design_ready = (q_cond_design_Jperkg > 0.0 && rho_suc_nom > 0.0);
}
#endif

double HeatPumpModule::step(double T_source_out_C,
    double T_load_in_C,
    double m_dot_load_kgps,
    double dt_s,
    double& COP,
    double& Q_out_kW,
    double& P_el_kW) {
    const std::string& fluid = cfg_.hp.fluid;
    COP = 0.0; Q_out_kW = 0.0; P_el_kW = 0.0;
    last_debug_ = HPDebug{}; last_debug_.fluid = fluid;
    static bool warned_failure = false;

    // 1. 蒸发侧压力和温度确定 (与原方案相同，但 Tevap_sat_K 需要被冷凝侧温度约束)
    double Tevap_sat_K = (T_source_out_C - cfg_.hp.evap_approach_K) + 273.15;
    Tevap_sat_K = std::max(cfg_.hp.evap_T_min_C + 273.15, Tevap_sat_K);

    //double T_superheat_K = Tevap_sat_K + std::max(3.0, cfg_.hp.superheat_K);
    double m_dot_src = std::max(1e-9, m_dot_src_kgps_);
    double cp_src = cfg_.fluid.cp;
    double cp_load = cfg_.fluid.cp; // 简化：负荷侧用同一物性

    bool used_cp = false;

#ifdef _WIN32
    if (cfg_.hp.use_coolprop && g_coolprop.try_load()) {
        try {
            // 缓存临界温度和三相点温度
            calibrate_design_point(cfg_, g_coolprop, q_cond_design_Jperkg_, rho_suc_nom_kgpm3_, design_ready_);
            static std::string lastFluid; static double Tcrit_cached = 0.0, Ttrip_cached = 0.0; static bool fluidInfoReady = false;
            if (fluid != lastFluid || !fluidInfoReady) {
                try { Tcrit_cached = g_coolprop.PropsSI("Tcrit", "", 0.0, "", 0.0, fluid.c_str()); }
                catch (...) { Tcrit_cached = 1e9; }
                try { Ttrip_cached = g_coolprop.PropsSI("T_triple", "", 0.0, "", 0.0, fluid.c_str()); }
                catch (...) { Ttrip_cached = 50.0; }
                lastFluid = fluid; fluidInfoReady = true;
            }
            double Tcrit = (Tcrit_cached > 0.0) ? Tcrit_cached : 1e9;
            double Ttrip = (Ttrip_cached > 0.0) ? Ttrip_cached : 50.0;
            const double Kmin = 2.0; // 最小安全温差 (Ttrip + Kmin)

            // --- [Start of Modified CoolProp Logic: Load-Driven Condenser] ---

            // 2. 冷凝饱和温度由负载侧决定 (T_load_set_C)
            // 使用目标供水温度来确定冷凝侧所需的饱和温度
            double Tcond_sat_K_base = (cfg_.hp.target_T_load_out_C + cfg_.hp.cond_approach_K) + 273.15;

            // 3. 边界保护: 确保 Tcond 在安全范围内
            double Tcond_sat_K = clampd(Tcond_sat_K_base, Ttrip + Kmin, Tcrit - 5.0);

            // 4. 最小温升保护: 确保 Tcond_sat_K - Tevap_sat_K >= min_temp_lift_K
            if (Tcond_sat_K - Tevap_sat_K < cfg_.hp.min_temp_lift_K) {
                Tcond_sat_K = Tevap_sat_K + cfg_.hp.min_temp_lift_K;
                Tcond_sat_K = clampd(Tcond_sat_K, Ttrip + Kmin, Tcrit - 5.0); // 再次钳位，防止溢出 Tcrit
            }

            // 5. 最终确定 P_cond
            double P_cond_Pa = g_coolprop.PropsSI("P", "T", Tcond_sat_K, "Q", 0.0, fluid.c_str());

            // 6. 最终确定 P_evap (可能需要基于 Tcond_sat_K 再次约束 Tevap_sat_K)
            Tevap_sat_K = std::min(Tevap_sat_K, Tcond_sat_K - cfg_.hp.min_temp_lift_K);
            Tevap_sat_K = std::max(Tevap_sat_K, Ttrip + Kmin); // 确保 Tevap 不低于 Ttrip
            double P_evap_Pa = g_coolprop.PropsSI("P", "T", Tevap_sat_K, "Q", 1.0, fluid.c_str());

            // --- [End of Modified CoolProp Logic] ---

            last_debug_.P_evap_kPa = P_evap_Pa / 1000.0;
            last_debug_.P_cond_kPa = P_cond_Pa / 1000.0;
            last_debug_.T_evap_sat_K = Tevap_sat_K;
            last_debug_.T_cond_sat_K = Tcond_sat_K;

            // 状态点 1: 蒸发器出口 (过热蒸汽)
            double T_superheat = Tevap_sat_K + std::max(Kmin, cfg_.hp.superheat_K);
            double h1 = g_coolprop.PropsSI("H", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());
            double s1 = g_coolprop.PropsSI("S", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());

            // 状态点 2s/2: 压缩机出口 (实际压缩)
            double h2s = g_coolprop.PropsSI("H", "P", P_cond_Pa, "S", s1, fluid.c_str());
            double h2 = h1 + (h2s - h1) / std::max(0.1, cfg_.hp.eta_isentropic);

            // 状态点 3: 冷凝器出口 (过冷液体, 使用负载回水温度 T_load_in_C)
            // T3_K = clamp(T_load_in_C + 273.15, T_trip+2K, T_cond_sat_K - max(3K, subcool_K))
            double T_subcool_min_K = std::max(3.0, cfg_.hp.subcool_K);
            double T3_K_max = Tcond_sat_K - T_subcool_min_K;
            double T3_K_base = T_load_in_C + 276.15; // 使用负载回水温度 T_load_in_C 作为基准
            double T3_K = clampd(T3_K_base, Ttrip + Kmin, T3_K_max);

            // 查找 h3 (随 T_load_in_C 变化)
            double h3 = g_coolprop.PropsSI("H", "T", T3_K, "P", P_cond_Pa, fluid.c_str());
            double h4 = h3; // 状态点 4: 膨胀阀出口 (等焓)

            last_debug_.used_coolprop = true;
            last_debug_.h1 = h1; last_debug_.h2s = h2s; last_debug_.h2 = h2; last_debug_.h3 = h3;

            // 单位制冷剂供热与压缩功 (J/kg)
            double q_cond = (h2 - h3);
            double w_comp = (h2 - h1);

            double rho_suc = g_coolprop.PropsSI("D", "T", T_superheat, "P", P_evap_Pa, fluid.c_str());

            // ===== compressor capacity from Q_nom (auto) =====
            double m_ref_nom = 0.0;
            if (design_ready_ && q_cond_design_Jperkg_ > 1e-9) {
                // m_ref_nom is derived from Q_nom at design point
                m_ref_nom = (cfg_.hp.Q_nom_kW * 1000.0) / q_cond_design_Jperkg_; // kg/s
            }

            // fallback safety: if not ready, keep old parameter (optional)
            if (m_ref_nom <= 0.0) {
                m_ref_nom = cfg_.hp.m_ref_nom_kgps; // 若你想彻底不手填，可把这一行改成一个小常数
            }

            // density correction
            double rho_nom = (design_ready_ ? rho_suc_nom_kgpm3_ : cfg_.hp.rho_suc_nom_kgpm3);
            double m_ref_max = m_ref_nom;
            if (rho_nom > 1e-6) m_ref_max *= (rho_suc / rho_nom);




            // ===== 水侧温度目标 → 所需热量 =====
            const double T_load_target_C = cfg_.hp.target_T_load_out_C; // 60℃
            double Q_need_kW = 0.0;
            if (m_dot_load_kgps > 1e-9) {
                double dT_need = T_load_target_C - T_load_in_C;
                if (dT_need > 0.0)
                    Q_need_kW = (m_dot_load_kgps * cp_load * dT_need) / 1000.0;
            }

            // ===== 由 Q_need 反算制冷剂流量 =====
            double m_ref_need = (q_cond > 1e-9) ? (Q_need_kW * 1000.0 / q_cond) : 0.0;

            // 实际制冷剂流量：定频上限裁剪
            double m_ref_use = std::min(m_ref_need, m_ref_max);

            // ===== 实际输出 =====
            double Q_cond_W = m_ref_use * q_cond;
            double P_el_W = m_ref_use * w_comp / std::max(0.1, cfg_.hp.eta_motor);

            Q_out_kW = Q_cond_W / 1000.0;
            P_el_kW = P_el_W / 1000.0;
            COP = (P_el_W > 1e-6) ? (Q_cond_W / P_el_W) : 0.0;

            used_cp = true;
        }
        catch (...) {
            if (!warned_failure) { std::cerr << "[Warning] CoolProp call threw. Fallback model used.\n"; warned_failure = true; }
        }
    }
#endif

    // ... Fallback 模型和边界条件检查逻辑保持不变 ...

    if (!used_cp) {
        // Fallback 模式（基于 T_load_in_C）
        double T_cond_K = (T_load_in_C + cfg_.hp.cond_approach_K) + 273.15;
        double T_evap_K = Tevap_sat_K + std::max(3.0, cfg_.hp.superheat_K);
        const double dT_min = 5.0;
        if (T_cond_K - T_evap_K < dT_min) T_evap_K = T_cond_K - dT_min;
        double COP_carnot = (T_cond_K > T_evap_K + 1e-6) ? (T_cond_K / (T_cond_K - T_evap_K)) : 1.0;
        double eff = cfg_.hp.eff_carnot * (0.9 * cfg_.hp.eta_isentropic + 0.1);
        COP = (std::max)(1.5, eff * COP_carnot);
        double Q_need_kW = 0.0;
        if (m_dot_load_kgps > 1e-9) {
            double dT = cfg_.hp.target_T_load_out_C - T_load_in_C;
            if (dT > 0.0)
                Q_need_kW = m_dot_load_kgps * cp_load * dT / 1000.0;
        }
        Q_out_kW = std::min(Q_need_kW, cfg_.hp.max_Q_out_kW);
        P_el_kW = (Q_out_kW > 1e-9 ? Q_out_kW / COP : 0.0);
    }

    if (m_dot_load_kgps > 1e-9) {
        const double T_load_target_C = cfg_.hp.target_T_load_out_C;
        if (T_load_in_C >= T_load_target_C - 1e-6) {
            Q_out_kW = 0.0;
            P_el_kW = 0.0;
            COP = 0.0;
            return T_source_out_C;
        }
    }


    {
        double guard_lo = cfg_.hp.min_source_return_C;
        double guard_span = 3.0;
        if (std::isfinite(guard_lo) && guard_span > 1e-9) {
            double guard_hi = guard_lo + guard_span;
            double scale_guard = 1.0;
            if (T_source_out_C <= guard_lo) scale_guard = 0.0;
            else if (T_source_out_C < guard_hi) scale_guard = (T_source_out_C - guard_lo) / (guard_hi - guard_lo);

            scale_guard = clampd(scale_guard, 0.0, 1.0);
            Q_out_kW *= scale_guard;
            P_el_kW *= scale_guard;
            // COP 理论上不变，但数值上保持一致
            if (P_el_kW > 1e-9) COP = Q_out_kW / P_el_kW;
            else COP = 0.0;
        }
    }

    /*
    double Q_out_base_kW = Q_out_kW;
    double Q_out_demand_kw = Q_target;
    Q_out_kW = std::max(0.0, std::min(Q_out_base_kW, Q_out_demand_kW));
    if (Q_out_base_kW > 1e-9 && Q_out_kW < Q_out_base_kW) {
        double ratio = Q_out_kW / Q_out_base_kW;
        P_el_kW *= ratio;
        COP = (P_el_kW > 1e-6) ? (Q_out_kW / P_el_kW) : COP;
    }
    if (Q_out_kW < 1e-9) { P_el_kW = 0.0; COP = 0.0; return T_source_out_C; }
    */
    // 负荷侧出水温
    double T_load_out_C = T_load_in_C;
    if (m_dot_load_kgps > 1e-9) {
        T_load_out_C = T_load_in_C + (Q_out_kW * 1000.0) / (m_dot_load_kgps * cp_load);
    }

    // 源侧回水
    double Q_geo_kW = Q_out_kW - P_el_kW;
    double dT_need = (Q_geo_kW * 1000.0) / (m_dot_src * cp_src);
    double T_return_C = T_source_out_C - dT_need;
    if (cfg_.hp.min_source_return_C > -1e9 && T_return_C < cfg_.hp.min_source_return_C) T_return_C = cfg_.hp.min_source_return_C;

    return T_return_C;

}