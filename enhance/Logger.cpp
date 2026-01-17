// ========================= Logger.cpp =========================
#include "Logger.h"
#include <iomanip>


bool Logger::open(const std::string& path) {
    // Locally close any existing streams without calling close()
    if (this->ofs_.is_open()) this->ofs_.close();
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
    this->ofs_.open(path, std::ios::out | std::ios::trunc);
    if (!this->ofs_) return false;
    this->ofs_ << "step,date,time,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,T_outdoor_C,T_load_supply_C,T_load_return_C,T_hp_load_out_C,model,T_tank_C,HP_on,Q_space_req_kW,Q_dhw_req_kW,Q_space_served_kW,Q_dhw_served_kW,Q_unmet_kW,flow_src_kgps,flow_load_kgps,dP_kPa,P_pump_kW,solar_on,solar_mode,solar_irr_Wm2,solar_mdot_kgps,solar_T_in_C,solar_T_out_C,solar_Q_kW,solar_Q_raw_kW\n";
    return true;
}

bool Logger::openDebug(const std::string& path) {
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
    this->ofs_dbg_.open(path, std::ios::out | std::ios::trunc);
    if (!this->ofs_dbg_) return false;
    this->ofs_dbg_ << "hour,fluid,used_coolprop,T_evap_sat_K,T_cond_sat_K,P_evap_kPa,P_cond_kPa,h1,h2s,h2,h3,T_source_in_C,T_geo_in_C,Q_geo_kW,flow_src_kgps,flow_load_kgps,T_outdoor_C\n";
    return true;
}


void Logger::writeHour(int hour,
	double T_source_out_C,
	double T_return_C,
	double COP,
	double Q_out_kW,
	double P_el_kW,
	double Q_geo_kW,
	double T_outdoor_C,
    double T_load_supply_C,
    double T_load_return_C,
    double T_hp_load_out_C,
	const std::string& model,
	const std::string& date_ymd,
	const std::string& time_hm,
	double T_tank_C,
	int    HP_on,
	double Q_space_req_kW,
	double Q_dhw_req_kW,
	double Q_space_served_kW,
	double Q_dhw_served_kW,
	double Q_unmet_kW,
	double flow_src_kgps,
	double flow_load_kgps,
	double dP_kPa,
	double P_pump_kW,
    int    solar_on,
    int    solar_mode,
    double solar_irr_Wm2,
    double solar_mdot_kgps,
    double solar_T_in_C,
    double solar_T_out_C,
    double solar_Q_kW,
    double solar_Q_raw_kW) {
	if (!this->ofs_) return;
	this->ofs_ << hour << ',' << date_ymd << ',' << time_hm << ','
		<< std::fixed << std::setprecision(4)
		<< T_source_out_C << ','
		<< T_return_C << ','
		<< COP << ','
		<< Q_out_kW << ','
		<< P_el_kW << ','
		<< Q_geo_kW << ','
        << T_outdoor_C << ','
        << T_load_supply_C << ','
        << T_load_return_C << ','
        << T_hp_load_out_C << ','
        << model << ','
		<< T_tank_C << ','
		<< HP_on << ','
		<< Q_space_req_kW << ','
		<< Q_dhw_req_kW << ','
		<< Q_space_served_kW << ','
		<< Q_dhw_served_kW << ','
		<< Q_unmet_kW << ','
		<< flow_src_kgps << ','
		<< flow_load_kgps << ','
        << dP_kPa << ','
        << P_pump_kW << ','
        << solar_on << ','
        << solar_mode << ','
        << solar_irr_Wm2 << ','
        << solar_mdot_kgps << ','
        << solar_T_in_C << ','
        << solar_T_out_C << ','
        << solar_Q_kW << ','
        << solar_Q_raw_kW << "\n";
}

void Logger::writeDebugHour(int hour,
    const std::string& fluid,
    bool used_coolprop,
    double T_evap_sat_K,
    double T_cond_sat_K,
    double P_evap_kPa,
    double P_cond_kPa,
    double h1,
    double h2s,
    double h2,
    double h3,
    double T_source_in_C,
    double T_geo_in_C,
    double Q_geo_kW,
    double flow_src_kgps,
    double flow_load_kgps,
    double T_outdoor_C) {
    if (!this->ofs_dbg_) return;
    this->ofs_dbg_ << hour << ',' << fluid << ',' << (used_coolprop?1:0) << ','
        << std::fixed << std::setprecision(6)
        << T_evap_sat_K << ',' << T_cond_sat_K << ','
        << P_evap_kPa << ',' << P_cond_kPa << ','
        << h1 << ',' << h2s << ',' << h2 << ',' << h3 << ','
        << std::setprecision(4)
        << T_source_in_C << ',' << T_geo_in_C << ',' << Q_geo_kW << ','
        << flow_src_kgps << ',' << flow_load_kgps << ','
        << T_outdoor_C << "\n";
}


void Logger::close() {
    if (this->ofs_.is_open()) this->ofs_.close();
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
}


