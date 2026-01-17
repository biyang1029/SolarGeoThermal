// ========================= SimulationController.h =========================
#pragma once
#include <string>
#include "DataConfig.h"
#include "GroundModule.h"
#include "HeatPumpModule.h"
#include "TankModule.h"
#include "SolarModule.h"
#include "Logger.h"
#include "LoadModel.h"


class SimulationController {
public:
	explicit SimulationController(const DataConfig& cfg) : cfg_(cfg) {}
	bool run(const std::string& csv_path = "results.csv");


private:
    DataConfig cfg_;
    GroundModule ground_;
    HeatPumpModule hp_;
    TankModule    tank_;
    SolarModule   solar_;
    Logger logger_;
    LoadModel load_;
    bool hp_on_ = false;
};
