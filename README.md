# SolarGeoThermal

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)]() [![Visual Studio](https://img.shields.io/badge/IDE-Visual%20Studio%202019/2022-purple)]() [![Platform](https://img.shields.io/badge/Platform-Windows-x64)]()


A C++ simulation project for a solar‑assisted geothermal system with a heat pump, buffer tank, and a solar thermal concentrator.

## Overview
This project simulates a coaxial ground heat exchanger coupled with a heat pump and a buffer tank.  
It adds a solar thermal concentrator that:
- Uses irradiance from `Weather_gansu.csv` (default column `DNI`)
- In heating season: solar heat goes directly to the tank
- In non‑heating season: solar heat is injected to the ground loop for recharge
- Includes temperature limits (outlet <= 100 C; flow increases to avoid overheating)
- Stops when irradiance is too low (pump off)

## Features
- Seasonal operation (heating / off‑season)
- Weather‑driven load model
- Heat pump + tank control
- Ground heat exchanger model
- Hydraulic pressure drop and pump power estimation
- Solar thermal module with irradiance‑based output
- Detailed hourly outputs + summary

## Project Structure
solargeo.sln # Visual Studio solution
enhance/ # C++ source (VS project)
enhance/Weather_gansu.csv # Weather input (with DNI/GHI columns)
enhance/plot_results_solar.py # Plot script (3‑case comparison)
tools/ # Utilities (optional)



## Build (Visual Studio 2022)
1. Open `solargeo.sln`
2. Select `x64 | Release`
3. Build / Rebuild
4. Run: `x64\Release\enhance.exe`

## Key Environment Variables
### Weather / Load
LOAD_WEATHER_CSV=Weather_gansu.csv
LOAD_COLUMN=T_outdoor_C



### Heating Season
HEAT_SEASON_ENABLE=1
HEAT_START_MM=10
HEAT_START_DD=1
HEAT_END_MM=3
HEAT_END_DD=1



### Tank / Heat Pump
TANK_SET_C=40
TANK_DB_K=5
HP_MAX_Q_OUT_KW=350
HP_MIN_SRC_RETURN_C=0



### Solar Thermal (new)
SOLAR_ENABLE=1
SOLAR_WEATHER_CSV=Weather_gansu.csv
SOLAR_IRR_COL=DNI
SOLAR_AREA_M2=400
SOLAR_ETA=0.65
SOLAR_MDOT_NOM_KGPS=1
SOLAR_MDOT_MAX_KGPS=50
SOLAR_T_OUT_MAX_C=100
SOLAR_DT_MAX_C=50
SOLAR_IRR_MIN_WM2=80



## Output Files
- `results.csv`  
  Hourly system output (source/return, COP, HP power, tank, load, pump, etc.)
- `results_solar.csv`  
  Additional hourly output for plotting:  
  solar in/out, ground in/out, load supply/return, flow & velocity, HP on/off, COP, power
- `debug.csv`  
  Heat pump cycle diagnostics (optional)
- `summary.csv`  
  End‑of‑run summary

## Plotting (3‑case comparison)
Use the script below to compare 3 cases on 5 plots:
python enhance/plot_results_solar.py case1/results_solar.csv case2/results_solar.csv case3/results_solar.csv


Outputs are saved to `enhance/plots_solar/`:
- Solar irradiance & solar power
- Solar inlet/outlet temperature
- Ground inlet/outlet temperature
- HP on/off & COP
- Load supply/return + tank temperature

## Notes
- For best results, avoid committing large output CSVs to Git.
- You can add `.gitignore` rules for `results.csv`, `results_solar*.csv`,
