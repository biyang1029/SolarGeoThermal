# Enhance / Ground Heat Exchanger + Heat Pump Simulation (C++ / Visual Studio)

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)]() [![Visual Studio](https://img.shields.io/badge/IDE-Visual%20Studio%202019/2022-purple)]() [![Platform](https://img.shields.io/badge/Platform-Windows-x64)]()

Enhance is a numerical simulation program for a coaxial (annulus + inner) ground heat exchanger coupled with a heat pump and a buffer tank. It supports seasonal operation (heating half-year / off the other half), weather‑driven building loads, variable time step, dual‑side flow reporting, hydraulic pressure drop and pump power estimation, and a PowerShell sweep script.

---

## Features
- Seasonal operation: default heating season Oct/15 – Apr/15; off‑season heat pump stops while the ground recovers naturally
- Variable time step: first 2000 steps = 10 min; then +6 s per step up to 1 h; then fixed 2 h
- Weather‑driven load: based on outdoor temperature and indoor setpoint; load is 0 when `Tout > 26°C` (can be disabled)
- Heat pump + tank: outlet setpoint with deadband ON/OFF; COP hard upper cap removed (keep lower bound and freeze guard)
- Solar thermal concentrator: reads irradiance from weather CSV, heats tank in heating season, recharges ground loop in off-season
- Physical coupling:
  - Inner/annulus convection via Re/Pr/Nu correlations (Dittus–Boelter with transition blending)
  - Series thermal resistances of inner/outer pipe walls, grout, soil
  - Top equivalent insulation zone (length / thermal conductivity configurable), bottom enhancement zone (length configurable)
- Hydraulics / pump: Darcy + minor losses, report total dP and pump power
- Logging: `results.csv` (hourly rows), `debug.csv` (cycle states), end‑run `summary.csv`
- Parameter sweep: `sweep.ps1` runs multiple cases and generates `runs_yyyymmdd_hhmmss/` with `summary_runs.csv` (incl. comparison to measured COP/hours)

---

## Project Layout
```
enhance/
├─ enhance.sln                   # Visual Studio solution
├─ sweep.ps1                     # PowerShell sweep script (run from repo root)
├─ README.md                     # This file
├─ enhance/                      # Source code (VS project)
│  ├─ *.cpp *.h                  # SimulationController / HeatPump / Ground / Logger / etc.
│  ├─ weather_gansu.csv          # Gansu weather sample
│  └─ eigen-3.4.0/               # Dependency (if present)
└─ x64/Release/enhance.exe       # Built executable
```

---

## Build & Run (Visual Studio)
1. Open `enhance.sln`, choose configuration `x64 | Release`
2. Rebuild solution or project
3. Executable: `x64\Release\enhance.exe`

Threads: default 12. You can override via environment variables (choose ONE):
```
$env:OMP_NUM_THREADS='12'
# or
$env:OMP_THREADS='12'
# or percentage of logical CPUs
$env:OMP_THREADS_PCT='50'
```

---

## Single Run – Common Environment Variables (PowerShell)
```
# Weather & load
$env:LOAD_WEATHER_CSV='enhance\weather_gansu.csv'
$env:LOAD_INDOOR_T='26'
$env:LOAD_UA_KW_PER_K='22'
$env:LOAD_BASE_KW='0'
$env:LOAD_CUTOFF_ENABLE='1'
$env:LOAD_HEAT_CUTOFF_C='26'

# Heating season (Oct/15 – Apr/15)
$env:HEAT_SEASON_ENABLE='1'
$env:HEAT_START_MM='10'; $env:HEAT_START_DD='15'
$env:HEAT_END_MM='4';    $env:HEAT_END_DD='15'

# Heat pump / tank
$env:TANK_VOL_M3='8'
$env:TANK_SET_C='40'
$env:TANK_DB_K='2'
$env:HP_MAX_Q_OUT_KW='350'
$env:HP_EFF_CARNOT='0.65'
$env:HP_EVAP_APPROACH_K='3'
$env:HP_COND_APPROACH_K='3'
$env:HP_MIN_SRC_RETURN_C='10'
$env:HP_MOD_RANGE_C='1'
$env:HP_MAX_SRC_DT_PER_H='0'  # 0 = disable source-side dT cap

# Geometry / materials / enhancement / top insulation
$env:D_OUTER_M='0.2'; $env:D_INNER_M='0.1'
$env:BORE_D_M='0.22'; $env:PIPE_THICK_M='0.008'
$env:PIPE_K_INNER='0.4'    # PE-RT II
$env:PIPE_K_OUTER='4.0'    # J55
$env:INSUL_TOP_LEN_M='100'; $env:INSUL_K_INNER='0.1'
$env:EHEP_BOTTOM_LEN_M='100'
$env:SOIL_K='2.5'; $env:SOIL_RHO='2600'; $env:SOIL_CP='1600'; $env:GROUT_K='4'

# Flows (two sides)
$env:FLOW_SRC_KGPS='10'      # source side (kg/s), used in hydraulics & ground
$env:FLOW_LOAD_KGPS='22.2'   # load side (kg/s), for reporting

# Run
./x64/Release/enhance.exe
```

---

## Sweep Script – `sweep.ps1`
Run from repository root:
```
powershell -NoProfile -ExecutionPolicy Bypass -File .\sweep.ps1
```
What it does:
- Locates `x64\Release\enhance.exe`
- Applies baseline env and sweeps across lists declared at the top (e.g., `$hpList`, `$uaList`, `$pipeKInnerList`)
- Archives each run into `runs_YYYYMMDD_HHMMSS/` and writes `summary_runs.csv`

`summary_runs.csv` columns:
- `hp_max_kW, UA_kWperK, flow_src_kgps, flow_load_kgps, set_load_out_C, tank_vol_m3,
   Q_load_kWh, Q_src_kWh, P_el_kWh, P_pump_kWh, dP_kPa_avg,
   avg_Q_load_kW, avg_Q_src_kW, COP_annual, COP_measured, COP_delta,
   HP_on_hours, HP_on_measured, HP_on_delta`
(Default comparison: `COP_measured = 4.3`, `HP_on_measured = 2890`.)

---

## Output Files
- `results.csv` (hourly rows):
  - `step,date,time,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,model,
     T_tank_C,HP_on,Q_space_req_kW,Q_dhw_req_kW,Q_space_served_kW,Q_dhw_served_kW,Q_unmet_kW,
     flow_src_kgps,flow_load_kgps,dP_kPa,P_pump_kW`
- `results_solar.csv` (hourly rows): solar/ground/load temperatures, flow rates, HP on/off, COP, and power terms (for plotting)
- `debug.csv`: cycle state and CoolProp diagnostics
- `summary.csv`: end‑of‑run energy/heat/peaks/economics
- `runs_*/results_*.csv` / `debug_*.csv`: archived sweep results
- `runs_*/summary_runs.csv`: sweep summary (including measured comparisons)

---

## Plotting (3-case comparison)
The script `enhance/plot_results_solar.py` draws five comparison figures from three cases.
Edit `CASE_FILES` at the top or pass three files on the command line:
```
python plot_results_solar.py path\to\case1.csv path\to\case2.csv path\to\case3.csv
```
Figures are saved to `enhance/plots_solar/`.

---

## Notes (progress & performance)
- The progress bar updates on integer percent and tries to refresh on a single line. Some hosts (e.g., Visual Studio “Output” window) convert carriage return to newline and you may see stacked lines. For the single‑line flicker, run the executable in a normal PowerShell/terminal.
- Default threads = 12; you can override via `OMP_NUM_THREADS` / `OMP_THREADS` / `OMP_THREADS_PCT`.

---

## Git / GitHub Tips
- Recommended ignore entries: `.vs/`, `x64/`, `Debug/`, `Release/`, `*.pdb`, `*.obj`, `enhance/results.csv`, `enhance/debug.csv`, `runs_*/`, `opt_runs*/`, `enhance/summary.csv`
- Visual Studio flow: Fetch/Pull → resolve conflicts & Commit → Push
