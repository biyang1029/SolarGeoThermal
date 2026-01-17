#!/usr/bin/env python3
import csv
import math
import os
import sys

import matplotlib.pyplot as plt


def safe_float(v):
    try:
        if v is None:
            return math.nan
        s = str(v).strip()
        if s == "":
            return math.nan
        return float(s)
    except Exception:
        return math.nan


CASE_FILES = [
    r"D:\\yangb\\material\\OneDrive\\Private\\1 postdoc\\1.artilcle\\21.geosolar\\Solargeo\\enhance\\case1\\results_solar.csv",
    r"D:\\yangb\\material\\OneDrive\\Private\\1 postdoc\\1.artilcle\\21.geosolar\\Solargeo\\enhance\\case2\\results_solar.csv",
    r"D:\\yangb\\material\\OneDrive\\Private\\1 postdoc\\1.artilcle\\21.geosolar\\Solargeo\\enhance\\case3\\results_solar.csv",
]

CASE_LABELS = [
    "case1",
    "case2",
    "case3",
]


def read_results(path):
    data = {}
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return data
        for name in reader.fieldnames:
            data[name] = []
        for row in reader:
            for name in reader.fieldnames:
                data[name].append(row.get(name, ""))
    return data


def get_series(data, name):
    if name not in data:
        return None
    return [safe_float(v) for v in data[name]]


def get_hour_axis(data):
    if "hour" in data:
        return [safe_float(v) for v in data["hour"]]
    if "step" in data:
        return [safe_float(v) for v in data["step"]]
    # Fallback to index
    if data:
        n = len(next(iter(data.values())))
        return list(range(n))
    return []


def ensure_plot_dir():
    out_dir = "plots_solar"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    return out_dir


def load_cases(paths, labels):
    cases = []
    for i, path in enumerate(paths):
        if not os.path.exists(path):
            print("Missing file:", path)
            continue
        data = read_results(path)
        if not data:
            print("No data in:", path)
            continue
        hours = get_hour_axis(data)
        label = labels[i] if labels and i < len(labels) else os.path.splitext(os.path.basename(path))[0]
        cases.append({"path": path, "label": label, "data": data, "hours": hours})
    return cases


def slice_series(case, name, n):
    series = get_series(case["data"], name)
    if series is None:
        return None
    return series[:n]


def plot_series(ax, cases, name, ylabel, n, step=False):
    plotted = 0
    for case in cases:
        y = slice_series(case, name, n)
        if y is None:
            continue
        x = case["hours"][:n]
        if step:
            ax.step(x, y, where="mid", label=case["label"])
        else:
            ax.plot(x, y, label=case["label"])
        plotted += 1
    if plotted:
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.grid(True, alpha=0.25)
    return plotted > 0


def main():
    paths = CASE_FILES[:]
    labels = CASE_LABELS[:]
    if len(sys.argv) >= 4:
        paths = sys.argv[1:4]
        labels = [os.path.splitext(os.path.basename(p))[0] for p in paths]

    cases = load_cases(paths, labels)
    if not cases:
        print("No valid case files loaded. Update CASE_FILES or pass 3 files as args.")
        return 2

    min_len = min(len(c["hours"]) for c in cases if c["hours"])
    if min_len <= 0:
        print("No time axis found.")
        return 3

    out_dir = ensure_plot_dir()

    # 1) Solar irradiance and solar power
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ok1 = plot_series(axes[0], cases, "solar_irr_Wm2", "Irradiance (W/m2)", min_len, step=False)
    ok2 = plot_series(axes[1], cases, "solar_Q_kW", "Solar Q (kW)", min_len, step=False)
    axes[1].set_xlabel("Hour")
    if ok1 or ok2:
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig1_solar_irr_Q.png"), dpi=160)
    plt.close(fig)

    # 2) Solar inlet/outlet temperature
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ok1 = plot_series(axes[0], cases, "solar_T_in_C", "Solar T_in (C)", min_len, step=False)
    ok2 = plot_series(axes[1], cases, "solar_T_out_C", "Solar T_out (C)", min_len, step=False)
    axes[1].set_xlabel("Hour")
    if ok1 or ok2:
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig2_solar_T.png"), dpi=160)
    plt.close(fig)

    # 3) Ground inlet/outlet temperature
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ok1 = plot_series(axes[0], cases, "T_ground_in_C", "Ground T_in (C)", min_len, step=False)
    ok2 = plot_series(axes[1], cases, "T_ground_out_C", "Ground T_out (C)", min_len, step=False)
    axes[1].set_xlabel("Hour")
    if ok1 or ok2:
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig3_ground_T.png"), dpi=160)
    plt.close(fig)

    # 4) HP on/off and COP
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ok1 = plot_series(axes[0], cases, "COP", "COP", min_len, step=False)
    ok2 = plot_series(axes[1], cases, "HP_on", "HP_on (0/1)", min_len, step=True)
    axes[1].set_xlabel("Hour")
    if ok1 or ok2:
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig4_hp_cop.png"), dpi=160)
    plt.close(fig)

    # 5) Load supply/return and tank temperature
    fig, axes = plt.subplots(3, 1, figsize=(10, 7), sharex=True)
    ok1 = plot_series(axes[0], cases, "T_load_supply_C", "Load supply T (C)", min_len, step=False)
    ok2 = plot_series(axes[1], cases, "T_load_return_C", "Load return T (C)", min_len, step=False)
    ok3 = plot_series(axes[2], cases, "T_tank_C", "Tank T (C)", min_len, step=False)
    axes[2].set_xlabel("Hour")
    if ok1 or ok2 or ok3:
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig5_load_tank_T.png"), dpi=160)
    plt.close(fig)

    print("Saved plots to:", out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
