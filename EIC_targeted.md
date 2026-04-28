# Borrelia-lipidomics_2026
This repository contains python code used to plot extracted ion chromatogram

"""
Plot MassHunter EIC panels: BMP examples

# BMP: 4 rows x 2 columns (ACCESSIBLE VERSION)
# Left:  Bb vs lipid extraction control
# Right: CM vs lipid extraction control
# - Big fonts for accessibility
# - CM y-axis forced to match Bb y-axis for each lipid
# - Lipid title is BOLD and uses LIPID_TITLE_FS
# - Big figure size

import re
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Replicate IDs
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

# Colors
B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 3.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# mpl.rcParams["svg.fonttype"] = "path"

# -----------------------------
# ACCESSIBILITY FONT SIZES
# -----------------------------
LABEL_FS  = 30
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# Lipid title (subplot title) ONLY
LIPID_TITLE_FS = 50


# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# Big figure size
PANEL_W = 11.0
PANEL_H = 7.5

# -----------------------------
# Files / RT windows
# -----------------------------
PLOTS = [
    {"bb_txt": "BMP_16_0_18_2.txt",     "cm_txt": "BMP_16_0_18_2_CM.txt",     "rt_min": 5.278, "rt_max": 5.492,
     "title": "(16:0/18:2)"},
    {"bb_txt": "BMP_18_1_20_0.txt",     "cm_txt": "BMP_18_1_20_0_CM.txt",     "rt_min": 6.308, "rt_max": 6.618,
     "title": "(18:1/20:0)"},
    {"bb_txt": "BMP_18_1_20_1.txt",     "cm_txt": "BMP_18_1_20_1_CM.txt",     "rt_min": 7.0,   "rt_max": 7.4,
     "title": "(18:1/20:1)"},
    {"bb_txt": "BMP_18_1_22_1_new.txt", "cm_txt": "BMP_18_1_22_1_CM.txt",     "rt_min": 5.6,   "rt_max": 6.2,
     "title": "(18:1/22:1)"},
]

# -----------------------------
# Parsing / stats
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, letter=None, show_letter=True):
    ax.set_xlim(rt_min, rt_max)

    # Lipid title: BOLD + custom size
    ax.set_title(title, fontsize=LIPID_TITLE_FS, fontweight="bold", pad=14)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if show_letter and letter:
        ax.text(-0.14, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_pair(axL, axR, bb_txt, cm_txt, rt_min, rt_max, title, letter):
    blocks_bb = parse_chrom_txt(bb_txt)
    blocks_cm = parse_chrom_txt(cm_txt)

    # Left: Bb vs control
    rt_b, mean_b, sem_b = extract_group(blocks_bb, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks_bb, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    axL.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    axL.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)
    axL.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="Lipid extraction control")
    axL.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(axL, rt_min, rt_max, title, letter=letter, show_letter=True)

    # Right: CM vs control
    rt_m, mean_m, sem_m = extract_group(blocks_cm, CM_REPS, rt_min, rt_max)
    rt_c2, mean_c2, sem_c2 = extract_group(blocks_cm, C_REPS, rt_min, rt_max)
    mean_c2 = np.interp(rt_m, rt_c2, mean_c2)
    sem_c2  = np.interp(rt_m, rt_c2, sem_c2)

    axR.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="Complete medium")
    axR.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)
    axR.plot(rt_m, mean_c2, color=C_LINE, lw=LINEWIDTH, label="Lipid extraction control")
    axR.fill_between(rt_m, mean_c2-sem_c2, mean_c2+sem_c2, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(axR, rt_min, rt_max, title, letter=None, show_letter=False)

    # Force same y-scale (CM matches Bb)
    axR.set_ylim(axL.get_ylim())


# -----------------------------
# Run
# -----------------------------
n = len(PLOTS)
fig, axes = plt.subplots(
    nrows=n, ncols=2,
    figsize=(PANEL_W * 2, PANEL_H * n),
    dpi=300,
    constrained_layout=True
)
axes = np.atleast_2d(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    plot_pair(
        axes[i, 0], axes[i, 1],
        cfg["bb_txt"], cfg["cm_txt"],
        cfg["rt_min"], cfg["rt_max"],
        cfg["title"], letters[i]
    )

# Legends on top row
hL, lL = axes[0, 0].get_legend_handles_labels()
hR, lR = axes[0, 1].get_legend_handles_labels()

axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.42, 0.98), borderaxespad=0)
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "BMP_Bb_vs_CM_accessible_sameY"
fig.savefig("BMP_Bb_vs_CM_accessible_sameY.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("BMP_Bb_vs_CM_accessible_sameY.svg", bbox_inches="tight", facecolor="white")
plt.show()







#Plotting Phosphatidylethanolamine EIC#
import re
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Replicate IDs
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

# Colors
B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 3.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# mpl.rcParams["svg.fonttype"] = "path"

# -----------------------------
# ACCESSIBILITY FONT SIZES
# -----------------------------
LABEL_FS  = 30
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# Lipid title (subplot title) ONLY
LIPID_TITLE_FS = 50

# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# Big figure (inches per subplot)
PANEL_W = 12.0
PANEL_H = 9.0

# -----------------------------
# PE lipids (3 rows)
# -----------------------------
PLOTS = [
    dict(bb_txt="PE_20_0_18_1.txt", cm_txt="PE_20_0_18_1_CM.txt", rt_min=9.700,  rt_max=10.104, title="(P-20:0/18:1)"),
    dict(bb_txt="PE_20_0_22_6.txt", cm_txt="PE_20_0_22_6_CM.txt", rt_min=10.093, rt_max=10.321, title="(P-20:0/22:6)"),
    dict(bb_txt="PE_20_1_22_6.txt", cm_txt="PE_20_1_22_6_CM.txt", rt_min=10.321, rt_max=10.740, title="(P-20:1/22:6)"),
]

# -----------------------------
# Parsing / stats
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, panel_letter=None):
    ax.set_xlim(rt_min, rt_max)

    # Lipid title: BOLD + custom font size
    ax.set_title(title, fontsize=LIPID_TITLE_FS, fontweight="bold", pad=16)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=14)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=14)

    sci_yaxis(ax)

    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if panel_letter:
        ax.text(-0.18, 1.07, panel_letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks_bb, rt_min, rt_max, title, panel_letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks_bb, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks_bb, C_REPS, rt_min, rt_max)

    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="Lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, panel_letter)


def plot_cm_vs_ctrl(ax, blocks_cm, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks_cm, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks_cm, C_REPS, rt_min, rt_max)

    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="Complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="Lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 3 rows x 2 columns
# -----------------------------
n = len(PLOTS)
fig, axes = plt.subplots(
    nrows=n, ncols=2,
    figsize=(PANEL_W * 2, PANEL_H * n),
    dpi=300,
    constrained_layout=True
)
axes = np.atleast_2d(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    blocks_bb = parse_chrom_txt(cfg["bb_txt"])
    blocks_cm = parse_chrom_txt(cfg["cm_txt"])

    plot_bb_vs_ctrl(axes[i, 0], blocks_bb, cfg["rt_min"], cfg["rt_max"], cfg["title"], panel_letter=letters[i])
    plot_cm_vs_ctrl(axes[i, 1], blocks_cm, cfg["rt_min"], cfg["rt_max"], cfg["title"])

    # Force same y-scale: CM axis limits equal to Bb axis limits
    axes[i, 1].set_ylim(axes[i, 0].get_ylim())

# Legends (top row)
hL, lL = axes[0, 0].get_legend_handles_labels()
hR, lR = axes[0, 1].get_legend_handles_labels()

axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.02, 0.98), borderaxespad=0)
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "PE_Bb_vs_CM_accessible_sameY"
fig.savefig("pE_Bb_vs_CM_accessible_sameY.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("pE_Bb_vs_CM_accessible_sameY.svg", bbox_inches="tight", facecolor="white")
plt.show()







#Plotting EIC plots of Phosphatidylglycerol (PG)#

import re, os, math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.45
LINEWIDTH = 2.5

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# mpl.rcParams["svg.fonttype"] = "path"

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# NEW: ONLY lipid title size (subplot title)
LIPID_TITLE_FS = 50   # <-- change this only to adjust lipid name size

# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# Data
# -----------------------------
PLOTS = [
    dict(txt="PG_14_0_16_0.txt", rt_min=5.054, rt_max=5.304, title="(14:0/16:0)"),
    dict(txt="PG_14_0_18_1.txt", rt_min=5.113, rt_max=5.457, title="(14:0/18:1)"),
    dict(txt="PG_16_0_16_0.txt", rt_min=5.540, rt_max=5.874, title="(16:0/16:0)"),
    dict(txt="PG_16_0_18_0.txt", rt_min=6.133, rt_max=6.454, title="(16:0/18:0)"),
    dict(txt="PG_16_0_18_1.txt", rt_min=5.6,   rt_max=6.0,   title="(16:0/18:1)"),
    dict(txt="PG_16_1_18_1.txt", rt_min=5.182, rt_max=5.646, title="(16:1/18:1)"),
    dict(txt="PG_16_1_18_2.txt", rt_min=4.874, rt_max=5.195, title="(16:1/18:2)"),
    dict(txt="PG_18_0_18_0.txt", rt_min=6.239, rt_max=6.524, title="(18:0/18:0)"),
    dict(txt="PG_18_0_18_1.txt", rt_min=6.239, rt_max=6.572, title="(18:0/18:1)"),
    dict(txt="PG_18_1_18_1.txt", rt_min=5.704, rt_max=6.094, title="(18:1/18:1)"),
    dict(txt="PG_18_1_18_2.txt", rt_min=5.347, rt_max=5.704, title="(18:1/18:2)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue
            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    # ONLY lipid name size is controlled here:
    ax.set_title(title, fontweight="bold",fontsize=LIPID_TITLE_FS, pad=12)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.18, 1.06, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)

    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b - sem_b, mean_b + sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c - sem_c, mean_c + sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)

    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m - sem_m, mean_m + sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c - sem_c, mean_c + sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 4 plots in a row (4 columns)
# -----------------------------
n = len(PLOTS)
ncols = 4
nrows = math.ceil(n / 2)

PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True
)
axes = np.array(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // 2
    c0 = (i % 2) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

    # same y-scale for Bb and CM for each lipid
    axes[r, c0 + 1].set_ylim(axes[r, c0].get_ylim())

# delete unused axes (if odd number of lipids)
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# Legends only on the top-left pair
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.40, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.78, 0.98), borderaxespad=0)

out_base = "PG_two_grids_4_across_ACCESSIBLE_sameY_lipidTitleOnly"
fig.savefig("PG_two_grids_4_across_ACCESSIBLE_sameY_lipidTitleOnly.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig("PG_two_grids_4_across_ACCESSIBLE_sameY_lipidTitleOnly.svg", bbox_inches="tight", facecolor="white")
plt.show()



EIC plots of ceramides (Cer)
import re, os, math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 2.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# mpl.rcParams["svg.fonttype"] = "path"

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# ONLY lipid title size (subplot title)
LIPID_TITLE_FS = 50  # <-- adjust lipid name size here

# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# Ceramide / Cer1P list
# -----------------------------
PLOTS = [
    dict(txt="CER1P_D18_1_16_0.txt", rt_min=10.466, rt_max=10.778, title="Cer-1-P (d18:1/16:0)"),
    dict(txt="CER_D18_1_20_0.txt",   rt_min=10.330, rt_max=10.537, title="(d18:1/20:0)"),
    dict(txt="CER_D18_1_22_0.txt",   rt_min=10.430, rt_max=10.707, title="(d18:1/22:0)"),
    dict(txt="CER_D18_1_23_0.txt",   rt_min=10.499, rt_max=10.706, title="(d18:1/23:0)"),
    dict(txt="CER_D18_1_24_0.txt",   rt_min=10.495, rt_max=10.841, title="(d18:1/24:0)"),
    dict(txt="CER_D18_1_24_1.txt",   rt_min=10.357, rt_max=10.773, title="(d18:1/24:1)"),
    dict(txt="CER_D18_2_24_0.txt",   rt_min=10.627, rt_max=10.981, title="(d18:2/24:0)"),
    dict(txt="CER_T18_0_26_0.txt",   rt_min=10.227, rt_max=10.680, title="(t18:0/26:0)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    # Lipid title: BOLD + custom size
    ax.set_title(title, fontsize=LIPID_TITLE_FS, fontweight="bold", pad=12)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.20, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 4 plots in a row (4 columns)
# -----------------------------
n = len(PLOTS)
lipids_per_row = 2
ncols = 4
nrows = math.ceil(n / lipids_per_row)

# Big figure per panel (increase if needed)
PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True
)
axes = np.array(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // lipids_per_row
    c0 = (i % lipids_per_row) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0],     blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

# delete unused panels
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# legends only on first Bb and first CM axes
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.72, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "CER_grid_4_across_ACCESSIBLE"
fig.savefig("CER_grid_4_across_ACCESSIBLE.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig("CER_grid_4_across_ACCESSIBLE.svg", bbox_inches="tight", facecolor="white")
plt.show()










Plotting EIC of Sphingomyelins (SM)
import re, os, math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 2.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
# mpl.rcParams["svg.fonttype"] = "path"

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# ONLY lipid title size (subplot title)
LIPID_TITLE_FS = 50  # <-- adjust lipid name size here

# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# Sphingomyelin (SM) list
# (titles edited: removed the word "Sphingomyelin")
# -----------------------------
PLOTS = [
    dict(txt="SM_16_1_24_1.txt",   rt_min=9.463, rt_max=10.321, title="(16:1/24:1)"),
    dict(txt="SM_38_3.txt",        rt_min=9.712, rt_max=10.161, title="(38:3)"),
    dict(txt="SM_40_3.txt",        rt_min=9.260, rt_max=9.995,  title="(40:3)"),
    dict(txt="SM_41_0.txt",        rt_min=9.419, rt_max=9.786,  title="(41:0)"),
    dict(txt="SM_41_1.txt",        rt_min=9.991, rt_max=10.562, title="(41:1)"),
    dict(txt="SM_D16_1_23_0.txt",  rt_min=9.874, rt_max=10.364, title="(d16:1/23:0)"),
    dict(txt="SM_D18_0_16_0.txt",  rt_min=8.775, rt_max=9.183,  title="(d18:0/16:0)"),
    dict(txt="SM_D18_0_22_0.txt",  rt_min=10.115,rt_max=10.401, title="(d18:0/22:0)"),
    dict(txt="SM_D18_1_14_0.txt",  rt_min=3.672, rt_max=4.081,  title="(d18:1/14:0)"),
    dict(txt="SM_D18_1_16_0.txt",  rt_min=8.653, rt_max=9.143,  title="(d18:1/16:0)"),
    dict(txt="SM_D18_1_18_0.txt",  rt_min=8.038, rt_max=8.500,  title="(d18:1/18:0)"),
    dict(txt="SM_D18_1_20_0.txt",  rt_min=9.425, rt_max=9.875,  title="(d18:1/20:0)"),
    dict(txt="SM_D18_1_24_1.txt",  rt_min=10.029,rt_max=10.478, title="(d18:1/24:1)"),
    dict(txt="SM_D18_2_16_0.txt",  rt_min=7.673, rt_max=8.286,  title="(d18:2/16:0)"),
    dict(txt="SM_D18_2_20_0.txt",  rt_min=9.508, rt_max=9.793,  title="(d18:2/20:0)"),
    dict(txt="SM_D19_0_23_1.txt",  rt_min=10.070,rt_max=10.560, title="(d19:0/23:1)"),
    dict(txt="SM_D38_4.txt",       rt_min=9.712, rt_max=10.325, title="(d38:4)"),
    dict(txt="SM_D44_5.txt",       rt_min=9.860, rt_max=10.351, title="(d44:5)"),
    dict(txt="SM_D18_2_23_0.txt",  rt_min=9.951, rt_max=10.200, title="(d18:2/23:0)"),
    dict(txt="SM_D18_2_18_0.txt",  rt_min=9.019, rt_max=9.632,  title="(d18:2/18:0)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    # Lipid title: BOLD + custom size
    ax.set_title(title, fontsize=LIPID_TITLE_FS, fontweight="bold", pad=12)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.20, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)

    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 6 plots in a row (6 columns) = 3 lipids per row (each lipid uses 2 columns)
# -----------------------------
n = len(PLOTS)

lipids_per_row = 3          # 3 lipids across
ncols = lipids_per_row * 2  # 2 panels per lipid => 6 columns
nrows = math.ceil(n / lipids_per_row)

# Big figure per panel (increase if needed)
PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True
)

axes = np.array(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // lipids_per_row
    c0 = (i % lipids_per_row) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0],     blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

# delete unused panels
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# legends only on first Bb and first CM axes
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.72, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "SM_grid_6_across_ACCESSIBLE"
fig.savefig("SM.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig("SM.svg", bbox_inches="tight", facecolor="white")
plt.show()








#EIC plots of ten randomly selected phosphatidylcholines (PC)#

import re, os, math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 2.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# ONLY lipid title size (subplot title)
LIPID_TITLE_FS = 50  # <-- adjust lipid name size here

# Make axes/ticks more visible
mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# Phosphatidylcholine (PC) list
# (titles edited: removed "PC ")
# -----------------------------
PLOTS = [
    dict(txt="PC_16_0_16_0.txt",      rt_min=9.859, rt_max=10.172, title="(16:0/16:0)"),
    dict(txt="PC_16_0_18_2.txt",      rt_min=9.225, rt_max=9.812,  title="(16:0/18:2)"),
    dict(txt="PC_16_0_18_1.txt",      rt_min=9.889, rt_max=10.398, title="(16:0/18:1)"),
    dict(txt="PC_16_0_16_1.txt",      rt_min=8.919, rt_max=9.624,  title="(16:0/16:1)"),
    dict(txt="PC_18_0_18_1.txt",      rt_min=9.816, rt_max=10.521, title="(18:0/18:1)"),
    dict(txt="PC_15_0_16_0.txt",      rt_min=9.079, rt_max=9.901,  title="(15:0/16:0)"),
    dict(txt="PC_18_0_18_3.txt",      rt_min=9.449, rt_max=10.115, title="(18:0/18:3)"),
    dict(txt="PC_16_0_17_0.txt",      rt_min=9.933, rt_max=10.324, title="(16:0/17:0)"),
    dict(txt="PC_15_MHDA_18_1.txt",   rt_min=10.007,rt_max=10.335, title="(15:MHDA/18:1)"),
    dict(txt="PC_16_0_20_4.txt",      rt_min=8.546, rt_max=9.093,  title="(16:0/20:4)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    # Lipid title: BOLD + custom size
    ax.set_title(title, fontsize=LIPID_TITLE_FS, fontweight="bold", pad=12)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.20, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 4 plots in a row (4 columns)
# -----------------------------
n = len(PLOTS)
lipids_per_row = 2
ncols = 4
nrows = math.ceil(n / lipids_per_row)

# Big figure per panel (increase if needed)
PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True
)
axes = np.array(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // lipids_per_row
    c0 = (i % lipids_per_row) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0],     blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

# delete unused panels
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# legends only on first Bb and first CM axes
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.72, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "PC_grid_4_across_ACCESSIBLE"
fig.savefig("PC_grid_4_across_ACCESSIBLE.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig("PC_grid_4_across_ACCESSIBLE.svg", bbox_inches="tight", facecolor="white")
plt.show()













#EIC plots of TGs and DG#
# Fix overlapping/ugly titles by:
# 1) Using SHORT titles (no "Triacylglycerol"/"Diacylglycerol" prefix)
# 2) Adding line breaks for long compositions
# 3) Wrapping + smaller font for long titles (auto)
# 4) More title padding

import re, os, math, textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 2.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# Title sizing: base + auto-shrink for long titles
LIPID_TITLE_FS = 50
TITLE_PAD = 18  # more vertical space above axes

mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# TG / DG list (short titles)
# -----------------------------
PLOTS = [
    dict(txt="TG_O_52_2.txt",           rt_min=11.57, rt_max=11.89, title="(O-52:2)"),
    dict(txt="TG_O_50_1.txt",           rt_min=11.45, rt_max=12.12, title="(O-50:1)"),
    dict(txt="TG_50_3.txt",             rt_min=11.16, rt_max=11.94, title="(50:3)"),
    dict(txt="TG_49_1.txt",             rt_min=11.41, rt_max=12.12, title="(49:1)"),
    dict(txt="TG_48_2.txt",             rt_min=11.29, rt_max=11.58, title="(48:2)"),
    dict(txt="TG_18_3_18_2_16_0.txt",   rt_min=11.16, rt_max=11.63, title="(18:3/18:2/\n16:0)"),
    dict(txt="TG_18_1E_16_0_16_0.txt",  rt_min=11.63, rt_max=12.04, title="(18:1e/16:0/\n16:0)"),
    dict(txt="TG_18_0_16_0_18_3.txt",   rt_min=11.30, rt_max=11.85, title="(18:0/16:0/\n18:3)"),
    dict(txt="DG_16_0_18_1.txt",        rt_min=9.30,  rt_max=9.90,  title="DG (16:0/18:1)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1

        i += 1
        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def pretty_title(title, wrap_width=18):
    # Wrap only if user didn't manually add newlines
    if "\n" not in title:
        title = "\n".join(textwrap.wrap(title, width=wrap_width))
    # Auto-shrink if still long
    n_chars = max((len(line) for line in title.split("\n")), default=len(title))
    if n_chars <= 14:
        fs = LIPID_TITLE_FS
    elif n_chars <= 18:
        fs = int(LIPID_TITLE_FS * 0.90)
    elif n_chars <= 24:
        fs = int(LIPID_TITLE_FS * 0.80)
    else:
        fs = int(LIPID_TITLE_FS * 0.70)
    return title, fs


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    t, tfs = pretty_title(title)
    ax.set_title(t, fontsize=tfs, fontweight="bold", pad=TITLE_PAD)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.20, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 4 plots in a row (4 columns)
# -----------------------------
n = len(PLOTS)
lipids_per_row = 2
ncols = 4
nrows = math.ceil(n / lipids_per_row)

PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True
)
axes = np.array(axes)
letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // lipids_per_row
    c0 = (i % lipids_per_row) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0],     blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

# delete unused panels
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# legends only on first Bb and first CM axes
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.72, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "TG_DG_grid_4_across_ACCESSIBLE_titles_fixed"
fig.savefig("tgdg.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig("tgdg.svg", bbox_inches="tight", facecolor="white")
plt.show()















#EIC plots of Cholesteryl esters (CE)#

import re, os, math, textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# -----------------------------
# Settings
# -----------------------------
B_REPS  = ["Bb1", "Bb2", "Bb3", "Bb4"]
CM_REPS = ["CM1", "CM2", "CM3", "CM4"]
C_REPS  = ["Waterblank1", "Waterblank2", "Waterblank3", "Waterblank4"]

B_LINE,  B_BAND  = "#eb7944", "#ce8259"
CM_LINE, CM_BAND = "#517355", "#a3e9ab"
C_LINE,  C_BAND  = "#3539eb", "#8a89d8"
BAND_ALPHA = 0.55
LINEWIDTH = 2.0

mpl.rcParams["font.family"] = "Myriad Pro"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# -----------------------------
# ACCESSIBILITY SIZES
# -----------------------------
LABEL_FS  = 40
TICK_FS_X = 30
TICK_FS_Y = 30
OFFSET_FS = 30
LETTER_FS = 70
LEGEND_FS = 30

# Title sizing: base + auto-shrink for long titles
LIPID_TITLE_FS = 50
TITLE_PAD = 18  # more vertical space above axes

mpl.rcParams["axes.linewidth"] = 1.6
mpl.rcParams["xtick.major.width"] = 1.4
mpl.rcParams["ytick.major.width"] = 1.4
mpl.rcParams["xtick.major.size"] = 7
mpl.rcParams["ytick.major.size"] = 7

# -----------------------------
# Cholesteryl esters (CE) list
# Titles are shortened: only the composition in parentheses.
# -----------------------------
PLOTS = [
    dict(txt="CE_24_1.txt", rt_min=10.6,  rt_max=10.9,  title="(24:1)"),
    dict(txt="CE_16_0.txt", rt_min=11.8,  rt_max=12.1,  title="(16:0)"),
    dict(txt="CE_16_1.txt", rt_min=11.6,  rt_max=11.9,  title="(16:1)"),
    dict(txt="CE_17_0.txt", rt_min=11.9,  rt_max=12.2,  title="(17:0)"),
    dict(txt="CE_17_1.txt", rt_min=11.7,  rt_max=12.0,  title="(17:1)"),
    dict(txt="CE_18_0.txt", rt_min=12.05, rt_max=12.25, title="(18:0)"),
    dict(txt="CE_18_1.txt", rt_min=11.85, rt_max=12.10, title="(18:1)"),
    dict(txt="CE_18_2.txt", rt_min=11.6,  rt_max=11.9,  title="(18:2)"),
    dict(txt="CE_18_3.txt", rt_min=10.65, rt_max=11.0,  title="(18:3)"),
    dict(txt="CE_20_1.txt", rt_min=10.3,  rt_max=10.9,  title="(20:1)"),
    dict(txt="CE_20_3.txt", rt_min=9.9,   rt_max=10.2,  title="(20:3)"),
    dict(txt="CE_20_4.txt", rt_min=11.5,  rt_max=11.85, title="(20:4)"),
    dict(txt="CE_22_1.txt", rt_min=10.6,  rt_max=10.9,  title="(22:1)"),
    dict(txt="CE_22_6.txt", rt_min=10.2,  rt_max=10.7,  title="(22:6)"),
    dict(txt="CE_24_0.txt", rt_min=10.6,  rt_max=10.9,  title="(24:0)"),
]

# -----------------------------
# Helpers
# -----------------------------
def parse_chrom_txt(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}\nWorking dir: {os.getcwd()}")

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_re = re.compile(r'([A-Za-z0-9]+)_[^"\s]*?\.d', re.IGNORECASE)

    blocks = {}
    i = 0
    while i < len(lines):
        m = header_re.search(lines[i])
        if not m:
            i += 1
            continue

        rep = m.group(1)

        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("#Point"):
            i += 1
        i += 1

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue

            if header_re.search(s) and not s.lstrip().startswith("#Point"):
                break
            if s.startswith("#"):
                i += 1
                continue

            parts = re.split(r"\s+|\t+|,", s)
            if len(parts) >= 3:
                try:
                    rows.append((float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
            i += 1

        if rows:
            blocks[rep] = (
                pd.DataFrame(rows, columns=["rt", "counts"])
                .sort_values("rt")
                .reset_index(drop=True)
            )

    return blocks


def sem(a, axis=0):
    a = np.asarray(a, dtype=float)
    n = np.sum(~np.isnan(a), axis=axis)
    return np.nanstd(a, axis=axis, ddof=1) / np.sqrt(n)


def extract_group(blocks, rep_ids, rt_min, rt_max):
    dfs = []
    for r in rep_ids:
        if r in blocks:
            df = blocks[r]
            df = df[(df["rt"] >= rt_min) & (df["rt"] <= rt_max)].copy()
            if len(df) > 1:
                dfs.append(df.sort_values("rt"))

    if not dfs:
        raise ValueError(f"No replicates found for {rep_ids}. Found keys: {sorted(blocks.keys())}")

    rt_grid = dfs[0]["rt"].to_numpy()
    y_mat = np.vstack([np.interp(rt_grid, d["rt"].to_numpy(), d["counts"].to_numpy()) for d in dfs])
    return rt_grid, y_mat.mean(axis=0), sem(y_mat, axis=0)


def sci_yaxis(ax):
    sf = ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(sf)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_size(OFFSET_FS)


def pretty_title(title, wrap_width=18):
    if "\n" not in title:
        title = "\n".join(textwrap.wrap(title, width=wrap_width))
    n_chars = max((len(line) for line in title.split("\n")), default=len(title))
    if n_chars <= 14:
        fs = LIPID_TITLE_FS
    elif n_chars <= 18:
        fs = int(LIPID_TITLE_FS * 0.90)
    elif n_chars <= 24:
        fs = int(LIPID_TITLE_FS * 0.80)
    else:
        fs = int(LIPID_TITLE_FS * 0.70)
    return title, fs


def style_ax(ax, rt_min, rt_max, title, letter=None):
    ax.set_xlim(rt_min, rt_max)

    t, tfs = pretty_title(title)
    ax.set_title(t, fontsize=tfs, fontweight="bold", pad=TITLE_PAD)

    ax.set_xlabel("retention time (min)", fontsize=LABEL_FS, labelpad=10)
    ax.set_ylabel("extracted ion counts", fontsize=LABEL_FS, labelpad=10)

    sci_yaxis(ax)
    ax.tick_params(axis="x", labelsize=TICK_FS_X)
    ax.tick_params(axis="y", labelsize=TICK_FS_Y)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if letter:
        ax.text(-0.20, 1.05, letter, transform=ax.transAxes,
                fontsize=LETTER_FS, fontweight="bold", va="bottom")


def plot_bb_vs_ctrl(ax, blocks, rt_min, rt_max, title, letter=None):
    rt_b, mean_b, sem_b = extract_group(blocks, B_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_b, rt_c, mean_c)
    sem_c  = np.interp(rt_b, rt_c, sem_c)

    ax.plot(rt_b, mean_b, color=B_LINE, lw=LINEWIDTH, label=r"$\it{B.\ burgdorferi}$")
    ax.fill_between(rt_b, mean_b-sem_b, mean_b+sem_b, color=B_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_b, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_b, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title, letter)


def plot_cm_vs_ctrl(ax, blocks, rt_min, rt_max, title):
    rt_m, mean_m, sem_m = extract_group(blocks, CM_REPS, rt_min, rt_max)
    rt_c, mean_c, sem_c = extract_group(blocks, C_REPS, rt_min, rt_max)
    mean_c = np.interp(rt_m, rt_c, mean_c)
    sem_c  = np.interp(rt_m, rt_c, sem_c)

    ax.plot(rt_m, mean_m, color=CM_LINE, lw=LINEWIDTH, label="complete medium")
    ax.fill_between(rt_m, mean_m-sem_m, mean_m+sem_m, color=CM_BAND, alpha=BAND_ALPHA, linewidth=0)
    ax.plot(rt_m, mean_c, color=C_LINE, lw=LINEWIDTH, label="lipid extraction control")
    ax.fill_between(rt_m, mean_c-sem_c, mean_c+sem_c, color=C_BAND, alpha=BAND_ALPHA, linewidth=0)

    style_ax(ax, rt_min, rt_max, title)


# -----------------------------
# Run: 4 plots in a row (4 columns)
# (2 lipids per row => 4 columns: Bb/ctrl + CM/ctrl for each lipid)
# -----------------------------
n = len(PLOTS)
lipids_per_row = 2
ncols = 4
nrows = math.ceil(n / lipids_per_row)

PANEL_W = 11.0
PANEL_H = 8.5

fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(PANEL_W * ncols, PANEL_H * nrows),
    dpi=300,
    constrained_layout=True,
    squeeze=False
)

letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

for i, cfg in enumerate(PLOTS):
    r = i // lipids_per_row
    c0 = (i % lipids_per_row) * 2

    blocks = parse_chrom_txt(cfg["txt"])

    plot_bb_vs_ctrl(axes[r, c0],     blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"], letter=letters[i])
    plot_cm_vs_ctrl(axes[r, c0 + 1], blocks, cfg["rt_min"], cfg["rt_max"], cfg["title"])

# delete unused panels
used_axes = n * 2
flat = axes.ravel()
for j in range(used_axes, nrows * ncols):
    fig.delaxes(flat[j])

# legends only on first Bb and first CM axes
hL, lL = axes[0, 0].get_legend_handles_labels()
axes[0, 0].legend(hL, lL, frameon=False, fontsize=LEGEND_FS,
                  loc="upper left", bbox_to_anchor=(0.72, 0.98), borderaxespad=0)

hR, lR = axes[0, 1].get_legend_handles_labels()
axes[0, 1].legend(hR, lR, frameon=False, fontsize=LEGEND_FS,
                  loc="upper right", bbox_to_anchor=(0.98, 0.98), borderaxespad=0)

out_base = "CE_grid_4_across_ACCESSIBLE_titles_fixed"
fig.savefig(f"{out_base}.png", dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(f"{out_base}.svg", bbox_inches="tight", facecolor="white")
plt.show()











