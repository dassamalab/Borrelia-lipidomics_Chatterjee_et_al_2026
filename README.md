#EIC plots of acylated cholesterol galactolipids ran on quadrupole time of flight liquid chromatography mass spectrometry#
"""
Ac Gal 16:0 

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), restricted to RT 4.543–4.992, EIC(804.6687)
RIGHT : MS2 (peak list) showing ONLY m/z 369.3561 and 804.6687 (annotated)

Inputs (per-sample MS1 EIC txt exports):
  263_1_AcGal_16_1.txt
  263_2_AcGal_16_1.txt
  263_3_AcGal_16_1.txt
  h2o_1_AcGal_16_1.txt
  h2o_2_AcGal_16_1.txt
  h2o_3_AcGal_16_1.txt

MS2 peak list txt export:
  Ac Gal 16_0_ms2.txt   (adjust MS2_TXT if different)

Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(804.6687)"

BB_FILES = [
    Path("263_1_AcGal_16_1.txt"),
    Path("263_2_AcGal_16_1.txt"),
    Path("263_3_AcGal_16_1.txt"),
]
H2O_FILES = [
    Path("h2o_1_AcGal_16_1.txt"),
    Path("h2o_2_AcGal_16_1.txt"),
    Path("h2o_3_AcGal_16_1.txt"),
]

MS2_TXT = Path("AcGal_16_0_ms2.txt")  # change if needed

RT_MIN, RT_MAX = 4.543, 4.992
N_GRID = 700

MS2_MZ_TARGETS = [369.3561, 804.6687]
MZ_TOL = 0.01  # Da

# Colors (match your style)
BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

# Output
OUT_PNG = "AcGal16_0_MS1left_MS2right_clean.png"
DPI = 600

# Typography (tuned to avoid overlap)
TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """
    Parse MassHunter TXT blocks:
      #"+ <desc> Scan <sample>.d "
      #Point X(...) Y(...)
      ...

    Returns list of dicts {desc, scan, df} with df columns ['x','y'].
    """
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()
        scan = m.group(2).strip().replace(".d", "")

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    blocks = parse_masshunter_blocks(txt_path)
    eics = [b for b in blocks if b["desc"].strip() == eic_target]
    if not eics:
        found = sorted({b["desc"] for b in blocks if b["desc"].startswith("EIC(")})
        raise RuntimeError(
            f"No EIC block '{eic_target}' found in {txt_path.name}.\nEIC blocks found: {found}"
        )
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    """
    Reads an MS2 peak list exported as text.
    Accepts:
      - 2 columns: mz intensity
      - 3 columns: index mz intensity
    Ignores lines starting with '#'.
    """
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    # existence checks
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # ---- MS1 traces (per-sample files) ----
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # ---- MS2 selected ions (peak list) ----
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # -----------------------------
    # Figure (1×2): MS1 left, MS2 right
    # -----------------------------
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.0, 4.2), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]}
    )

    # ---------- MS1 panel ----------
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)

    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(Ac Gal 16:0)", fontsize=TITLE_FS, fontweight="bold", pad=10)

    ax1.tick_params(axis="both", labelsize=TICK_FS)

    # scientific notation ×10^n, smaller offset text
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)

    beautify_axis(ax1)

    # place legend away from peak
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---------- MS2 panel ----------
    # Set x-limits with padding so labels don't touch the edges
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    # y-limit: add headroom for labels
    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        # nudge text inward to avoid edges
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)
    fig.savefig("AcGal_16_0.png", dpi=600, bbox_inches="tight")
    fig.savefig("AcGal_16_0.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()





"""
Ac Gal 18:0 — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), restricted to RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(832.6961)"

BB_FILES = [
    Path("263_1_AcGal_18_0.txt"),
    Path("263_2_AcGal_18_0.txt"),
    Path("263_3_AcGal_18_0.txt"),
]
H2O_FILES = [
    Path("h2o_1_AcGal_18_0.txt"),
    Path("h2o_2_AcGal_18_0.txt"),
    Path("h2o_3_AcGal_18_0.txt"),
]

MS2_TXT = Path("AcGal_18_0_ms2.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.546, 4.998
N_GRID = 700

MS2_MZ_TARGETS = [369.3542, 832.7052]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

# Output base name (no extension)
OUT_BASE = "AcGal18_0_MS1left_MS2right_clean"
DPI = 600

# Typography
TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()
        scan = m.group(2).strip().replace(".d", "")

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    blocks = parse_masshunter_blocks(txt_path)
    eics = [b for b in blocks if b["desc"].strip() == eic_target]
    if not eics:
        found = sorted({b["desc"] for b in blocks if b["desc"].startswith("EIC(")})
        raise RuntimeError(
            f"No EIC block '{eic_target}' found in {txt_path.name}.\nEIC blocks found: {found}"
        )
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    """
    Reads MS2 peak list exported as text.
    Accepts:
      - 2 columns: mz intensity
      - 3 columns: index mz intensity
    Ignores lines starting with '#'.
    """
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(Ac Gal 18:0)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("AcGal_18_0.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("AcGal_18_0.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()




"""
Ac Gal 18:1 — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), restricted to RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(830.6841)"

BB_FILES = [
    Path("263_1_AcGal_18_1.txt"),
    Path("263_2_AcGal_18_1.txt"),
    Path("263_3_AcGal_18_1.txt"),
]
H2O_FILES = [
    Path("h2o_1_AcGal_18_1.txt"),
    Path("h2o_2_AcGal_18_1.txt"),
    Path("h2o_3_AcGal_18_1.txt"),
]

MS2_TXT = Path("AcGal_18_1_ms2.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.543, 4.905
N_GRID = 700

MS2_MZ_TARGETS = [369.3511, 830.6832]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

# Output base name (no extension)
OUT_BASE = "AcGal18_1_MS1left_MS2right_clean"
DPI = 600

# Typography
TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()
        scan = m.group(2).strip().replace(".d", "")

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    blocks = parse_masshunter_blocks(txt_path)
    eics = [b for b in blocks if b["desc"].strip() == eic_target]
    if not eics:
        found = sorted({b["desc"] for b in blocks if b["desc"].startswith("EIC(")})
        raise RuntimeError(
            f"No EIC block '{eic_target}' found in {txt_path.name}.\nEIC blocks found: {found}"
        )
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    """
    Reads MS2 peak list exported as text.
    Accepts:
      - 2 columns: mz intensity
      - 3 columns: index mz intensity
    Ignores lines starting with '#'.
    """
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(Ac Gal 18:1)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("AcGal_18_1.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("AcGal_18_1.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


"""
Ac Gal 18:2 — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), restricted to RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(828.6705)"

BB_FILES = [
    Path("263_1_AcGal_18_2.txt"),
    Path("263_2_AcGal_18_2.txt"),
    Path("263_3_AcGal_18_2.txt"),
]
H2O_FILES = [
    Path("h2o_1_AcGal_18_2.txt"),
    Path("h2o_2_AcGal_18_2.txt"),
    Path("h2o_3_AcGal_18_2.txt"),
]

MS2_TXT = Path("AcGal_18_2_ms2.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.510, 4.881
N_GRID = 700

MS2_MZ_TARGETS = [369.3539, 828.6755]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

# Output base name (no extension)
OUT_BASE = "AcGal18_2_MS1left_MS2right_clean"
DPI = 600

# Typography
TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()
        scan = m.group(2).strip().replace(".d", "")

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    blocks = parse_masshunter_blocks(txt_path)
    eics = [b for b in blocks if b["desc"].strip() == eic_target]
    if not eics:
        found = sorted({b["desc"] for b in blocks if b["desc"].startswith("EIC(")})
        raise RuntimeError(
            f"No EIC block '{eic_target}' found in {txt_path.name}.\nEIC blocks found: {found}"
        )
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    """
    Reads MS2 peak list exported as text.
    Accepts:
      - 2 columns: mz intensity
      - 3 columns: index mz intensity
    Ignores lines starting with '#'.
    """
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(Ac Gal 18:2)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("AcGal_18_2.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("AcGal_18_2.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()










EIC plots of Monogalactosyldiacylglycerols
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
MGDG(16:0/18:2) — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(772.5889)"  # will be matched as substring within desc (e.g., '+ESI EIC(772.5889)')

BB_FILES = [
    Path("263_1_MGDG_16_0_18_2.txt"),
    Path("263_2_MGDG_16_0_18_2.txt"),
    Path("263_3_MGDG_16_0_18_2.txt"),
]
H2O_FILES = [
    Path("H2O_1_MGDG_16_0_18_2.txt"),
    Path("H2O_2_MGDG_16_0_18_2.txt"),
    Path("H2O_3_MGDG_16_0_18_2.txt"),
]

MS2_TXT = Path("MS2_MGDG_16_0_18_2.txt")  # peak list export (text)

RT_MIN, RT_MAX = 3.840, 4.572
N_GRID = 700

MS2_MZ_TARGETS = [313.2785, 337.2742, 575.4975, 593.5097, 772.5923]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

OUT_BASE = "MGDG_16_0_18_2_MS1left_MS2right_clean"
DPI = 600

TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()                 # e.g., '+ESI EIC(772.5889) Scan Frag=150.0V'
        scan = m.group(2).strip().replace(".d", "")  # e.g., '263_1'

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    """
    Return the EIC dataframe matching eic_target as a substring of the descriptor.
    This handles headers like '+ESI EIC(772.5889) ...'
    """
    blocks = parse_masshunter_blocks(txt_path)

    # substring match instead of strict equality
    eics = [b for b in blocks if eic_target in b["desc"]]
    if not eics:
        found = sorted({b["desc"] for b in blocks})
        raise RuntimeError(
            f"No EIC block containing '{eic_target}' found in {txt_path.name}.\n"
            f"Block descriptors found (first 5): {found[:5]}"
        )

    # If multiple matches exist, take the first
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(MGDG 16:0/18:2)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig(f"{OUT_BASE}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{OUT_BASE}.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


if __name__ == "__main__":
    main()



"""
MGDG(16:0/18:2) — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(800.6211)"  # will be matched as substring within desc (e.g., '+ESI EIC(772.5889)')

BB_FILES = [
    Path("263_1_MGDG_18_1_18_1.txt"),
    Path("263_2_MGDG_18_1_18_1.txt"),
    Path("263_3_MGDG_18_1_18_1.txt"),
]
H2O_FILES = [
    Path("H2O_1_MGDG_18_1_18_1.txt"),
    Path("H2O_2_MGDG_18_1_18_1.txt"),
    Path("H2O_3_MGDG_18_1_18_1.txt"),
]

MS2_TXT = Path("MS2_MGDG_18_1_18_1.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.167, 4.637
N_GRID = 700

MS2_MZ_TARGETS = [339.291, 603.5358, 621.5431, 800.6201]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

OUT_BASE = "MGDG_18_1_18_1_MS1left_MS2right_clean"
DPI = 600

TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()                 # e.g., '+ESI EIC(772.5889) Scan Frag=150.0V'
        scan = m.group(2).strip().replace(".d", "")  # e.g., '263_1'

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    """
    Return the EIC dataframe matching eic_target as a substring of the descriptor.
    This handles headers like '+ESI EIC(772.5889) ...'
    """
    blocks = parse_masshunter_blocks(txt_path)

    # substring match instead of strict equality
    eics = [b for b in blocks if eic_target in b["desc"]]
    if not eics:
        found = sorted({b["desc"] for b in blocks})
        raise RuntimeError(
            f"No EIC block containing '{eic_target}' found in {txt_path.name}.\n"
            f"Block descriptors found (first 5): {found[:5]}"
        )

    # If multiple matches exist, take the first
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(MGDG 18:1/18:1)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("MGDG_18_1_18_1.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("MGDG_18_1_18_1.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


if __name__ == "__main__":
    main()




    ''''''''''''''''''''''''''''''''''
    """
MGDG(16:0/18:2) — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(800.6211)"  # will be matched as substring within desc (e.g., '+ESI EIC(772.5889)')

BB_FILES = [
    Path("263_1_MGDG_18_1_18_1.txt"),
    Path("263_2_MGDG_18_1_18_1.txt"),
    Path("263_3_MGDG_18_1_18_1.txt"),
]
H2O_FILES = [
    Path("H2O_1_MGDG_18_1_18_1.txt"),
    Path("H2O_2_MGDG_18_1_18_1.txt"),
    Path("H2O_3_MGDG_18_1_18_1.txt"),
]

MS2_TXT = Path("MS2_MGDG_18_1_18_1.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.167, 4.637
N_GRID = 700

MS2_MZ_TARGETS = [339.291, 603.5358, 621.5431, 800.6201]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

OUT_BASE = "MGDG_18_1_18_1_MS1left_MS2right_clean"
DPI = 600

TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()                 # e.g., '+ESI EIC(772.5889) Scan Frag=150.0V'
        scan = m.group(2).strip().replace(".d", "")  # e.g., '263_1'

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    """
    Return the EIC dataframe matching eic_target as a substring of the descriptor.
    This handles headers like '+ESI EIC(772.5889) ...'
    """
    blocks = parse_masshunter_blocks(txt_path)

    # substring match instead of strict equality
    eics = [b for b in blocks if eic_target in b["desc"]]
    if not eics:
        found = sorted({b["desc"] for b in blocks})
        raise RuntimeError(
            f"No EIC block containing '{eic_target}' found in {txt_path.name}.\n"
            f"Block descriptors found (first 5): {found[:5]}"
        )

    # If multiple matches exist, take the first
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(MGDG 18:1/18:1)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("MGDG_18_1_18_1.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("MGDG_18_1_18_1.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


if __name__ == "__main__":
    main()

''''''''''''''''''''''''''''''''''
"""
MGDG(16:0/18:1) — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(774.6044)"  # will be matched as substring within desc (e.g., '+ESI EIC(772.5889)')

BB_FILES = [
    Path("263_1_MGDG_16_0_18_1.txt"),
    Path("263_2_MGDG_16_0_18_1.txt"),
    Path("263_3_MGDG_16_0_18_1.txt"),
]
H2O_FILES = [
    Path("H2O_1_MGDG_16_0_18_1.txt"),
    Path("H2O_2_MGDG_16_0_18_1.txt"),
    Path("H2O_3_MGDG_16_0_18_1.txt"),
]

MS2_TXT = Path("MS2_MGDG_16_0_18_1.txt")  # peak list export (text)

RT_MIN, RT_MAX = 4.079, 4.752
N_GRID = 700

MS2_MZ_TARGETS = [313.2726, 339.2924, 577.5245, 595.5358, 774.6107]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

OUT_BASE = "MGDG_16_0_18_1_MS1left_MS2right_clean"
DPI = 600

TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()                 # e.g., '+ESI EIC(772.5889) Scan Frag=150.0V'
        scan = m.group(2).strip().replace(".d", "")  # e.g., '263_1'

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    """
    Return the EIC dataframe matching eic_target as a substring of the descriptor.
    This handles headers like '+ESI EIC(772.5889) ...'
    """
    blocks = parse_masshunter_blocks(txt_path)

    # substring match instead of strict equality
    eics = [b for b in blocks if eic_target in b["desc"]]
    if not eics:
        found = sorted({b["desc"] for b in blocks})
        raise RuntimeError(
            f"No EIC block containing '{eic_target}' found in {txt_path.name}.\n"
            f"Block descriptors found (first 5): {found[:5]}"
        )

    # If multiple matches exist, take the first
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(MGDG 16:0/18:1)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("MGDG_16_0_18_1.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("MGDG_16_0_18_1.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


if __name__ == "__main__":
    main()

"""
MGDG(16:0/16:1) — clean, non-overlapping 1×2 figure

LEFT  : MS1 EIC mean ± SEM (Bb vs H2O), RT_MIN–RT_MAX, EIC_TARGET
RIGHT : MS2 (peak list) showing ONLY specified m/z values (annotated)

Saves BOTH PNG and SVG.
Requires: numpy, pandas, matplotlib
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# -----------------------------
# User settings
# -----------------------------
EIC_TARGET = "EIC(746.5668)"  # will be matched as substring within desc (e.g., '+ESI EIC(772.5889)')

BB_FILES = [
    Path("263_1_MGDG_16_0_16_1.txt"),
    Path("263_2_MGDG_16_0_16_1.txt"),
    Path("263_3_MGDG_16_0_16_1.txt"),
]
H2O_FILES = [
    Path("H2O_1_MGDG_16_0_16_1.txt"),
    Path("H2O_2_MGDG_16_0_16_1.txt"),
    Path("H2O_3_MGDG_16_0_16_1.txt"),
]

MS2_TXT = Path("4_MS2_MGDG_16_0_16_1.txt")  # peak list export (text)

RT_MIN, RT_MAX = 3.742, 4.290
N_GRID = 700

MS2_MZ_TARGETS = [311.2601, 313.2732, 549.4869, 567.506, 746.5650]
MZ_TOL = 0.01  # Da

BB_COLOR = "#f26c2a"   # orange
H2O_COLOR = "#0000ff"  # blue

OUT_BASE = "MGDG_16_0_16_1_MS1left_MS2right_clean"
DPI = 600

TITLE_FS = 26
AXLABEL_FS = 18
YLAB_FS = 20
TICK_FS = 12
LEGEND_FS = 15
SMALL_FS = 11


# -----------------------------
# Matplotlib style
# -----------------------------
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "svg.fonttype": "none",  # keep editable text in SVG
    }
)


# -----------------------------
# MassHunter EIC block parser
# -----------------------------
HEADER_RE = re.compile(r'#"\+\s*(.*?)\s+Scan\s+(.+?)\s*"')


def parse_masshunter_blocks(txt_path: Path) -> list[dict]:
    """Parse MassHunter TXT blocks and return list of dicts {desc, scan, df}."""
    lines = txt_path.read_text(errors="ignore").splitlines()
    blocks: list[dict] = []

    i = 0
    while i < len(lines):
        m = HEADER_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue

        desc = m.group(1).strip()                 # e.g., '+ESI EIC(772.5889) Scan Frag=150.0V'
        scan = m.group(2).strip().replace(".d", "")  # e.g., '263_1'

        # find "#Point"
        i += 1
        while i < len(lines) and not lines[i].startswith("#Point"):
            i += 1
        if i >= len(lines):
            break
        i += 1  # skip "#Point ..."

        rows = []
        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("#"):
                break
            parts = s.split()
            if len(parts) >= 3:
                rows.append((float(parts[1]), float(parts[2])))
            i += 1

        if rows:
            df = pd.DataFrame(rows, columns=["x", "y"])
            blocks.append({"desc": desc, "scan": scan, "df": df})

        i += 1

    return blocks


def get_eic_df(txt_path: Path, eic_target: str) -> pd.DataFrame:
    """
    Return the EIC dataframe matching eic_target as a substring of the descriptor.
    This handles headers like '+ESI EIC(772.5889) ...'
    """
    blocks = parse_masshunter_blocks(txt_path)

    # substring match instead of strict equality
    eics = [b for b in blocks if eic_target in b["desc"]]
    if not eics:
        found = sorted({b["desc"] for b in blocks})
        raise RuntimeError(
            f"No EIC block containing '{eic_target}' found in {txt_path.name}.\n"
            f"Block descriptors found (first 5): {found[:5]}"
        )

    # If multiple matches exist, take the first
    return eics[0]["df"]


# -----------------------------
# MS1 processing
# -----------------------------
def interp_trace(df: pd.DataFrame, rt_grid: np.ndarray, rt_min: float, rt_max: float) -> np.ndarray:
    w = df[(df["x"] >= rt_min) & (df["x"] <= rt_max)].sort_values("x")
    if w.empty:
        return np.zeros_like(rt_grid)
    return np.interp(rt_grid, w["x"].to_numpy(), w["y"].to_numpy(), left=0.0, right=0.0)


def mean_and_sem(traces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = traces.mean(axis=0)
    sem = traces.std(axis=0, ddof=1) / np.sqrt(traces.shape[0]) if traces.shape[0] > 1 else np.zeros_like(mean)
    return mean, sem


# -----------------------------
# MS2 processing (peak list TXT)
# -----------------------------
def read_ms2_peaklist(txt_path: Path) -> pd.DataFrame:
    rows = []
    for line in txt_path.read_text(errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) == 2:
            mz, inten = parts
        elif len(parts) >= 3:
            mz, inten = parts[-2], parts[-1]
        else:
            continue
        try:
            rows.append((float(mz), float(inten)))
        except ValueError:
            continue

    if not rows:
        raise RuntimeError(f"No numeric MS2 peaks parsed from: {txt_path.name}")

    return pd.DataFrame(rows, columns=["mz", "intensity"])


def intensity_near_mz(ms2: pd.DataFrame, mz0: float, tol: float) -> float:
    hit = ms2[(ms2["mz"] >= mz0 - tol) & (ms2["mz"] <= mz0 + tol)]
    return float(hit["intensity"].max()) if not hit.empty else 0.0


# -----------------------------
# Plot helpers
# -----------------------------
def beautify_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# -----------------------------
# Main
# -----------------------------
def main() -> None:
    for p in BB_FILES + H2O_FILES + [MS2_TXT]:
        if not p.exists():
            raise FileNotFoundError(f"Missing file: {p.resolve()}")

    rt_grid = np.linspace(RT_MIN, RT_MAX, N_GRID)

    # MS1 traces
    bb_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in BB_FILES])
    h2o_traces = np.vstack([interp_trace(get_eic_df(f, EIC_TARGET), rt_grid, RT_MIN, RT_MAX) for f in H2O_FILES])

    bb_mean, bb_sem = mean_and_sem(bb_traces)
    h2o_mean, h2o_sem = mean_and_sem(h2o_traces)

    # MS2 selected ions
    ms2 = read_ms2_peaklist(MS2_TXT)
    ms2_sel = pd.DataFrame(
        {
            "mz": MS2_MZ_TARGETS,
            "intensity": [intensity_near_mz(ms2, mz0, MZ_TOL) for mz0 in MS2_MZ_TARGETS],
        }
    )

    # Figure
    fig, (ax1, ax2) = plt.subplots(
        1, 2,
        figsize=(12.0, 4.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 1.0]},
    )

    # ---- MS1 (left) ----
    ax1.plot(rt_grid, bb_mean, color=BB_COLOR, lw=2.4, label=r"$\it{B.\ burgdorferi}$", zorder=3)
    ax1.fill_between(rt_grid, bb_mean - bb_sem, bb_mean + bb_sem, color=BB_COLOR, alpha=0.28, linewidth=0, zorder=2)

    ax1.plot(rt_grid, h2o_mean, color=H2O_COLOR, lw=2.4, label="lipid extraction control", zorder=3)
    ax1.fill_between(
        rt_grid, h2o_mean - h2o_sem, h2o_mean + h2o_sem, color=H2O_COLOR, alpha=0.18, linewidth=0, zorder=1
    )

    ax1.set_xlim(RT_MIN, RT_MAX)
    ax1.set_xlabel("retention time (min)", fontsize=AXLABEL_FS, labelpad=6)
    ax1.set_ylabel("extracted ion counts", fontsize=YLAB_FS, labelpad=10)
    ax1.set_title("(MGDG 16:0/16:1)", fontsize=TITLE_FS, fontweight="bold", pad=10)
    ax1.tick_params(axis="both", labelsize=TICK_FS)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    ax1.yaxis.get_offset_text().set_size(SMALL_FS)
    beautify_axis(ax1)
    ax1.legend(
        frameon=False,
        fontsize=LEGEND_FS,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        borderaxespad=0.2,
        handlelength=2.2,
    )

    # ---- MS2 (right) ----
    xmin = min(MS2_MZ_TARGETS) - 35
    xmax = max(MS2_MZ_TARGETS) + 35
    ax2.set_xlim(xmin, xmax)

    ax2.vlines(ms2_sel["mz"], 0, ms2_sel["intensity"], color="black", lw=2.2, zorder=2)
    ax2.scatter(ms2_sel["mz"], ms2_sel["intensity"], color="black", s=30, zorder=3)

    ymax = float(ms2_sel["intensity"].max()) * 1.15 if ms2_sel["intensity"].max() > 0 else 1.0
    ax2.set_ylim(0, ymax)

    for mz, inten in ms2_sel.to_numpy():
        dx = -12 if mz < np.mean(MS2_MZ_TARGETS) else 12
        ax2.annotate(
            f"{mz:.4f}",
            xy=(mz, inten),
            xytext=(dx, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=SMALL_FS,
        )

    ax2.set_xlabel("m/z", fontsize=AXLABEL_FS, labelpad=6)
    ax2.set_ylabel("MS/MS intensity", fontsize=AXLABEL_FS, labelpad=10)
    ax2.set_title("MS2 (selected ions)", fontsize=16, pad=10)
    ax2.tick_params(axis="both", labelsize=TICK_FS)
    beautify_axis(ax2)

    # Save BOTH formats
    fig.savefig("MGDG_16_0_16_1.png", dpi=DPI, bbox_inches="tight")
    fig.savefig("MGDG_16_0_16_1.svg", bbox_inches="tight")  # dpi ignored for SVG
    plt.show()


if __name__ == "__main__":
    main()
    
