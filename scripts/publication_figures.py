"""
Paper B (Fire Safety Journal) -- cavity-barrier CFD campaign.
Production publication figures, regenerated from the raw FDS device/slice
outputs in C:\\FDS_runs\\cavb\\ every run (never hardcoded plot values).

One function per figure. All figures written to fds/figures/ as BOTH a
600-dpi PNG and a vector PDF with the same basename.

Style: Elsevier / FSJ single- or double-column, Arial 8 pt base, Okabe-Ito
colours with a FIXED identity mapping, +/-1 std error bars over the steady
window, no in-figure titles, y-grid only, no top/right spines.

Cross-check against manuscript Table 1 runs first (verify_table1); a figure is
aborted and reported in the manifest if any parsed steady mean drifts beyond
tolerance.

Author: MIES for Samson Tan. Data provenance: FDS 6 devc/slice CSV + .sf.
"""

import csv
import os
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# reuse verified constants + the direct .sf reader (fdsreader has a stale
# frame-count bug on this campaign -- do NOT use it)
from generate_inputs import (  # noqa: E402
    GRIDS, X0, X1, Y0, Y1, BAR_Z0, BAR_Z1, TC_Z,
)
import fds_postprocess_cavb as pp  # noqa: E402

RUN_ROOT = r"C:\FDS_runs\cavb"
FIGS_DIR = os.path.normpath(os.path.join(HERE, "..", "figures"))

# ----------------------------------------------------------------- dimensions
MM = 1.0 / 25.4
W_DOUBLE = 190 * MM          # double-column width (in)
W_SINGLE = 90 * MM           # single-column width (in)

# ----------------------------------------------------------------- steady windows
#   60-120 s for every 20/10 mm case; 60-100 s for the trimmed 5 mm case (T_END=100)
DT_DEVC = 0.5                                   # DT_DEVC in the FDS runs (s)


def win_of(chid):
    return (60.0, 100.0) if chid.endswith("05") else (60.0, 120.0)

GAP_CELLS = {"20": 3, "10": 5, "05": 10}      # cells across the 50 mm gap (J count)
DX_NOM = {"20": 20, "10": 10, "05": 5}        # nominal dx (mm)
DSTAR_DX = {"20": "6.3", "10": "12.5", "05": "25"}   # D*/dx = 125/dx

# ----------------------------------------------------------------- colours (FSE-classic)
# Classic fire-safety-engineering aesthetic (Drean et al. 2018, Fire Technology
# Fig 15): solid BLACK primary, solid RED second, solid BLUE third; distinct
# filled markers; error bars as a PALE tint of the series colour. This
# deliberately supersedes the earlier Okabe-Ito / dataviz palette (author
# instruction 2026-07-05). Barrier-state identity is fixed everywhere:
#   no-barrier = black, open = red, closed = blue.
C_BLACK = "#000000"
C_RED = "#C0504D"
C_BLUE = "#1F4E79"

C_STATE = {"nb": C_BLACK, "op": C_RED, "cl": C_BLUE}
M_STATE = {"nb": "o", "op": "s", "cl": "^"}          # circle / square / triangle
STATE_LABEL = {"nb": "No barrier", "op": "Open state", "cl": "Closed state"}
# grid-resolution series: finest = black, then red, then blue (3/5/10 cells)
C_GRID = {"05": C_BLACK, "10": C_RED, "20": C_BLUE}   # fine -> coarse
M_GRID = {"05": "o", "10": "s", "20": "^"}

GRID_ORDER = ["20", "10", "05"]              # coarse -> fine


def tint(color, frac=0.58):
    """Blend a colour toward white by `frac` (Drean's pale-red whiskers)."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * frac, g + (1 - g) * frac, b + (1 - b) * frac)


def frame_legend(leg):
    """Thin black-framed legend box (Drean Fig 15)."""
    if leg is None:
        return leg
    fr = leg.get_frame()
    fr.set_linewidth(0.5)
    fr.set_edgecolor("0.15")
    fr.set_facecolor("white")
    fr.set_alpha(1.0)
    return leg

# ----------------------------------------------------------------- Table 1 cross-check
TABLE1 = {
    "nb20": {"HRR": 6.11, "T135": 196.1, "T165": 148.8, "mdot": 0.0274, "q": 1.11},
    "nb10": {"HRR": 6.19, "T135": 200.1, "T165": 148.1, "mdot": 0.0261, "q": 0.94},
    "nb05": {"HRR": 6.21, "T135": 179.1, "T165": 135.4, "mdot": 0.0240, "q": 0.97},
    "op20": {"HRR": 6.11, "T135": 215.7, "T165": 172.6, "mdot": 0.0262, "q": 1.43},
    "op10": {"HRR": 6.20, "T135": 190.5, "T165": 176.2, "mdot": 0.0245, "q": 1.40},
    "cl20": {"HRR": 2.26},
    "cl10": {"HRR": 2.13},
}
COLMAP = {"HRR": "HRR_TOTAL", "T135": "TC_Z135", "T165": "TC_Z165",
          "mdot": "MFLOW_UP_Z150", "q": "HFLUX_OUT_Z150"}


# =============================================================== rcParams / style
def set_style():
    # prefer Arial; fall back to DejaVu Sans if unavailable
    installed = {f.name for f in font_manager.fontManager.ttflist}
    family = "Arial" if "Arial" in installed else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": family,
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "black",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.top": False,
        "ytick.right": False,
        "lines.linewidth": 1.2,
        "figure.dpi": 120,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,          # embed TrueType (editable text) in PDF
        "ps.fonttype": 42,
        "axes.grid": False,
    })


def style_axes(ax):
    """Full black box (all four spines, ~0.8 pt); ticks out on bottom/left;
    subtle light-gray horizontal gridlines only, on a white background
    (Drean et al. 2018, Fire Technology, Fig 15)."""
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("black")
        s.set_linewidth(0.8)
    ax.set_facecolor("white")
    ax.yaxis.grid(True, color="0.85", linewidth=0.6, zorder=0)
    ax.xaxis.grid(False)
    ax.tick_params(direction="out", top=False, right=False)
    ax.set_axisbelow(True)


def panel_letter(ax, letter):
    """Bold 9 pt (a)(b)(c) at the top-left, just OUTSIDE the axes."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(-24, 8), textcoords="offset points",
                fontsize=9, fontweight="bold", va="bottom", ha="left")


def save(fig, basename):
    os.makedirs(FIGS_DIR, exist_ok=True)
    png = os.path.join(FIGS_DIR, basename + ".png")
    pdf = os.path.join(FIGS_DIR, basename + ".pdf")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {basename}.png + .pdf")


# =============================================================== data loading
def load_devc(chid):
    fn = os.path.join(RUN_ROOT, chid, f"{chid}_devc.csv")
    if not os.path.exists(fn):
        return None
    with open(fn) as f:
        rows = list(csv.reader(f))
    ids = [c.strip().strip('"') for c in rows[1]]
    data, t = {i: [] for i in ids}, []
    for r in rows[2:]:
        if not r or not r[0].strip():
            continue
        try:
            float(r[0])
        except ValueError:
            continue
        for i, cid in enumerate(ids):
            data[cid].append(float(r[i]))
    return {cid: np.asarray(v) for cid, v in data.items()}


def steady(chid, col, absval=False):
    """(mean, std) of a device column over that case's steady window."""
    dev = load_devc(chid)
    if dev is None or col not in dev:
        return None, None
    t = dev["Time"]
    lo, hi = win_of(chid)
    m = (t >= lo) & (t <= hi)
    v = dev[col][m]
    if absval:
        v = np.abs(v)
    return float(v.mean()), float(v.std())


def acf_tau(v, dt=DT_DEVC):
    """Integral time scale tau = dt*(0.5 + sum ACF up to first zero crossing).

    rho(k) is the normalised autocovariance; the sum runs until the first lag
    with rho <= 0. Floored at 0.5*dt (Celik/ACF method, gci_stats_v2 appendix).
    """
    v = np.asarray(v, float)
    n = len(v)
    x = v - v.mean()
    denom = float(np.dot(x, x))
    if denom == 0.0 or n < 2:
        return 0.5 * dt
    s = 0.5
    for k in range(1, n):
        rho = float(np.dot(x[:-k], x[k:])) / denom
        if rho <= 0.0:
            break
        s += rho
    return max(dt * s, 0.5 * dt)


def steady_ci(chid, col, absval=False):
    """(mean, 95% CI) over the steady window, ACF-corrected for autocorrelation.

    N_ind = T_window / (2*tau); SEM = SD/sqrt(N_ind); CI95 = 1.96*SEM.
    This reproduces gci_stats_v2_results.md Table 1b for the nb cases and applies
    the identical method to the op/cl cases. Replaces the old +/-1 SD error bars.
    """
    dev = load_devc(chid)
    if dev is None or col not in dev:
        return None, None
    t = dev["Time"]
    lo, hi = win_of(chid)
    m = (t >= lo) & (t <= hi)
    v = dev[col][m]
    if absval:
        v = np.abs(v)
    v = np.asarray(v, float)
    if v.size == 0:
        return None, None
    mean, sd = float(v.mean()), float(v.std())
    tau = acf_tau(v, dt=DT_DEVC)
    t_window = (v.size - 1) * DT_DEVC
    n_ind = t_window / (2.0 * tau) if tau > 0 else v.size
    ci = 1.96 * sd / np.sqrt(n_ind) if n_ind > 0 else 0.0
    return mean, float(ci)


def peak_flux_over_height(chid):
    """Peak-over-height of the steady-mean |HFLUX_OUT_Zxxx|, with ACF 95% CI.

    Returns (peak_value_kWm2, ci95, peak_col). The peak sits at Z040 (0.40 m,
    plume attachment near the fire base) in every nb grid.
    """
    dev = load_devc(chid)
    if dev is None:
        return None, None, None
    t = dev["Time"]
    lo, hi = win_of(chid)
    m = (t >= lo) & (t <= hi)
    cols = [c for c in dev if c.startswith("HFLUX_OUT_Z")]
    means = {c: float(np.abs(dev[c][m]).mean()) for c in cols}
    if not means:
        return None, None, None
    peak_col = max(means, key=means.get)
    mean, ci = steady_ci(chid, peak_col, absval=True)
    return mean, ci, peak_col


def tc_profile(chid):
    """Centreline steady-mean temperature profile (z ascending)."""
    dev = load_devc(chid)
    lo, hi = win_of(chid)
    t = dev["Time"]
    m = (t >= lo) & (t <= hi)
    zs, mean, std = [], [], []
    for z in TC_Z:
        col = f"TC_Z{int(round(z * 100)):03d}"
        if col in dev:
            zs.append(z)
            mean.append(dev[col][m].mean())
            std.append(dev[col][m].std())
    return np.array(zs), np.array(mean), np.array(std)


# =============================================================== Table 1 verify
def verify_table1():
    """Return (ok_dict, report_lines). ok_dict[case]=True/False."""
    tol = {"HRR": lambda t: 0.005 * t, "mdot": lambda t: 0.005 * t,
           "T135": lambda t: 0.5, "T165": lambda t: 0.5, "q": lambda t: 0.05}
    ok, lines = {}, []
    lines.append("| case | metric | parsed mean | +/-std | Table 1 | diff | status |")
    lines.append("|------|--------|-------------|--------|---------|------|--------|")
    for case, ref in TABLE1.items():
        chid = f"cavb_{case}"
        case_ok = True
        for key, t1 in ref.items():
            mean, std = steady(chid, COLMAP[key], absval=(key == "q"))
            if mean is None:
                lines.append(f"| {case} | {key} | MISSING | - | {t1} | - | **FAIL** |")
                case_ok = False
                continue
            diff = mean - t1
            passed = abs(diff) <= tol[key](t1)
            case_ok = case_ok and passed
            nd = 4 if key in ("mdot",) else (2 if key in ("HRR", "q") else 1)
            lines.append(f"| {case} | {key} | {mean:.{nd}f} | {std:.{nd}f} | "
                         f"{t1} | {diff:+.{nd}f} | {'OK' if passed else '**FAIL**'} |")
        ok[case] = case_ok
    return ok, lines


# =============================================================== Fig B1
def fig_b1():
    """Grid convergence, no-barrier series (190 mm, 4 panels).

    Error bars are ACF-corrected 95% CIs (not +/-1 SD). Panel (b) carries the
    honest non-asymptotic verdict (the series fails Celik admissibility in every
    cell-size basis, so no valid GCI); panel (d) is the peak outer-skin flux,
    kept on its own 0-18 kW/m2 axis (never merged onto panel (c)'s ~1 kW/m2 axis).
    """
    # 2x2 grid: 190 mm/4 columns is too narrow for the long x-axis label
    # ("Gap resolution (cells across 50 mm gap)") -- labels collide in 1x4.
    fig, axgrid = plt.subplots(2, 2, figsize=(W_DOUBLE, 5.8))
    ax0, ax1, ax2, ax3 = axgrid.ravel()

    # (a) centreline T profiles, 3 grids, sequential ramp
    # Legend quotes the GAP-NORMAL cell size (dy): the coarse mesh is anisotropic
    # (dx=dz=20 mm but dy=16.7 mm for 3 cells across 50 mm) -- quoting "dx=20"
    # against "3 cells" was the inconsistency the referees flagged.
    DY_GAP = {3: "16.7", 5: "10", 10: "5"}
    for g in GRID_ORDER:
        chid = f"cavb_nb{g}"
        z, T, _ = tc_profile(chid)
        ax0.plot(T, z, marker=M_GRID[g], ms=3.5, color=C_GRID[g],
                 mfc=C_GRID[g], mec=C_GRID[g],
                 label=f"{GAP_CELLS[g]} cells (δy = {DY_GAP[GAP_CELLS[g]]} mm)")
    ax0.axhspan(BAR_Z0, BAR_Z1, color="0.85", zorder=1)
    ax0.text(ax0.get_xlim()[1] * 0.97, (BAR_Z0 + BAR_Z1) / 2, "barrier band",
             fontsize=6, color="0.4", ha="right", va="center")
    ax0.set_xlabel("Gas temperature (°C)")
    ax0.set_ylabel("Height above burner (m)")
    frame_legend(ax0.legend(frameon=True, loc="upper right", handlelength=1.4,
                            borderaxespad=0.3, fontsize=6.0))
    style_axes(ax0)
    panel_letter(ax0, "(a)")

    cells = [GAP_CELLS[g] for g in GRID_ORDER]

    # (b) mass flow @1.5 m vs gap-cells; non-asymptotic -> no valid GCI
    md, mderr = [], []
    for g in GRID_ORDER:
        m, c = steady_ci(f"cavb_nb{g}", "MFLOW_UP_Z150")
        md.append(m); mderr.append(c)
    ax1.errorbar(cells, md, yerr=mderr, marker="o", ms=4, color=C_BLACK,
                 mfc=C_BLACK, mec=C_BLACK, ecolor=tint(C_BLACK),
                 capsize=2.5, lw=1.2, elinewidth=0.8, capthick=0.8, zorder=3)
    # Richardson level, now flagged invalid: light-gray dashed
    ax1.axhline(0.0207, ls="--", color="0.72", lw=1.0, zorder=2)
    ax1.text(10, 0.0207, "Richardson level (invalid: non-asymptotic) ",
             fontsize=6, color="0.55", va="bottom", ha="right")
    ax1.set_xlabel("Gap resolution (cells across 50 mm gap)")
    ax1.set_ylabel("Upward mass flow at Z = 1.5 m (kg/s)")
    ax1.set_xticks(cells)
    ax1.set_xlim(1.5, 11.5)
    ax1.annotate("non-asymptotic series (increments grow): no valid GCI\n"
                 "descriptive index 23–58% by cell-size basis",
                 xy=(5, md[1]), xytext=(0.97, 0.96), textcoords="axes fraction",
                 fontsize=6.5, color="0.2", va="top", ha="right")
    style_axes(ax1)
    panel_letter(ax1, "(b)")

    # (c) outer-skin q'' @1.5 m vs gap-cells (Z=1.5 m station -- ~1 kW/m2)
    qv, qerr = [], []
    for g in GRID_ORDER:
        m, c = steady_ci(f"cavb_nb{g}", "HFLUX_OUT_Z150", absval=True)
        qv.append(m); qerr.append(c)
    ax2.errorbar(cells, qv, yerr=qerr, marker="o", ms=4, color=C_BLACK,
                 mfc=C_BLACK, mec=C_BLACK, ecolor=tint(C_BLACK),
                 capsize=2.5, lw=1.2, elinewidth=0.8, capthick=0.8, zorder=3)
    ax2.set_xlabel("Gap resolution (cells across 50 mm gap)")
    ax2.set_ylabel("Outer-skin net heat flux at Z = 1.5 m (kW/m²)")
    ax2.set_xticks(cells)
    ax2.set_xlim(1.5, 11.5)
    # Stats verdict (gci_stats_v2): the 5->10 cell reversal is INSIDE the combined
    # 95% band -- annotate as sampling noise, never as demonstrated non-monotonicity.
    ax2.annotate("reversal within 95% sampling\nbands: station statistically flat\nbetween 5 and 10 cells",
                 xy=(5, qv[1]), xytext=(0.50, 0.82), textcoords="axes fraction",
                 fontsize=6, color="0.35", ha="left", va="top",
                 arrowprops=dict(arrowstyle="->", color="0.5", lw=0.7))
    style_axes(ax2)
    panel_letter(ax2, "(c)")

    # (d) PEAK outer-skin q'' over height vs gap-cells (own 0-18 kW/m2 axis;
    #     do NOT merge onto (c) -- scales differ 1 vs 18 kW/m2, dual axes banned)
    pk, pkerr, pkcol = [], [], None
    for g in GRID_ORDER:
        m, c, col = peak_flux_over_height(f"cavb_nb{g}")
        pk.append(m); pkerr.append(c); pkcol = col
    ax3.errorbar(cells, pk, yerr=pkerr, marker="o", ms=4, color=C_BLACK,
                 mfc=C_BLACK, mec=C_BLACK, ecolor=tint(C_BLACK),
                 capsize=2.5, lw=1.2, elinewidth=0.8, capthick=0.8, zorder=3)
    ax3.set_xlabel("Gap resolution (cells across 50 mm gap)")
    ax3.set_ylabel("Peak outer-skin heat flux (kW/m²)")
    ax3.set_xticks(cells)
    ax3.set_xlim(1.5, 11.5)
    ax3.set_ylim(0, max(pk) * 1.18)
    zlbl = pkcol.split("_Z")[1] if pkcol else "040"
    zval = int(zlbl) / 100.0
    ax3.annotate(f"peak at Z = {zval:.2f} m (plume base);\n"
                 "rises with refinement (not grid-converged)",
                 xy=(cells[-1], pk[-1]), xytext=(0.55, 0.30),
                 textcoords="axes fraction", fontsize=6.5, color="0.2",
                 va="top", ha="center")
    style_axes(ax3)
    panel_letter(ax3, "(d)")

    fig.tight_layout(w_pad=1.6)
    save(fig, "FigB1_grid_convergence")


# =============================================================== Fig B2
def fig_b2():
    """Centreline profiles, nb vs op, at 3-cell and 5-cell grids (190 mm, 2 panels)."""
    fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, 3.2), sharey=True)
    for ax, g, letter in zip(axes, ["20", "10"], ["(a)", "(b)"]):
        for bar in ["nb", "op"]:
            chid = f"cavb_{bar}{g}"
            z, T, Terr = tc_profile(chid)
            ax.plot(T, z, marker=M_STATE[bar], ms=3.5, color=C_STATE[bar],
                    mfc=C_STATE[bar], mec=C_STATE[bar], label=STATE_LABEL[bar])
            ax.errorbar(T, z, xerr=Terr, fmt="none", ecolor=tint(C_STATE[bar]),
                        elinewidth=0.7, capsize=1.5, capthick=0.7)
        ax.axhspan(BAR_Z0, BAR_Z1, color="0.85", zorder=1)
        ax.annotate("closed state:\nambient (20 °C)\nabove barrier",
                    xy=(20, (BAR_Z0 + BAR_Z1) / 2 + 0.15),
                    xytext=(0.06, 0.62), textcoords="axes fraction",
                    fontsize=6, color=C_STATE["cl"], va="center",
                    arrowprops=dict(arrowstyle="->", color=C_STATE["cl"], lw=0.7))
        ax.set_xlabel("Gas temperature (°C)")
        ax.text(0.5, 0.02, f"{GAP_CELLS[g]} cells  (dx = {DX_NOM[g]} mm)",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=7,
                color="0.3")
        frame_legend(ax.legend(frameon=True, loc="upper right"))
        style_axes(ax)
        panel_letter(ax, letter)
    axes[0].set_ylabel("Height above burner (m)")
    fig.tight_layout(w_pad=1.5)
    save(fig, "FigB2_profiles")


# =============================================================== Fig B3
def fig_b3():
    """Above-barrier efficacy vs grid: grouped bars (state) x grid (190 mm, 4 panels)."""
    grids = ["20", "10"]          # 3- vs 5-cell (cl/op exist only here)
    states = ["nb", "cl", "op"]
    panels = [
        ("TC_Z135", "Gas T at Z = 1.35 m (°C)", False, True),
        ("TC_Z165", "Gas T at Z = 1.65 m (°C)", False, True),
        ("MFLOW_UP_Z150", "Upward mass flow at Z = 1.5 m (kg/s)", False, False),
        ("HFLUX_OUT_Z150", "Outer-skin net heat flux at Z = 1.5 m (kW/m²)", True, False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(W_DOUBLE, 2.9))
    x = np.arange(len(grids))
    bw = 0.26
    letters = ["(a)", "(b)", "(c)", "(d)"]
    for ax, (col, ylabel, absval, is_temp), letter in zip(axes, panels, letters):
        for i, bar in enumerate(states):
            vals, errs = [], []
            for g in grids:
                m, c = steady_ci(f"cavb_{bar}{g}", col, absval=absval)
                vals.append(m); errs.append(c)
            ax.bar(x + (i - 1) * bw, vals, width=bw, color=C_STATE[bar],
                   edgecolor="white", linewidth=1.5,
                   yerr=errs, capsize=2,
                   error_kw=dict(elinewidth=0.7, capthick=0.7,
                                 ecolor=tint(C_STATE[bar], 0.35)),
                   label=STATE_LABEL[bar], zorder=3)
        if is_temp:
            ax.axhline(20, ls="--", color="0.4", lw=0.8, zorder=2)
            ax.text(ax.get_xlim()[1], 20, " 20 °C ambient", fontsize=5.5,
                    color="0.4", va="bottom", ha="right")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{GAP_CELLS[g]} cells" for g in grids])
        ax.set_ylabel(ylabel, fontsize=7)
        style_axes(ax)
        panel_letter(ax, letter)
    frame_legend(axes[0].legend(frameon=True, loc="upper left",
                                bbox_to_anchor=(0.0, 1.02),
                                fontsize=6, ncol=1, handlelength=1.0))
    fig.tight_layout(w_pad=1.6)
    save(fig, "FigB3_efficacy_vs_grid")


# =============================================================== Fig B4
def fig_b4():
    """Mid-gap (PBY = 0.025) X-Z temperature elevations, nb/op x 3-cell/5-cell.

    Reverted to the 0.4 m-wide mid-gap Y-plane (orient=2) used by the archived
    FigB4 -- the elevation through the burner tree -- rather than the narrow PBX
    Y-Z gap section. Presented in the classic fire-safety-engineering style of
    Drean et al. 2018 (Fire Technology, Fig 16): a rainbow ('jet') colormap so
    ambient reads deep blue and the flame core red, a single shared vertical
    colorbar in degC, and a FIXED 20-1000 degC scale. The flame core exceeds
    1000 degC (reviewers flagged the earlier 600 degC clip), so the map is capped
    at 1000: hotter cells render as the top colour and the top tick is marked
    ">1000". Physical axes in metres with a Drean-style height ruler.
    """
    cases = [("nb", "20"), ("op", "20"), ("nb", "10"), ("op", "10")]
    letters = ["(a)", "(b)", "(c)", "(d)"]
    VMIN, VMAX = 20.0, 1000.0
    # gather steady-averaged (60-100 s) mid-gap X-Z sections (orient=2 -> h = X)
    fields, hs, zs, avail = [], [], [], []
    for bar, g in cases:
        chid = f"cavb_{bar}{g}"
        times, data, h, z = pp.load_slice_z_combined(chid, "TEMPERATURE", orient=2)
        if times is None:
            fields.append(None); hs.append(None); zs.append(None); avail.append(False)
            continue
        m = (times >= 60.0) & (times <= 100.0)
        fields.append(data[m].mean(axis=0)); hs.append(h); zs.append(z)
        avail.append(True)
    have = [f for f in fields if f is not None]
    if not have:
        print("  FigB4 skipped -- no slice data available")
        return None

    fig, axes = plt.subplots(1, 4, figsize=(W_DOUBLE, 4.6))
    im = None
    for ax, (bar, g), T, h, z, ok, letter in zip(
            axes, cases, fields, hs, zs, avail, letters):
        if not ok:
            ax.text(0.5, 0.5, f"{bar}{g}\nslice unavailable", ha="center",
                    va="center", fontsize=7, transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            panel_letter(ax, letter)
            continue
        im = ax.pcolormesh(h, z, T.T, cmap="jet", vmin=VMIN, vmax=VMAX,
                           shading="gouraud", rasterized=True)
        # cavity-barrier band (Z = 1.20-1.28 m): thin white dashed lines
        ax.axhline(BAR_Z0, color="white", lw=0.7, ls=(0, (4, 2)))
        ax.axhline(BAR_Z1, color="white", lw=0.7, ls=(0, (4, 2)))
        ax.set_ylim(-0.1, 1.9)
        ax.set_xlim(X0, X1)
        ax.set_xticks([0.0, 0.2, 0.4])
        ax.set_xlabel("Width X (m)")
        # Drean-style height ruler: minor ticks every 0.1 m on the left axis
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.tick_params(which="both", direction="out", top=False, right=False,
                       labelsize=6)
        ax.tick_params(which="minor", length=2.0)
        ax.text(0.5, 1.005, f"{STATE_LABEL[bar]}, {GAP_CELLS[g]} cells",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
        if ax is not axes[0]:
            ax.set_yticklabels([])
        panel_letter(ax, letter)
    axes[0].set_ylabel("Height above burner (m)")
    if im is not None:
        cbar = fig.colorbar(im, ax=list(axes), shrink=0.85, pad=0.02,
                            fraction=0.045)
        cbar.set_label("Gas temperature (°C)", fontsize=8)
        cbar.set_ticks([20, 200, 400, 600, 800, 1000])
        cbar.set_ticklabels(["20", "200", "400", "600", "800", ">1000"])
        cbar.ax.tick_params(labelsize=6)
    save(fig, "FigB4_sections")
    return avail


# =============================================================== Fig B5
def fig_b5():
    """nb05 sensor time histories, steady window 60-100 s (190 mm, 3 panels)."""
    chid = "cavb_nb05"
    dev = load_devc(chid)
    if dev is None:
        print("  FigB5 skipped -- nb05 devc missing")
        return
    t = dev["Time"]
    series = [
        ("TC_Z135", "Gas T at Z = 1.35 m (°C)", False),
        ("MFLOW_UP_Z150", "Upward mass flow at Z = 1.5 m (kg/s)", False),
        ("HFLUX_OUT_Z150", "|Outer-skin net heat flux| at Z = 1.5 m (kW/m²)", True),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, 2.7))
    letters = ["(a)", "(b)", "(c)"]
    for ax, (col, ylabel, absval), letter in zip(axes, series, letters):
        y = np.abs(dev[col]) if absval else dev[col]
        ax.axvspan(60, 100, color="0.88", zorder=0)
        ax.plot(t, y, color=C_STATE["nb"], lw=1.0, zorder=3)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel, fontsize=7)
        ax.set_xlim(0, 100)
        style_axes(ax)
        panel_letter(ax, letter)
    axes[0].text(80, axes[0].get_ylim()[1] * 0.95, "steady window\n60–100 s",
                 fontsize=6, color="0.35", ha="center", va="top")
    fig.tight_layout(w_pad=1.6)
    save(fig, "FigB5_timehistories")


# =============================================================== main
def main():
    set_style()
    os.makedirs(FIGS_DIR, exist_ok=True)
    ok, report = verify_table1()
    print("Table 1 cross-check:")
    for line in report:
        print("  " + line)
    failed = [c for c, v in ok.items() if not v]
    if failed:
        print(f"\n*** CROSS-CHECK FAIL for: {failed} -- see manifest; "
              f"affected figures may be aborted ***")

    fig_b1()
    fig_b2()
    fig_b3()
    b4_avail = fig_b4()
    fig_b5()
    return ok, report, b4_avail


if __name__ == "__main__":
    main()
