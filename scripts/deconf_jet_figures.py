"""
Paper B -- candidate figures 2026-07-07 (Sam-approved build, show-first):

  FigB7_deconf_verification   4-panel centreline-profile composite for the
                              de-confounding campaign (one panel per completed
                              4.6.x subsection, baseline vs variant).
  FigB6_constriction_jet      the section-4.5 mechanism: steady W-velocity in
                              the Y-Z gap section (PBX = 0.20) for the open
                              state at both stable grids, plus the vertical
                              jet profile through the residual gap.

Style, colours and the direct .sf reader are imported from the existing
publication modules (FSE-classic per Sam 2026-07-05; never fdsreader).
Data provenance: C:\\FDS_runs\\cavb\\<chid>\\ devc CSVs + slice files.
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import fds_postprocess_cavb as pp  # noqa: E402  (load_devc, slice reader)
from generate_inputs import BAR_Z0, BAR_Z1, TC_Z, Y0, Y1  # noqa: E402
from publication_figures import (  # noqa: E402
    C_BLACK, C_RED, C_BLUE, W_DOUBLE, frame_legend, panel_letter,
    save, set_style, style_axes, tint,
)

# steady windows per case (5 mm cases stop at 100 s)
WIN = {"nb05": (60, 100), "nb05if": (60, 100)}
WIN_DEFAULT = (60.0, 120.0)


def window(chid_short):
    return WIN.get(chid_short, WIN_DEFAULT)


def steady_mean(dev, col, win):
    t = dev["Time"]
    m = (t >= win[0]) & (t <= win[1])
    return float(np.mean(dev[col][m]))


def profile(case):
    """Steady centreline TC profile for cavb_<case>: (z, T)."""
    dev = pp.load_devc(f"cavb_{case}")
    if dev is None:
        return None, None
    win = window(case)
    zs, vals = [], []
    for z in TC_Z:
        col = f"TC_Z{int(z * 100):03d}"
        if col in dev:
            zs.append(z)
            vals.append(steady_mean(dev, col, win))
    return np.array(zs), np.array(vals)


def hrr_mean(case):
    fn = os.path.join(pp.RUN_ROOT, f"cavb_{case}", f"cavb_{case}_hrr.csv")
    _, cols = __import__("gci_stats_v2")._load_csv(fn)
    if cols is None or "HRR" not in cols:
        return None
    t, win = cols["Time"], window(case)
    m = (t >= win[0]) & (t <= win[1])
    return float(np.mean(cols["HRR"][m]))


# =============================================== FigB6: de-confounding composite

def fig_deconf():
    """One panel per completed de-confounding subsection, baseline (black)
    vs variant (red), third reference series (blue) where the text argues
    against one. Profile presentation matches Figure 3 (T on x, Z on y,
    barrier band shaded)."""
    panels = [
        # (letter, title, [(case, colour, marker, label, ls)], annotation)
        ("(a)", "4.6.1 Domain extent",
         [("nb10", C_BLACK, "o", "nb, baseline domain", "-"),
          ("nb10dm", C_RED, "s", "nb, margins doubled", "-")],
         "all metrics within\ncombined 95% bands"),
        ("(b)", "4.6.2 Interface relocation",
         [("nb05", C_BLACK, "o", "nb, original splits", "-"),
          ("nb05if", C_RED, "s", "nb, relocated splits", "-")],
         "T shifts resolved −15 °C;\nmdot, q, HRR\nwithin bands"),
        ("(c)", "4.6.3 Closed-state 5% edge leak",
         [("cl10", C_BLACK, "o", "cl, airtight", "-"),
          ("cl10lk", C_RED, "s", "cl, 5% edge slots", "-"),
          ("nb10", C_BLUE, "^", "no barrier", "--")],
         None),   # per-curve HRR labels added from data below
        ("(d)", "4.6.4 De-snapped strip",
         [("op20", C_BLACK, "o", "op, snapped 16.7 mm strip", "-"),
          ("op20ds", C_RED, "s", "op, exact 20 mm strip", "-"),
          ("op10", C_BLUE, "^", "op, middle grid", "--")],
         None),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(W_DOUBLE, 4.0), sharey=True)
    for ax, (letter, title, series, note) in zip(axes, panels):
        style_axes(ax)
        for case, colr, mk, lab, ls in series:
            z, T = profile(case)
            if z is None:
                ax.text(0.5, 0.5, f"{case} unavailable", transform=ax.transAxes,
                        ha="center", fontsize=7)
                continue
            ax.plot(T, z, marker=mk, ms=3, color=colr, ls=ls, lw=1.1,
                    label=lab, zorder=3)
        ax.axhspan(BAR_Z0, BAR_Z1, color="0.88", zorder=1)
        ax.set_xlabel("Gas temperature (°C)")
        ax.set_ylim(0.1, 1.85)
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.text(0.5, 1.005, title, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=7)
        panel_letter(ax, letter)
        frame_legend(ax.legend(loc="upper right", fontsize=5.6,
                               handlelength=1.6, borderaxespad=0.3))
        if note:
            # bottom-left: the only region all these panels keep clear
            ax.text(0.04, 0.05, note, transform=ax.transAxes, fontsize=5.6,
                    color="0.25", ha="left", va="bottom")
    # panel (b): the resolved -15 degC shift is invisible at full scale --
    # above-barrier zoom inset so the panel shows what its note claims
    # upper-middle-right: curves there hug the left axis and the legend sits
    # above 0.85, so the box overlaps nothing but the barrier-band shading
    axb = axes[1].inset_axes([0.44, 0.50, 0.51, 0.32])
    for case, colr, mk in [("nb05", C_BLACK, "o"), ("nb05if", C_RED, "s")]:
        z, T = profile(case)
        if z is not None:
            m = z >= 1.3
            axb.plot(T[m], z[m], marker=mk, ms=2.6, color=colr, lw=1.0)
    for s in axb.spines.values():
        s.set_linewidth(0.6)
    axb.tick_params(labelsize=5, length=2, direction="out", top=False,
                    right=False)
    axb.set_ylim(1.3, 1.8)
    axb.set_title("above barrier", fontsize=5.6, pad=1.5)

    # panel (c): achieved-HRR label at each curve's foot, computed not
    # hardcoded; z-staggered so the red and blue labels cannot collide
    for case, colr, dz in [("cl10", C_BLACK, 0.02), ("cl10lk", C_RED, 0.02),
                           ("nb10", C_BLUE, 0.105)]:
        h = hrr_mean(case)
        z, T = profile(case)
        if h is None or z is None:
            continue
        i0 = int(np.argmin(z))
        axes[2].text(T[i0] + 25, z[i0] + dz, f"{h:.2f} kW", color=colr,
                     fontsize=5.6, ha="left", va="bottom")
    axes[0].set_ylabel("Height above burner (m)")
    save(fig, "FigB7_deconf_verification")


# =============================================== FigB7: constriction jet

STRIP_Y0, STRIP_Y1 = 0.00, 0.02      # inner-leaf strip (residual gap 30 mm)
ZOOM = (0.90, 1.70)                  # Z window around the barrier
JET_Y = 0.035                        # residual-gap centreline for the profile


def w_section(chid):
    """Steady-averaged W-velocity Y-Z section at PBX=0.20.

    Returns (Wmean[y,z], y, z, Wseries[t,y,z]); the raw windowed series is kept
    so profile panels can carry ACF-corrected sampling bands (slices are stored
    every 2 s -> only ~30 frames in the window; the mean is band-limited)."""
    times, data, y, z = pp.load_slice_z_combined(chid, "W-VELOCITY", orient=1)
    if times is None:
        return None, None, None, None
    win = window(chid.replace("cavb_", ""))
    m = (times >= win[0]) & (times <= win[1])
    if not m.any():
        m = times >= times.max() - 30.0
    return data[m].mean(axis=0), y, z, data[m]


def fig_jet():
    """(a, b) steady W-velocity in the gap section for the open state, coarse
    then middle grid, jet colormap on a fixed scale, strip drawn; (c) vertical
    W profile through the residual-gap centre (Y = 0.035 m) with the no-barrier
    middle grid dashed for reference."""
    cases = [("cavb_op20", "Open state, 3 cells"),
             ("cavb_op10", "Open state, 5 cells")]
    fields = [(chid, lab) + w_section(chid)[:3] for chid, lab in cases]
    ok = [f for f in fields if f[2] is not None]
    if not ok:
        print("  FigB7 skipped -- no W-VELOCITY PBX slices")
        return

    # fixed colour scale across panels, symmetric-ish around the data
    vmax = max(np.nanmax(f[2]) for f in ok)
    vmin = min(np.nanmin(f[2]) for f in ok)
    vmax = np.ceil(vmax * 2) / 2
    vmin = np.floor(vmin * 2) / 2

    fig = plt.figure(figsize=(W_DOUBLE, 3.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.5], wspace=0.42)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    im = None
    for ax, letter, (chid, lab, W, y, z) in zip(axes[:2], ["(a)", "(b)"], fields):
        if W is None:
            ax.text(0.5, 0.5, f"{chid}\nslice unavailable", ha="center",
                    va="center", fontsize=7, transform=ax.transAxes)
            panel_letter(ax, letter)
            continue
        im = ax.pcolormesh(y * 1000.0, z, W.T, cmap="jet", vmin=vmin,
                           vmax=vmax, shading="gouraud", rasterized=True)
        # the strip: solid light-grey block with a black edge
        ax.add_patch(Rectangle((STRIP_Y0 * 1000, BAR_Z0),
                               (STRIP_Y1 - STRIP_Y0) * 1000, BAR_Z1 - BAR_Z0,
                               facecolor="0.82", edgecolor="black", lw=0.6,
                               zorder=4))
        ax.set_xlim(Y0 * 1000, Y1 * 1000)
        ax.set_ylim(*ZOOM)
        ax.set_xticks([0, 10, 20, 30, 40, 50])
        ax.set_xlabel("Gap depth Y (mm)")
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.tick_params(which="both", direction="out", top=False, right=False,
                       labelsize=6)
        ax.tick_params(which="minor", length=2.0)
        ax.text(0.5, 1.005, lab, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=7)
        panel_letter(ax, letter)
    axes[0].set_ylabel("Height above burner (m)")
    axes[1].set_yticklabels([])
    if im is not None:
        cbar = fig.colorbar(im, ax=axes[:2], shrink=0.9, pad=0.03,
                            fraction=0.06)
        cbar.set_label("Vertical velocity W (m/s)", fontsize=8)
        cbar.ax.tick_params(labelsize=6)

    # (c) vertical W profile through the residual-gap centre, with ACF-corrected
    # 95% sampling bands from the windowed slice series (30 frames at 2 s):
    # the below-barrier waviness then reads as what it is -- sampling noise --
    # while the jet-peak separation between grids stays resolved.
    from gci_stats_v2 import _stats_core
    ax = axes[2]
    style_axes(ax)
    series = [("cavb_op20", C_BLUE, "^", "Open, 3 cells", "-"),
              ("cavb_op10", C_RED, "s", "Open, 5 cells", "-"),
              ("cavb_nb10", C_BLACK, None, "No barrier, 5 cells", "--")]
    for chid, colr, mk, lab, ls in series:
        W, y, z, Wt = w_section(chid)
        if W is None:
            continue
        # linear interpolation to EXACTLY Y = 35 mm on every grid (the nearest
        # slice node differs per grid: 33.3 mm coarse, 30/40 mm middle), so all
        # three profiles sample the same physical line
        y = np.asarray(y)
        j1 = int(np.searchsorted(y, JET_Y))
        j0 = j1 - 1
        wgt = (JET_Y - y[j0]) / (y[j1] - y[j0])
        series_t = (1 - wgt) * Wt[:, j0, :] + wgt * Wt[:, j1, :]
        prof = series_t.mean(axis=0)
        tdum = np.arange(Wt.shape[0]) * 2.0          # slice frames are 2 s apart
        ci = np.array([(_stats_core(series_t[:, k], tdum, tdum[-1]) or
                        {"ci95": np.nan})["ci95"] for k in range(len(z))])
        ax.fill_betweenx(z, prof - ci, prof + ci, color=tint(colr, 0.75),
                         lw=0, zorder=2)
        ax.plot(prof, z, color=colr, ls=ls, lw=1.2, marker=mk, ms=2.6,
                markevery=4, label=lab, zorder=3)
    ax.axhspan(BAR_Z0, BAR_Z1, color="0.88", zorder=1)
    ax.axvline(0.0, color="0.75", lw=0.7, zorder=2)
    ax.set_ylim(*ZOOM)
    ax.set_xlabel("W at Y = 35 mm (m/s)")
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.text(0.5, 1.005, "Residual-gap centreline", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=7)
    panel_letter(ax, "(c)")
    frame_legend(ax.legend(loc="upper left", fontsize=6.0, handlelength=1.8,
                           borderaxespad=0.3))
    save(fig, "FigB6_constriction_jet")


if __name__ == "__main__":
    set_style()
    fig_jet()
    fig_deconf()
    print("done")
