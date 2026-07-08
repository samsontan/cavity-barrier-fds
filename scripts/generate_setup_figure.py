"""
Paper B — numerical-setup figure (Fig 1).

Four panels:
  (a) Y-Z section through the cavity: domain extent, leaves, 50 mm gap, base
      burner, barrier band, open/sealed boundaries, and instrumentation.
  (b) 3D isometric of the slab cavity (0.4 m wide x 50 mm deep x 1.8 m tall).
  (c) The three barrier states drawn at gap scale (none / closed / open).
  (d) Gap discretisation at the three grids (3 / 5 / 10 cells) — the paper's thesis:
      a D*/dx-compliant coarse mesh spans the controlling gap with three cells.

Output: fds/figures/FigB0_setup.png (300 dpi).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

FIGS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures"))

# Geometry (m)
GAP = 0.05
ZB0, ZB1 = 0.0, 1.8          # board span
ZD0, ZD1 = -0.3, 2.1         # domain
BAR_Z0, BAR_Z1 = 1.20, 1.28
BURN_ZTOP = 0.06
OPEN_STRIP = 0.02
LEAF_T = 0.012               # drawn thickness (exaggerated for clarity)

C_INNER, C_OUTER = "#8a7a5c", "#9aa0a6"
C_GAP = "#eef3f8"
C_FIRE = "#e8472b"
C_BAR = "#c0379a"
C_OPEN = "#2c7fb8"


def panel_section(ax):
    # gap region
    ax.add_patch(Rectangle((0, ZD0), GAP, ZD1 - ZD0, fc=C_GAP, ec="none", zorder=0))
    # leaves (drawn with exaggerated thickness outside the gap)
    ax.add_patch(Rectangle((-LEAF_T, ZB0), LEAF_T, ZB1 - ZB0, fc=C_INNER, ec="k", lw=0.6, label="Inner leaf (12 mm fibre cement)"))
    ax.add_patch(Rectangle((GAP, ZB0), LEAF_T, ZB1 - ZB0, fc=C_OUTER, ec="k", lw=0.6, label="Outer skin (3 mm aluminium)"))
    # burner
    ax.add_patch(Rectangle((0, ZB0), GAP, BURN_ZTOP, fc=C_FIRE, ec="k", lw=0.5, label="Burner (6.2 kW, 31 kW/m)"))
    # barrier band (shown as closed example, hatched)
    ax.add_patch(Rectangle((0, BAR_Z0), GAP, BAR_Z1 - BAR_Z0, fc=C_BAR, ec="k", lw=0.5, alpha=0.85, label="Cavity barrier (Z=1.20-1.28 m)"))
    # open boundaries: arrows top, bottom, and mouths
    for z, dz in [(ZD1, 0.12), (ZD0, -0.12)]:
        ax.annotate("", xy=(GAP / 2, z + dz), xytext=(GAP / 2, z),
                    arrowprops=dict(arrowstyle="-|>", color=C_OPEN, lw=1.4))
    ax.text(GAP / 2, ZD1 + 0.16, "OPEN", color=C_OPEN, ha="center", fontsize=7, fontweight="bold")
    ax.text(GAP / 2, ZD0 - 0.18, "OPEN (makeup air)", color=C_OPEN, ha="center", fontsize=7, fontweight="bold")
    # instrumentation: TC tree (centreline), flux probes (outer skin), mass-flow planes
    tc_z = [0.2, 0.4, 0.6, 0.8, 1.0, 1.15, 1.35, 1.5, 1.65, 1.75]
    ax.scatter([GAP / 2] * len(tc_z), tc_z, s=9, c="k", marker="o", zorder=5, label="Centreline TC tree")
    hf_z = [0.4, 0.8, 1.0, 1.35, 1.5, 1.7]
    ax.scatter([GAP] * len(hf_z), hf_z, s=18, c="none", edgecolors="#d62728", marker="s", zorder=5, label="Outer-skin heat-flux probes")
    for z in [0.8, 1.5]:
        ax.plot([0, GAP], [z, z], color="#2ca02c", lw=1.0, ls=(0, (4, 2)), zorder=4)
    ax.annotate("ṁ↑ planes\n(0.8, 1.5 m)", xy=(GAP, 0.8), xytext=(GAP + 0.03, 0.45),
                color="#2ca02c", fontsize=6, va="center", ha="left",
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=0.6))
    # dimension annotations
    ax.annotate("", xy=(0, -0.05), xytext=(GAP, -0.05), arrowprops=dict(arrowstyle="<->", lw=0.8))
    ax.text(GAP / 2, -0.10, "50 mm gap", ha="center", fontsize=7)
    ax.annotate("", xy=(-0.03, ZB0), xytext=(-0.03, ZB1), arrowprops=dict(arrowstyle="<->", lw=0.8))
    ax.text(-0.045, 0.9, "1.8 m boards", rotation=90, va="center", fontsize=7)
    ax.text(-0.045, ZD1 - 0.15, "domain\n-0.3 to +2.1 m", fontsize=6, color="0.4", va="top")
    ax.set_xlim(-0.07, GAP + 0.055)
    ax.set_ylim(ZD0 - 0.3, ZD1 + 0.3)
    ax.set_aspect(0.045 / 1.0 * 6)   # stretch x so the thin gap is visible
    ax.set_title("(a) Cavity section (Y-Z)", fontsize=9, fontweight="bold")
    ax.set_ylabel("Height Z (m)", fontsize=8)
    ax.set_xticks([])
    ax.legend(fontsize=5.4, loc="center left", bbox_to_anchor=(1.05, 0.80), framealpha=0.95)


def panel_3d(ax):
    W, D, H = 0.4, GAP, 1.8

    def box(x0, x1, y0, y1, z0, z1, color, alpha):
        v = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
             [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
        faces = [[v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]],
                 [v[0], v[1], v[5], v[4]], [v[2], v[3], v[7], v[6]],
                 [v[1], v[2], v[6], v[5]], [v[0], v[3], v[7], v[4]]]
        ax.add_collection3d(Poly3DCollection(faces, facecolor=color, edgecolor="k", linewidths=0.3, alpha=alpha))

    box(0, W, -0.012, 0, 0, H, C_INNER, 0.85)        # inner leaf
    box(0, W, GAP, GAP + 0.003, 0, H, C_OUTER, 0.5)  # outer skin (translucent)
    box(0, W, 0, GAP, 0, BURN_ZTOP, C_FIRE, 0.95)    # burner
    box(0, W, 0, GAP, BAR_Z0, BAR_Z1, C_BAR, 0.9)    # barrier
    ax.set_xlim(0, W); ax.set_ylim(0, 0.4); ax.set_zlim(0, H)
    ax.set_box_aspect((W, 0.4, H * 0.5))
    ax.view_init(elev=18, azim=-60)
    ax.set_xlabel("X 0.4 m", fontsize=7); ax.set_ylabel("Y", fontsize=7); ax.set_zlabel("Z 1.8 m", fontsize=7)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.set_title("(b) Isometric (gap exaggerated)", fontsize=9, fontweight="bold")


def panel_states(ax):
    # Vertical stack of 3 horizontal gap-sections; labels on the left (no collision).
    # Gap drawn horizontally from GX0 to GX1; inner leaf left, outer skin right.
    GX0, GX1 = 0.38, 0.95
    states = [("None", None), ("Closed", "full"), ("Open", "strip")]
    rows = [0.78, 0.45, 0.12]            # top -> bottom
    h = 0.20
    for (name, kind), y in zip(states, rows):
        ax.add_patch(Rectangle((GX0, y), GX1 - GX0, h, fc=C_GAP, ec="k", lw=0.6))
        ax.add_patch(Rectangle((GX0 - 0.04, y), 0.04, h, fc=C_INNER, ec="k", lw=0.4))   # inner leaf
        ax.add_patch(Rectangle((GX1, y), 0.04, h, fc=C_OUTER, ec="k", lw=0.4))          # outer skin
        if kind == "full":
            ax.add_patch(Rectangle((GX0, y), GX1 - GX0, h, fc=C_BAR, ec="k", lw=0.5))
        elif kind == "strip":
            sx = GX0 + (GX1 - GX0) * (OPEN_STRIP / GAP)   # 20 mm of 50 mm gap
            ax.add_patch(Rectangle((GX0, y), sx - GX0, h, fc=C_BAR, ec="k", lw=0.5))
            ax.annotate("30 mm residual", xy=((sx + GX1) / 2, y + h / 2),
                        xytext=(GX1 + 0.06, y + h / 2), fontsize=5.5, color=C_OPEN,
                        va="center", ha="left",
                        arrowprops=dict(arrowstyle="->", color=C_OPEN, lw=0.6))
        ax.text(0.30, y + h / 2, name, ha="right", va="center", fontsize=7.5, fontweight="bold")
    # inner/outer leaf labels once, on the top row
    ax.text(GX0 - 0.02, rows[0] + h + 0.05, "inner", fontsize=5, ha="center", color="0.4")
    ax.text(GX1 + 0.02, rows[0] + h + 0.05, "outer", fontsize=5, ha="center", color="0.4")
    ax.set_xlim(-0.05, 1.30); ax.set_ylim(0.0, 1.12)
    ax.axis("off")
    ax.set_title("(c) Barrier states", fontsize=9, fontweight="bold")


def panel_grid(ax):
    # Vertical stack of 3 horizontal gap strips with vertical grid lines; labels on the left.
    GX0, GX1 = 0.42, 0.99
    grids = [("δx = 20 mm  (3 cells)", 3), ("δx = 10 mm  (5 cells)", 5), ("δx = 5 mm  (10 cells)", 10)]
    rows = [0.74, 0.44, 0.14]
    h = 0.16
    for (label, n), y in zip(grids, rows):
        ax.add_patch(Rectangle((GX0, y), GX1 - GX0, h, fc="none", ec="k", lw=1.0))
        for e in np.linspace(GX0, GX1, n + 1):
            ax.plot([e, e], [y, y + h], color="0.35", lw=0.6)
        ax.text(0.36, y + h / 2, label, ha="right", va="center", fontsize=6.8)
    ax.text(0.68, 1.02, "Gap resolved by 3 / 5 / 10 cells   (D*/δx = 6.3 / 12.5 / 25)",
            fontsize=6.4, ha="center", color="#b30000")
    ax.set_xlim(-0.05, 1.10); ax.set_ylim(0.05, 1.12)
    ax.axis("off")
    ax.set_title("(d) Gap discretisation", fontsize=9, fontweight="bold")


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig = plt.figure(figsize=(11, 6.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.0, 1.0], height_ratios=[1, 1],
                          hspace=0.35, wspace=0.45)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[:, 1], projection="3d")
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 2])
    panel_section(ax_a)
    panel_3d(ax_b)
    panel_states(ax_c)
    panel_grid(ax_d)
    # No in-figure suptitle -- journals set the title in the caption.
    fn = os.path.join(FIGS, "FigB0_setup.png")
    fig.savefig(fn, dpi=600, bbox_inches="tight")
    fig.savefig(os.path.join(FIGS, "FigB0_setup.pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fn} (+ .pdf)")


if __name__ == "__main__":
    main()
