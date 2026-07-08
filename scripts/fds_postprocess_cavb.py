"""
Paper B cavity-barrier campaign -- post-processor.

Reads completed runs from C:\\FDS_runs\\cavb\\<chid>\\ and produces the
publication figure set in fds\\figures\\:

  FigB1_grid_convergence.png   centreline dT profile + mass flow + peak flux vs gap/dx
  FigB2_profiles_<grid>.png    barrier-state comparison of centreline profiles per grid
  FigB3_efficacy_vs_grid.png   above-barrier metrics (nb vs cl vs op) at each grid
  FigB4_sections_<grid>.png    mid-gap temperature sections (PBY) nb/cl/op, steady avg
  FigB5_timehistories.png      sensor time histories at the converged grid

Direct Fortran binary .sf reader (fdsreader stale-frame bug bypass) per the
fds-cfd-fire skill; meshes are combined along Z (this campaign's split axis).
Statistics window: 120-180 s.
"""

import csv
import os
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from generate_inputs import GRIDS, X0, X1, Y0, Y1, BAR_Z0, BAR_Z1, TC_Z  # noqa: E402

RUN_ROOT = r"C:\FDS_runs\cavb"
FIGS_DIR = os.path.normpath(os.path.join(HERE, "..", "figures"))
T_AVG = (60.0, 120.0)           # steady-state statistics window

GRID_OF = {"20": 20, "10": 10, "05": 5}
GAP_MM = 50.0
BARRIER_LABEL = {"nb": "No barrier", "cl": "Closed-state barrier", "op": "Open-state barrier"}
BARRIER_COLOR = {"nb": "tab:red", "cl": "tab:blue", "op": "tab:orange"}


# ---------------------------------------------------------------- devc loading

def load_devc(chid):
    """Return dict {column_id: np.array}, plus 'Time'. None if missing."""
    fn = os.path.join(RUN_ROOT, chid, f"{chid}_devc.csv")
    if not os.path.exists(fn):
        return None
    with open(fn) as f:
        rows = list(csv.reader(f))
    ids = [c.strip('"') for c in rows[1]]
    data = np.array([[float(v) for v in r] for r in rows[2:] if r])
    return {cid: data[:, i] for i, cid in enumerate(ids)}


def steady(dev, col):
    """Mean of a device column over the statistics window."""
    t = dev["Time"]
    m = (t >= T_AVG[0]) & (t <= T_AVG[1])
    if not m.any():
        m = t >= t.max() - 30.0
    return float(np.mean(dev[col][m]))


def tc_profile(dev):
    zs, vals = [], []
    for z in TC_Z:
        col = f"TC_Z{int(z*100):03d}"
        if col in dev:
            zs.append(z)
            vals.append(steady(dev, col))
    return np.array(zs), np.array(vals)


# ---------------------------------------------------------------- slice reader

def _rec(f):
    h = f.read(4)
    if len(h) < 4:
        return None
    n = struct.unpack("<i", h)[0]
    d = f.read(n)
    f.read(4)
    return d


def _read_sf_file(filepath):
    with open(filepath, "rb") as f:
        _rec(f); _rec(f); _rec(f)
        idx = _rec(f)
        if idx is None:
            return None, None
        imin, imax, jmin, jmax, kmin, kmax = struct.unpack("<6i", idx)
        ni, nj, nk = imax - imin + 1, jmax - jmin + 1, kmax - kmin + 1
        times, frames = [], []
        while True:
            t_rec = _rec(f)
            if t_rec is None or len(t_rec) < 4:
                break
            t = struct.unpack("<f", t_rec[:4])[0]
            d_rec = _rec(f)
            if d_rec is None:
                break
            arr = np.frombuffer(d_rec, dtype="<f4").reshape((ni, nj, nk), order="F").squeeze()
            times.append(t)
            frames.append(arr)
    if not frames:
        return None, None
    return np.array(times), np.stack(frames, axis=0)


def parse_smv_slices(smv_path):
    mapping = {}
    with open(smv_path) as f:
        lines = [ln.rstrip() for ln in f]
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("SLCF"):
            parts = line.split()
            try:
                mesh_num = int(parts[1])
                amp, bang = line.index("&"), line.index("!")
                imin, imax, jmin, jmax, kmin, kmax = map(int, line[amp + 1:bang].split())
                tail = line[bang + 1:].split()
                slcf_num, orient = int(tail[0]), int(tail[2])
            except (ValueError, IndexError):
                i += 1
                continue
            plane_idx = {3: kmin, 2: jmin, 1: imin}[orient]
            if i + 2 < len(lines):
                key = (lines[i + 2].strip().upper(), orient, plane_idx)
                mapping.setdefault(key, []).append((mesh_num, slcf_num, lines[i + 1].strip()))
            i += 3
        else:
            i += 1
    return mapping


def load_slice_z_combined(chid, quantity, orient):
    """Load a vertical slice (orient 1 or 2) and combine meshes along Z.

    Returns (times, data[n_t, nh, nz_combined], h_nodes, z_nodes) or Nones.
    """
    grid_label = chid[-2:]
    meshes = GRIDS[grid_label][0]
    fds_dir = os.path.join(RUN_ROOT, chid)
    smv = os.path.join(fds_dir, f"{chid}.smv")
    if not os.path.exists(smv):
        return None, None, None, None
    mapping = parse_smv_slices(smv)
    entries = []
    for (qty, o, _pi), ents in mapping.items():
        if qty == quantity.upper() and o == orient:
            entries.extend(ents)
    if not entries:
        return None, None, None, None
    entries.sort(key=lambda e: e[0])         # mesh order = ascending Z

    all_times, blocks, z_blocks = None, [], []
    for k, (mesh_num, _snum, fname) in enumerate(entries):
        fn = os.path.join(fds_dir, fname)
        if not os.path.exists(fn):
            return None, None, None, None
        times, data = _read_sf_file(fn)
        if times is None:
            return None, None, None, None
        I, J, K, z0, z1 = meshes[mesh_num - 1]
        z_nodes = np.linspace(z0, z1, K + 1)
        if k > 0:                              # drop shared boundary plane
            data = data[:, :, 1:]
            z_nodes = z_nodes[1:]
        blocks.append(data)
        z_blocks.append(z_nodes)
        all_times = times
    data = np.concatenate(blocks, axis=2)
    z = np.concatenate(z_blocks)
    I, J, K, _, _ = meshes[0]
    h = np.linspace(X0, X1, I + 1) if orient == 2 else np.linspace(Y0, Y1, J + 1)
    return all_times, data, h, z


def steady_slice(chid, quantity, orient=2):
    times, data, h, z = load_slice_z_combined(chid, quantity, orient)
    if times is None:
        return None, None, None
    m = (times >= T_AVG[0]) & (times <= T_AVG[1])
    if not m.any():
        m = times >= times.max() - 30.0
    return data[m].mean(axis=0), h, z


# ---------------------------------------------------------------- figures

def completed_runs():
    out = []
    for grid in ["20", "10", "05"]:
        for bar in ["nb", "cl", "op"]:
            chid = f"cavb_{bar}{grid}"
            outf = os.path.join(RUN_ROOT, chid, f"{chid}.out")
            if os.path.exists(outf):
                with open(outf, errors="replace") as f:
                    if "completed successfully" in f.read()[-3000:]:
                        out.append(chid)
    return out


def fig_b1(done):
    runs = [c for c in done if c.startswith("cavb_nb")]
    if len(runs) < 2:
        print("  FigB1 skipped (need >=2 nb runs)")
        return
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2))
    gap_dx, mflow, pkflux = [], [], []
    for chid in sorted(runs, key=lambda c: -GRID_OF[c[-2:]]):
        dev = load_devc(chid)
        dx = GRID_OF[chid[-2:]]
        zs, T = tc_profile(dev)
        axes[0].plot(T, zs, marker="o", ms=3,
                     label=f"dx={dx} mm (gap/dx={GAP_MM/dx:.0f}, D*/dx={125/dx:.0f})")
        gap_dx.append(GAP_MM / dx)
        mflow.append(steady(dev, "MFLOW_UP_Z150"))
        flux_cols = [c for c in dev if c.startswith("HFLUX_OUT")]
        pkflux.append(max(abs(steady(dev, c)) for c in flux_cols))
    axes[0].axhspan(BAR_Z0, BAR_Z1, color="0.85")
    axes[0].set_xlabel("Gas temperature (°C)")
    axes[0].set_ylabel("Height (m)")
    axes[0].set_title("(a) Cavity centreline temperature", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[1].plot(gap_dx, mflow, "ks-")
    axes[1].set_xlabel("Cavity gap / dx")
    axes[1].set_ylabel("Upward mass flow at Z=1.5 m (kg/s)")
    axes[1].set_title("(b) Mass flow convergence", fontsize=9)
    axes[2].plot(gap_dx, pkflux, "ks-")
    axes[2].set_xlabel("Cavity gap / dx")
    axes[2].set_ylabel("Peak outer-skin net heat flux (kW/m²)")
    axes[2].set_title("(c) Wall heat flux convergence", fontsize=9)
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.suptitle("Grid convergence at fixed D*-based 'compliance' -- no-barrier cavity, 6.2 kW",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    fn = os.path.join(FIGS_DIR, "FigB1_grid_convergence.png")
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fn}")


def fig_b2(done):
    for grid in ["20", "10", "05"]:
        runs = [f"cavb_{b}{grid}" for b in ["nb", "cl", "op"] if f"cavb_{b}{grid}" in done]
        if len(runs) < 2:
            continue
        fig, ax = plt.subplots(figsize=(5, 5.5))
        for chid in runs:
            bar = chid[5:7]
            zs, T = tc_profile(load_devc(chid))
            ax.plot(T, zs, marker="o", ms=3, color=BARRIER_COLOR[bar],
                    label=BARRIER_LABEL[bar])
        ax.axhspan(BAR_Z0, BAR_Z1, color="0.85", label="Barrier position")
        ax.set_xlabel("Gas temperature (°C)")
        ax.set_ylabel("Height (m)")
        dx = GRID_OF[grid]
        ax.set_title(f"Centreline temperature, dx={dx} mm (gap/dx={GAP_MM/dx:.0f})",
                     fontsize=10, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fn = os.path.join(FIGS_DIR, f"FigB2_profiles_dx{grid}.png")
        fig.savefig(fn, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {fn}")


def fig_b3(done):
    grids = [g for g in ["20", "10", "05"]
             if all(f"cavb_{b}{g}" in done for b in ["nb", "cl", "op"])]
    if not grids:
        print("  FigB3 skipped (no grid has all three states)")
        return
    metrics = {
        "T_Z135": ("TC_Z135", "Gas T at Z=1.35 m (°C)"),
        "T_Z165": ("TC_Z165", "Gas T at Z=1.65 m (°C)"),
        "MFLOW": ("MFLOW_UP_Z150", "Mass flow Z=1.5 m (kg/s)"),
        "FLUX150": ("HFLUX_OUT_Z150", "Outer-skin flux Z=1.5 m (kW/m²)"),
    }
    fig, axes = plt.subplots(1, len(metrics), figsize=(3.0 * len(metrics), 4.0))
    x = np.arange(len(grids))
    for ax, (mk, (col, ylabel)) in zip(axes, metrics.items()):
        for i, bar in enumerate(["nb", "cl", "op"]):
            vals = [steady(load_devc(f"cavb_{bar}{g}"), col) for g in grids]
            if "FLUX" in mk:
                vals = [abs(v) for v in vals]
            ax.bar(x + (i - 1) * 0.26, vals, width=0.24,
                   color=BARRIER_COLOR[bar], label=BARRIER_LABEL[bar])
        ax.set_xticks(x)
        ax.set_xticklabels([f"dx={GRID_OF[g]}\ngap/dx={GAP_MM/GRID_OF[g]:.0f}" for g in grids],
                           fontsize=7)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(alpha=0.3, axis="y")
    axes[0].legend(fontsize=7)
    fig.suptitle("Above-barrier exposure vs grid resolution -- does the coarse grid mis-state efficacy?",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    fn = os.path.join(FIGS_DIR, "FigB3_efficacy_vs_grid.png")
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fn}")


def fig_b4(done):
    for grid in ["20", "10", "05"]:
        runs = [f"cavb_{b}{grid}" for b in ["nb", "cl", "op"] if f"cavb_{b}{grid}" in done]
        if len(runs) < 2:
            continue
        fig, axes = plt.subplots(1, len(runs), figsize=(2.6 * len(runs), 7))
        if len(runs) == 1:
            axes = [axes]
        im = None
        for ax, chid in zip(axes, runs):
            T, h, z = steady_slice(chid, "TEMPERATURE", orient=2)
            if T is None:
                ax.set_title(f"{chid}: slice missing", fontsize=8)
                continue
            im = ax.pcolormesh(h, z, T.T, cmap="inferno", vmin=20, vmax=600, shading="gouraud")
            ax.axhspan(BAR_Z0, BAR_Z1, color="cyan", alpha=0.35)
            ax.set_title(BARRIER_LABEL[chid[5:7]], fontsize=9)
            ax.set_xlabel("X (m)")
            ax.set_ylim(-0.3, 2.1)
        axes[0].set_ylabel("Height (m)")
        if im:
            fig.colorbar(im, ax=axes, label="Temperature (°C)", shrink=0.7)
        dx = GRID_OF[grid]
        fig.suptitle(f"Mid-gap temperature (steady {T_AVG[0]:.0f}-{T_AVG[1]:.0f} s), dx={dx} mm",
                     fontsize=10, fontweight="bold")
        fn = os.path.join(FIGS_DIR, f"FigB4_sections_dx{grid}.png")
        fig.savefig(fn, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {fn}")


def fig_b5(done):
    grid = "05" if any(c.endswith("05") for c in done) else \
           "10" if any(c.endswith("10") for c in done) else "20"
    runs = [f"cavb_{b}{grid}" for b in ["nb", "cl", "op"] if f"cavb_{b}{grid}" in done]
    if not runs:
        return
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for chid in runs:
        dev = load_devc(chid)
        bar = chid[5:7]
        axes[0].plot(dev["Time"], dev["TC_Z135"], color=BARRIER_COLOR[bar],
                     label=BARRIER_LABEL[bar])
        axes[1].plot(dev["Time"], dev["MFLOW_UP_Z150"], color=BARRIER_COLOR[bar])
        axes[2].plot(dev["Time"], np.abs(dev["HFLUX_OUT_Z150"]), color=BARRIER_COLOR[bar])
    axes[0].set_ylabel("Gas T at Z=1.35 m (°C)")
    axes[1].set_ylabel("Upward mass flow Z=1.5 m (kg/s)")
    axes[2].set_ylabel("|Outer-skin flux| Z=1.5 m (kW/m²)")
    for ax in axes:
        ax.set_xlabel("Time (s)")
        ax.grid(alpha=0.3)
        ax.axvspan(*T_AVG, color="0.92")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Sensor time histories above barrier, dx={GRID_OF[grid]} mm",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    fn = os.path.join(FIGS_DIR, "FigB5_timehistories.png")
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {fn}")


def main():
    os.makedirs(FIGS_DIR, exist_ok=True)
    done = completed_runs()
    print(f"Completed runs: {done}")
    if not done:
        print("Nothing to post-process yet.")
        return
    fig_b1(done)
    fig_b2(done)
    fig_b3(done)
    fig_b4(done)
    fig_b5(done)


if __name__ == "__main__":
    main()
