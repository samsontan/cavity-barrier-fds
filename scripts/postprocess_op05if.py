"""
Paper B -- cavb_op05if postprocessor (the section 4.5 disambiguation run).

Classifies the run outcome from the .out file, then produces every number the
pre-drafted branches in _LOOP/section_4.5_branches.md need:

  BRANCH A (completed stably to T_END=100 s):
      steady means 60-100 s with ACF-corrected 95% CIs (gci_stats_v2 method)
      for T135, T165, mdot150, q150, qout(40..170 cm), HRR;
      two-halves stationarity test; band comparison vs cavb_op10 (middle grid,
      60-120 s) and vs cavb_nb05if (fine-grid no-barrier, same decomposition).

  BRANCH B (numerical instability):
      divergence time, offending mesh, last stable Total Time.

Writes op05if_postprocess_results.md next to this script and prints the same.
Safe to run while the case is still computing (reports RUNNING and exits 2).
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
from gci_stats_v2 import _load_csv, series_stats, two_halves_test  # noqa: E402

RUNS = r"C:\FDS_runs\cavb"
OUT_MD = os.path.join(HERE, "op05if_postprocess_results.md")

# Steady windows: 5 mm cases stop at 100 s, 10 mm baseline used 60-120 s.
CASES = {
    "op05if": (60.0, 100.0),
    "nb05if": (60.0, 100.0),
    "op10":   (60.0, 120.0),
}

DEVC_COLS = {
    "T135":    "TC_Z135",
    "T165":    "TC_Z165",
    "mdot150": "MFLOW_UP_Z150",
    "q150":    "HFLUX_OUT_Z150",
}
HFLUX_OUT_HEIGHTS = [40, 80, 100, 135, 150, 170]


# ------------------------------------------------------------- outcome classify

def classify(chid="cavb_op05if"):
    """Return (state, info) where state is one of
    completed / diverged / stopped / running / missing."""
    fn = os.path.join(RUNS, chid, f"{chid}.out")
    if not os.path.exists(fn):
        return "missing", {}
    txt = open(fn, errors="replace").read()
    times = list(re.finditer(r"Total Time:\s*([0-9.Ee+-]+)", txt))
    t_last = float(times[-1].group(1)) if times else 0.0
    t_pos = times[-1].start() if times else -1
    info = {"t_last": t_last}
    # A restarted case APPENDS to the .out, so stale STOP markers from earlier
    # segments survive -- classify by the LAST marker, and only if nothing
    # (no timestep) was written after it.
    markers = [
        ("completed", txt.rfind("STOP: FDS completed successfully")),
        ("diverged", max(txt.rfind("Numerical Instability"),
                         txt.rfind("STOP: Numerical Instability"))),
        ("stopped", max(txt.rfind("STOP: FDS stopped by user"),
                        txt.rfind("STOP: Set-up only"))),
    ]
    state, pos = max(markers, key=lambda m: m[1])
    if pos < 0 or t_pos > pos:
        return "running", info
    if state == "diverged":
        m = re.search(r"Numerical Instability[^\n]*", txt[pos:])
        info["message"] = m.group(0).strip() if m else "Numerical Instability"
        mesh = re.search(r"Instability[^\n]*Mesh\s+(\d+)", txt[pos:])
        if mesh:
            info["mesh"] = int(mesh.group(1))
    return state, info


# ------------------------------------------------------------------ case loader

def load_case(case):
    """Load devc + hrr for cavb_<case> over its steady window; mirrors
    gci_stats_v2.load_run but keyed by this script's CASES table."""
    win = CASES[case]
    _, dcols = _load_csv(os.path.join(RUNS, f"cavb_{case}", f"cavb_{case}_devc.csv"))
    if dcols is None:
        return None
    t = dcols["Time"]
    mask = (t >= win[0]) & (t <= win[1])
    td = t[mask]
    out = {"window": win}
    for key, col in DEVC_COLS.items():
        out[key] = (dcols[col][mask], td) if col in dcols else None
    for h in HFLUX_OUT_HEIGHTS:
        col = f"HFLUX_OUT_Z{h:03d}"
        out[f"qout_{h}"] = (dcols[col][mask], td) if col in dcols else None
    _, hcols = _load_csv(os.path.join(RUNS, f"cavb_{case}", f"cavb_{case}_hrr.csv"))
    if hcols is not None and "HRR" in hcols:
        th = hcols["Time"]
        hm = (th >= win[0]) & (th <= win[1])
        out["HRR"] = (hcols["HRR"][hm], th[hm])
    else:
        out["HRR"] = None
    return out


def band_compare(sa, sb):
    """|mean_a - mean_b| vs combined 95% band; None if either side missing."""
    if sa is None or sb is None:
        return None
    diff = sa["mean"] - sb["mean"]
    comb = (sa["ci95"] ** 2 + sb["ci95"] ** 2) ** 0.5
    return dict(diff=diff, comb_ci95=comb, within=abs(diff) <= comb)


# ------------------------------------------------------------------------ main

PARAMS = (["T135", "T165", "mdot150", "q150", "HRR"]
          + [f"qout_{h}" for h in HFLUX_OUT_HEIGHTS])

UNITS = {"T135": "C", "T165": "C", "mdot150": "kg/s", "q150": "kW/m2",
         "HRR": "kW", **{f"qout_{h}": "kW/m2" for h in HFLUX_OUT_HEIGHTS}}


def main():
    lines = []
    w = lines.append
    state, info = classify()
    w("# cavb_op05if postprocess results")
    w("")
    w(f"Outcome: **{state.upper()}** (last Total Time {info.get('t_last', 0):.2f} s)")
    if state == "running":
        w("")
        w("Run still computing -- rerun this script after completion.")
        print("\n".join(lines))
        return 2
    if state == "stopped":
        w("")
        w("Graceful stop-file exit (restart dumps written). NOT a crash and NOT")
        w("final -- the logon resume task should pick it up; do not relaunch manually.")
        print("\n".join(lines))
        return 2

    if state == "diverged":
        w("")
        w("## BRANCH B -- diverged with interfaces relocated")
        w("")
        w(f"- Divergence at t = {info['t_last']:.2f} s")
        if "mesh" in info:
            w(f"- Offending mesh: {info['mesh']}")
        w(f"- FDS message: {info.get('message', 'Numerical Instability')}")
        w("")
        w("Fill Branch B of _LOOP/section_4.5_branches.md with the time above.")
    elif state == "completed":
        w("")
        w("## BRANCH A -- ran stably to completion")
        w("")
        op = load_case("op05if")
        op10 = load_case("op10")
        nb = load_case("nb05if")
        stats = {}
        w("| param | mean | 95% CI | two-halves resolved? |")
        w("|---|---|---|---|")
        for p in PARAMS:
            s = series_stats(op[p], op["window"])
            stats[p] = s
            if s is None:
                w(f"| {p} | -- missing -- | | |")
                continue
            th = two_halves_test(op[p], op["window"])
            drift = "yes (DRIFT?)" if (th and th["resolved"]) else "no (stationary)"
            w(f"| {p} [{UNITS[p]}] | {s['mean']:.4g} | +/-{s['ci95']:.3g} | {drift} |")
        w("")
        w("### vs cavb_op10 (middle grid) -- 'within combined sampling bands' test")
        w("")
        w("| param | op05if - op10 | combined 95% band | within? |")
        w("|---|---|---|---|")
        for p in PARAMS:
            c = band_compare(stats[p], series_stats(op10[p], op10["window"]) if op10 else None)
            if c is None:
                w(f"| {p} | -- | -- | -- |")
                continue
            w(f"| {p} | {c['diff']:+.4g} | {c['comb_ci95']:.3g} | "
              f"{'YES' if c['within'] else 'NO'} |")
        w("")
        w("### vs cavb_nb05if (fine-grid no-barrier, same decomposition)")
        w("")
        w("| param | op05if | nb05if | delta |")
        w("|---|---|---|---|")
        for p in PARAMS:
            sn = series_stats(nb[p], nb["window"]) if nb else None
            so = stats[p]
            if so is None or sn is None:
                w(f"| {p} | -- | -- | -- |")
                continue
            w(f"| {p} | {so['mean']:.4g} | {sn['mean']:.4g} | "
              f"{so['mean'] - sn['mean']:+.4g} |")
        w("")
        w("Branch A slot map: HRR -> [X.XX] kW; T165 -> [XXX.X] C (vs nb05if T165);")
        w("mdot150 -> [0.0XXX] kg/s; q150 (and peak qout_40) -> [X.XX] kW/m2.")

    txt = "\n".join(lines)
    with open(OUT_MD, "w") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"\n[written] {OUT_MD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
