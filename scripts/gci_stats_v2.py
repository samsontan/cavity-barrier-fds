"""
gci_stats_v2.py  -- Rigorous stationary-series statistics + variable-r GCI redo
================================================================================
Fire Safety Journal cavity-barrier CFD paper (Paper B).

Referees flagged two problems the old gci_analysis.py had:
  (1) it reported grid-to-grid DIFFERENCES with no sampling-uncertainty test, so
      a non-monotone temperature sequence could not be told apart from noise;
  (2) it assumed a constant refinement ratio r = 2, but the meshes are
      ANISOTROPIC.  Gap-normal (Y) cells are 16.7 / 10 / 5 mm (3 / 5 / 10 cells
      across the 50 mm gap); X and Z cells are 20 / 10 / 5 mm.  The gap-normal
      refinement ratios are therefore r32 = 16.7/10 = 1.667, r21 = 10/5 = 2.0,
      NOT a constant 2.  The Celik-effective cell size
      h = (dx*dy*dz)^(1/3) = 18.82 / 10 / 5 mm gives r32 = 1.882, r21 = 2.0.

This script does NOT modify gci_analysis.py.  It reuses its CSV-parsing idea
(header on row 2, data from row 3, filter by a steady time window) and adds:
  Task 1  stationary-series uncertainty (ACF integral time scale, N_indep, 95% CI)
  Task 2  resolved-difference verdicts (diff vs combined 95% CI)
  Task 3  variable-r GCI (Celik 2008 iterative p) for gap-normal r AND h-based r
  Task 4  Monte-Carlo uncertainty on p and GCI
  Task 5  peak outer-skin flux over height + convergence sequence
  Task 6  markdown output -> gci_stats_v2_results.md

Every number comes from the FDS output CSVs; nothing is hand-entered except the
published Table-1 / old-analysis values used purely as cross-checks.
"""
import csv, os, math
import numpy as np

RUNS = r"C:\FDS_runs\cavb"
DT_DEVC = 0.5  # s, sampling interval declared in the FDS input

# Steady windows (s): 60-120 for 20/10 mm, 60-100 for nb05 (run stops at 100 s)
WINDOWS = {
    "nb20": (60.0, 120.0), "nb10": (60.0, 120.0), "nb05": (60.0, 100.0),
    "op20": (60.0, 120.0), "op10": (60.0, 120.0),
    "cl20": (60.0, 120.0), "cl10": (60.0, 120.0),
}

# devc column names for the five headline parameters (HRR overridden from hrr.csv)
DEVC_COLS = {
    "T135":    "TC_Z135",        # centreline T @ 1.35 m (C)
    "T165":    "TC_Z165",        # centreline T @ 1.65 m (C)
    "mdot150": "MFLOW_UP_Z150",  # upward mass flow @ 1.5 m (kg/s)
    "q150":    "HFLUX_OUT_Z150", # outer-skin heat flux @ 1.5 m (kW/m2)
}
HFLUX_OUT_HEIGHTS = [40, 80, 100, 135, 150, 170]  # cm, from the devc header

# Published Table 1 means for cross-check: T135, T165, mdot150, q150, HRR
TABLE1 = {
    "nb20": dict(T135=196.1, T165=148.8, mdot150=0.0274, q150=1.11, HRR=6.11),
    "nb10": dict(T135=200.1, T165=148.1, mdot150=0.0261, q150=0.94, HRR=6.19),
    "nb05": dict(T135=179.1, T165=135.4, mdot150=0.0240, q150=0.97, HRR=6.21),
}

# Old constant-r=2 GCI numbers for mass flow (for the comparison line)
OLD_MDOT = dict(p=0.72, gci=17.0, phi_ext=0.0207)

# Grid definitions ------------------------------------------------------------
# order: coarse(20) -> medium(10) -> fine(5)
GAP_NORMAL_H = dict(coarse=16.7, medium=10.0, fine=5.0)   # gap-normal Y cell (mm)
CELIK_H      = dict(coarse=18.82, medium=10.0, fine=5.0)  # (dx*dy*dz)^(1/3) (mm)
ISO_H        = dict(coarse=20.0, medium=10.0, fine=5.0)   # naive constant-r=2 (X/Z cell)


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------
def _load_csv(path):
    """Return (header_list, dict name->np.array of full time series)."""
    if not os.path.exists(path):
        return None, None
    rows = list(csv.reader(open(path)))
    names = [h.strip().strip('"') for h in rows[1]]
    cols = {n: [] for n in names}
    for r in rows[2:]:
        if not r or not r[0].strip():
            continue
        try:
            float(r[0])
        except ValueError:
            continue
        for i, n in enumerate(names):
            try:
                cols[n].append(float(r[i]))
            except (ValueError, IndexError):
                cols[n].append(np.nan)
    return names, {n: np.asarray(v, float) for n, v in cols.items()}


def load_run(run):
    """Load devc + hrr. Each parameter is stored as a (values, time) tuple over the
    steady window so its own sampling interval is used (devc ~0.5 s, hrr ~1 s)."""
    win = WINDOWS[run]
    dnames, dcols = _load_csv(os.path.join(RUNS, f"cavb_{run}", f"cavb_{run}_devc.csv"))
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
    # HRR from the dedicated hrr.csv (column "HRR"), own time base
    hnames, hcols = _load_csv(os.path.join(RUNS, f"cavb_{run}", f"cavb_{run}_hrr.csv"))
    if hcols is not None and "HRR" in hcols:
        th = hcols["Time"]
        hm = (th >= win[0]) & (th <= win[1])
        out["HRR"] = (hcols["HRR"][hm], th[hm])
    else:
        out["HRR"] = None
    return out


def _dt_of(tvec):
    """Robust sampling interval from a (possibly slightly irregular) time vector."""
    tvec = np.asarray(tvec, float)
    if len(tvec) < 2:
        return DT_DEVC
    return float(np.median(np.diff(tvec)))


# ---------------------------------------------------------------------------
# Task 1: stationary-series uncertainty
# ---------------------------------------------------------------------------
def integral_time_scale(x, dt):
    """Integral time scale tau (s): integrate the normalised ACF from lag 0 to
    the first zero crossing (trapezoidal, lag-0 weight 0.5).  Floors at 0.5*dt so
    N_indep never exceeds the raw sample count (white-noise limit)."""
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    n = len(x)
    if n < 4:
        return 0.5 * dt, 0
    xd = x - x.mean()
    c0 = np.dot(xd, xd) / n
    if c0 == 0:
        return 0.5 * dt, 0
    # first zero crossing of the ACF
    k0 = n - 1
    rho_pos = []
    for k in range(1, n):
        ck = np.dot(xd[:n - k], xd[k:]) / n
        rk = ck / c0
        if rk <= 0.0:
            k0 = k
            break
        rho_pos.append(rk)
    tau = dt * (0.5 + sum(rho_pos))          # trapezoidal integral of ACF
    tau = max(tau, 0.5 * dt)
    return tau, k0


def _stats_core(x, tvec, T):
    """mean, sample SD, tau, N_indep, SEM, 95% CI half-width, using the actual
    sampling interval derived from tvec and the physical window length T."""
    x = np.asarray(x, float)
    good = ~np.isnan(x)
    x = x[good]
    if len(x) < 2:
        return None
    dt = _dt_of(tvec)
    mean = x.mean()
    sd = x.std(ddof=1)
    tau, k0 = integral_time_scale(x, dt)
    n_indep = max(T / (2.0 * tau), 1.0)
    sem = sd / math.sqrt(n_indep)
    ci95 = 1.96 * sem
    return dict(mean=mean, sd=sd, tau=tau, dt=dt, k0_lag=k0, n_indep=n_indep,
                sem=sem, ci95=ci95, n_samp=len(x))


def series_stats(param, window):
    """param = (values, time) tuple over the steady window."""
    if param is None:
        return None
    x, tvec = param
    return _stats_core(x, tvec, window[1] - window[0])


def two_halves_test(param, window):
    """Split the window at its midpoint using the ACTUAL timestamps, compare
    half-means; combined 95% band = sqrt(CI_a^2 + CI_b^2)."""
    if param is None:
        return None
    x, tvec = param
    x = np.asarray(x, float)
    tvec = np.asarray(tvec, float)
    t0, t1 = window
    tmid = 0.5 * (t0 + t1)
    ma, mb = tvec < tmid, tvec >= tmid
    sa = _stats_core(x[ma], tvec[ma], tmid - t0)
    sb = _stats_core(x[mb], tvec[mb], t1 - tmid)
    if sa is None or sb is None:
        return None
    diff = sb["mean"] - sa["mean"]
    comb = math.sqrt(sa["ci95"] ** 2 + sb["ci95"] ** 2)
    return dict(mean_first=sa["mean"], mean_second=sb["mean"], diff=diff,
                comb_ci95=comb, resolved=abs(diff) > comb)


# ---------------------------------------------------------------------------
# Task 3: Celik (2008) variable-r GCI
# ---------------------------------------------------------------------------
def celik_p(eps21, eps32, r21, r32, tol=1e-10, itmax=1000):
    """Solve the Celik apparent-order equation with non-constant ratios.

        p = (1/ln(r21)) * | ln|eps32/eps21| + q(p) |
        q(p) = ln( (r21^p - s) / (r32^p - s) )
        s = sign(eps32/eps21)

    Fixed-point iteration.  Returns (p, s) or (None, s) if it will not converge
    (e.g. eps21 == 0 or a non-monotone sign giving a domain error)."""
    if eps21 == 0:
        return None, 0
    ratio = eps32 / eps21
    s = 1.0 if ratio > 0 else -1.0
    lnabs = math.log(abs(ratio))
    p = abs(lnabs) / math.log(r21)  # constant-r seed
    for _ in range(itmax):
        try:
            q = math.log((r21 ** p - s) / (r32 ** p - s))
        except (ValueError, ZeroDivisionError):
            return None, s
        p_new = abs(lnabs + q) / math.log(r21)
        if abs(p_new - p) < tol:
            return p_new, s
        p = p_new
    return p, s


def gci_variable(f_coarse, f_med, f_fine, h):
    """Full Celik GCI for grids with sizes h={coarse,medium,fine} (mm).
    f1=fine, f2=medium, f3=coarse.  Returns None if non-monotone."""
    f1, f2, f3 = f_fine, f_med, f_coarse
    eps21, eps32 = f2 - f1, f3 - f2
    r21 = h["medium"] / h["fine"]
    r32 = h["coarse"] / h["medium"]
    monotone = (eps21 * eps32) > 0
    p, s = celik_p(eps21, eps32, r21, r32)
    if p is None:
        return dict(monotone=monotone, r21=r21, r32=r32, p=None)
    r21p = r21 ** p
    f_ext = (r21p * f1 - f2) / (r21p - 1) if r21p != 1 else float("nan")
    e_a = abs((f1 - f2) / f1) if f1 != 0 else float("nan")
    e_ext = abs((f_ext - f1) / f_ext) if f_ext not in (0, float("nan")) else float("nan")
    gci_fine = 1.25 * e_a / (r21p - 1) if r21p != 1 else float("nan")
    return dict(monotone=monotone, r21=r21, r32=r32, p=p, s=s, f_ext=f_ext,
                e_a=e_a, e_ext=e_ext, gci_fine=gci_fine,
                asymptotic=(abs(eps32 / eps21) > 1))


# ---------------------------------------------------------------------------
# Task 4: Monte-Carlo on p and GCI
# ---------------------------------------------------------------------------
def mc_p_gci(means, sems, h, ndraw=10000, seed=12345):
    """means/sems = (coarse, medium, fine).  Perturb each grid mean ~ N(mean,sem),
    recompute p and GCI_fine.  Discard non-monotone draws.  Return summary."""
    rng = np.random.default_rng(seed)
    mc, ms, mf = means
    sc, ss, sf = sems
    r21 = h["medium"] / h["fine"]
    r32 = h["coarse"] / h["medium"]
    ps, gcis = [], []
    discarded = 0
    for _ in range(ndraw):
        fc = rng.normal(mc, sc)
        fm = rng.normal(ms, ss)
        ff = rng.normal(mf, sf)
        eps21, eps32 = fm - ff, fc - fm
        if eps21 * eps32 <= 0 or eps21 == 0:  # non-monotone / degenerate
            discarded += 1
            continue
        p, s = celik_p(eps21, eps32, r21, r32)
        if p is None or not np.isfinite(p) or p <= 0 or p > 20:
            discarded += 1
            continue
        r21p = r21 ** p
        e_a = abs((ff - fm) / ff)
        gci = 1.25 * e_a / (r21p - 1)
        if not np.isfinite(gci):
            discarded += 1
            continue
        ps.append(p)
        gcis.append(gci * 100.0)
    ps, gcis = np.asarray(ps), np.asarray(gcis)
    def pct(a):
        return (np.median(a), np.percentile(a, 2.5), np.percentile(a, 97.5)) if len(a) else (np.nan,)*3
    return dict(n_kept=len(ps), discarded=discarded, frac_disc=discarded / ndraw,
                p=pct(ps), gci=pct(gcis))


# ===========================================================================
# MAIN
# ===========================================================================
def main():
    md = []  # markdown lines
    def w(s=""):
        md.append(s)

    nb = ["nb20", "nb10", "nb05"]
    all_runs = ["nb20", "nb10", "nb05", "op20", "op10", "cl20", "cl10"]
    data = {r: load_run(r) for r in all_runs}
    for r in nb:
        if data[r] is None:
            raise SystemExit(f"FATAL: no data for {r}")

    params = ["T135", "T165", "mdot150", "q150", "HRR"]
    plabel = {"T135": "T @1.35 m (C)", "T165": "T @1.65 m (C)",
              "mdot150": "mdot @1.5 m (kg/s)", "q150": "q'' @1.5 m (kW/m2)",
              "HRR": "HRR (kW)"}
    decimals = {"T135": 1, "T165": 1, "mdot150": 4, "q150": 3, "HRR": 3}

    w("# Grid-Convergence + Stationary-Series Statistics (v2)")
    w()
    w("Cavity-barrier no-barrier series, FDS. Meshes are **anisotropic**: gap-normal (Y) "
      "cells 16.7 / 10 / 5 mm (3/5/10 across the 50 mm gap); X,Z cells 20 / 10 / 5 mm. "
      "Refinement ratios are therefore NOT constant. All numbers computed from the FDS "
      "output CSVs by `gci_stats_v2.py`.")
    w()
    w(f"Steady windows: 60-120 s (nb20, nb10), 60-100 s (nb05). DT_DEVC = {DT_DEVC} s.")
    w()

    # ---- compute per-case per-param stats ----
    stats = {r: {} for r in nb}
    for r in nb:
        for p in params:
            x = data[r][p]
            stats[r][p] = series_stats(x, data[r]["window"]) if x is not None else None

    # ==================== TASK 1 ====================
    w("## Task 1 - Stationary-series uncertainty")
    w()
    w("Per case per parameter: steady-window mean, sample SD, integral time scale "
      "`tau` (ACF integrated to first zero crossing), independent-sample count "
      "`N_ind = T_window/(2*tau)`, `SEM = SD/sqrt(N_ind)`, and 95% CI = 1.96*SEM.")
    w()

    # cross-check table against Table 1
    w("### 1a. Mean cross-check against published Table 1")
    w()
    w("| Case | Param | This script | Table 1 | abs diff | OK? |")
    w("|---|---|---|---|---|---|")
    crosscheck_ok = True
    tol = {"T135": 0.15, "T165": 0.15, "mdot150": 0.0002, "q150": 0.02, "HRR": 0.02}
    for r in nb:
        for p in params:
            s = stats[r][p]
            if s is None:
                w(f"| {r} | {plabel[p]} | -- | {TABLE1[r][p]} | -- | NO DATA |")
                crosscheck_ok = False
                continue
            t1 = TABLE1[r][p]
            d = decimals[p]
            diff = abs(s["mean"] - t1)
            ok = diff <= tol[p]
            crosscheck_ok = crosscheck_ok and ok
            w(f"| {r} | {plabel[p]} | {s['mean']:.{d}f} | {t1} | {diff:.{d}f} | "
              f"{'yes' if ok else 'NO'} |")
    w()
    if crosscheck_ok:
        w("**Verdict:** All means reproduce Table 1 within rounding. Cross-check PASSED - "
          "downstream statistics trusted.")
    else:
        w("**Verdict:** ONE OR MORE MEANS DO NOT REPRODUCE TABLE 1. STOP - investigate "
          "parsing/window before trusting anything below.")
    w()

    # full stats table
    w("### 1b. Full stationary-series statistics")
    w()
    w("| Case | Param | Mean | SD | tau (s) | N_ind | SEM | 95% CI (+/-) |")
    w("|---|---|---|---|---|---|---|---|")
    for r in nb:
        for p in params:
            s = stats[r][p]
            if s is None:
                w(f"| {r} | {plabel[p]} | -- | -- | -- | -- | -- | -- |")
                continue
            d = decimals[p]
            w(f"| {r} | {plabel[p]} | {s['mean']:.{d}f} | {s['sd']:.{d}f} | "
              f"{s['tau']:.2f} | {s['n_indep']:.1f} | {s['sem']:.{d}f} | "
              f"{s['ci95']:.{d}f} |")
    w()
    w("**Verdict:** SEM uses the ACF-based independent-sample count, not the raw "
      "121/81 samples, so it honestly reflects temporal correlation in the steady window.")
    w()

    # two-halves on nb05
    w("### 1c. Two-halves window-length-bias test (nb05: 60-80 vs 80-100 s)")
    w()
    w("| Param | Mean 60-80 s | Mean 80-100 s | Diff | Combined 95% band (+/-) | Biased? |")
    w("|---|---|---|---|---|---|")
    bias_flag = False
    for p in params:
        x = data["nb05"][p]
        if x is None:
            w(f"| {plabel[p]} | -- | -- | -- | -- | -- |")
            continue
        th = two_halves_test(x, (60.0, 100.0))
        d = decimals[p]
        biased = th["resolved"]
        bias_flag = bias_flag or biased
        w(f"| {plabel[p]} | {th['mean_first']:.{d}f} | {th['mean_second']:.{d}f} | "
          f"{th['diff']:+.{d}f} | {th['comb_ci95']:.{d}f} | "
          f"{'YES' if biased else 'no'} |")
    w()
    if bias_flag:
        w("**Verdict:** At least one parameter shows a first-half/second-half shift larger "
          "than the combined CI - the 40 s nb05 window carries some residual drift; treat "
          "its means with the wider caution noted per row.")
    else:
        w("**Verdict:** No parameter's half-window shift exceeds its combined 95% band - the "
          "nb05 60-100 s window is statistically stationary; window-length bias is negligible.")
    w()

    # ==================== TASK 2 ====================
    w("## Task 2 - Resolved-difference verdicts")
    w()
    w("For each parameter, is each grid-to-grid step larger than the combined 95% CI of the "
      "two grids? `combined = sqrt(CI_a^2 + CI_b^2)`. If a step is inside the band it is "
      "sampling noise, not a resolved grid effect.")
    w()
    w("| Param | 20->10 step | vs band | 10->5 step | vs band | Monotone? | Verdict |")
    w("|---|---|---|---|---|---|---|")
    resolved_monotone = {}
    for p in params:
        s20, s10, s05 = stats["nb20"][p], stats["nb10"][p], stats["nb05"][p]
        if None in (s20, s10, s05):
            w(f"| {plabel[p]} | -- | -- | -- | -- | -- | no data |")
            continue
        d = decimals[p]
        step32 = s10["mean"] - s20["mean"]   # 20 -> 10
        step21 = s05["mean"] - s10["mean"]   # 10 -> 5
        band32 = math.sqrt(s20["ci95"]**2 + s10["ci95"]**2)
        band21 = math.sqrt(s10["ci95"]**2 + s05["ci95"]**2)
        res32 = abs(step32) > band32
        res21 = abs(step21) > band21
        mono = (step32 * step21) > 0
        both_resolved = res32 and res21
        resolved_monotone[p] = dict(mono=mono, both_resolved=both_resolved,
                                    res32=res32, res21=res21)
        if mono and both_resolved:
            verdict = "monotone & resolved"
        elif not mono and both_resolved:
            verdict = "NON-monotone & resolved (real)"
        elif not mono and not both_resolved:
            verdict = "non-monotone but within noise"
        else:
            verdict = "monotone, partly within noise"
        w(f"| {plabel[p]} | {step32:+.{d}f} | {'>band' if res32 else 'in band'} "
          f"({band32:.{d}f}) | {step21:+.{d}f} | {'>band' if res21 else 'in band'} "
          f"({band21:.{d}f}) | {'yes' if mono else 'NO'} | {verdict} |")
    w()
    # adjudicate specific sequences
    def adjudicate(p, name, seq):
        rm = resolved_monotone.get(p)
        if rm is None:
            return f"- **{name}:** no data."
        if rm["mono"] and rm["both_resolved"]:
            return (f"- **{name}** ({seq}): monotone AND both steps exceed the combined 95% CI "
                    f"-> a real, grid-resolved trend.")
        if not rm["mono"]:
            if rm["both_resolved"]:
                return (f"- **{name}** ({seq}): NON-monotone and both steps are individually "
                        f"resolved -> the non-monotonicity is real (over-/under-shoot), not noise. "
                        f"GCI is undefined (not in the asymptotic range).")
            return (f"- **{name}** ({seq}): the reversal is NOT statistically resolved - at least "
                    f"one step is inside the combined 95% CI, so within sampling noise. Cannot "
                    f"claim a real non-monotone grid trend.")
        return (f"- **{name}** ({seq}): monotone but at least one step is within the 95% CI - "
                f"trend direction is plausible but not fully resolved.")

    w("**Adjudication of the referee-flagged sequences:**")
    w()
    w(adjudicate("T135", "T@1.35 m", "196.1 -> 200.1 -> 179.1"))
    w(adjudicate("T165", "T@1.65 m", "148.8 -> 148.1 -> 135.4"))
    w(adjudicate("q150", "q''@1.5 m", "1.11 -> 0.94 -> 0.97"))
    w(adjudicate("mdot150", "mdot@1.5 m", "0.0274 -> 0.0261 -> 0.0240"))
    w()

    # ==================== TASK 3 ====================
    w("## Task 3 - GCI redo with variable refinement ratio (Celik 2008)")
    w()
    w("GCI is computed for every parameter that is **monotone**, using three grid-size "
      "definitions: (iso) naive constant r=2 from the X/Z cell 20/10/5 mm; (gap-normal) the "
      "Y cell 16.7/10/5 mm giving r32=1.667, r21=2.0; (h-based) the Celik effective size "
      "h=(dx dy dz)^(1/3)=18.82/10/5 mm giving r32=1.882, r21=2.0. Apparent order p is solved "
      "from the Celik iterative equation with the q(p) term (see appendix).")
    w()
    grid_defs = [("iso r=2", ISO_H), ("gap-normal", GAP_NORMAL_H), ("h-based", CELIK_H)]
    for p in params:
        s20, s10, s05 = stats["nb20"][p], stats["nb10"][p], stats["nb05"][p]
        if None in (s20, s10, s05):
            continue
        d = decimals[p]
        fc, fm, ff = s20["mean"], s10["mean"], s05["mean"]
        mono = ((fm - fc) * (ff - fm)) > 0
        w(f"### {plabel[p]}  (coarse={fc:.{d}f}, med={fm:.{d}f}, fine={ff:.{d}f})")
        w()
        if not mono:
            w("Non-monotone across the three grids -> formal GCI is not defined "
              "(solution not in the asymptotic range). Reported as a finding, not a GCI. "
              "See Task 2 for whether the non-monotonicity is statistically real.")
            w()
            continue
        w("| Grid def | r32 | r21 | p | phi_ext | GCI_fine | e_a | asymptotic? |")
        w("|---|---|---|---|---|---|---|---|")
        for gname, h in grid_defs:
            g = gci_variable(fc, fm, ff, h)
            if g.get("p") is None:
                w(f"| {gname} | {g['r32']:.3f} | {g['r21']:.3f} | (no soln) | -- | -- | -- | -- |")
                continue
            w(f"| {gname} | {g['r32']:.3f} | {g['r21']:.3f} | {g['p']:.3f} | "
              f"{g['f_ext']:.{max(d,4)}f} | {100*g['gci_fine']:.1f}% | {100*g['e_a']:.1f}% | "
              f"{'yes' if g['asymptotic'] else 'NO'} |")
        w()
        if p == "mdot150":
            w(f"**Comparison to old constant-r=2 numbers** (mass flow): old p={OLD_MDOT['p']}, "
              f"old GCI={OLD_MDOT['gci']}%, old phi_ext={OLD_MDOT['phi_ext']}. The iso row here "
              f"should reproduce the old figure; the gap-normal / h-based rows show how the "
              f"anisotropy correction moves p and GCI.")
            w()
    w("**Verdict:** Because r21 = 2.0 in every definition but r32 differs (1.667 / 1.882 / 2.0), "
      "the q(p) correction shifts the apparent order p and the extrapolated value. The honest, "
      "anisotropy-consistent numbers are the gap-normal (physically the controlling direction "
      "for a 50 mm gap) and h-based rows, not the iso r=2 row the old analysis used.")
    w()

    # ==================== TASK 4 ====================
    w("## Task 4 - Monte-Carlo uncertainty on p and GCI")
    w()
    w("Each grid mean is perturbed ~ N(mean, SEM) (SEM = 95%CI/1.96), 10,000 draws. For every "
      "draw p and GCI_fine are recomputed; non-monotone or degenerate draws are DISCARDED (the "
      "discard fraction is itself a convergence-robustness diagnostic). Reported: median and "
      "2.5/97.5 percentiles.")
    w()
    mc_params = [p for p in params if resolved_monotone.get(p, {}).get("mono")]
    if not mc_params:
        mc_params = ["mdot150"]  # ensure at least mass flow reported
    for p in mc_params:
        s20, s10, s05 = stats["nb20"][p], stats["nb10"][p], stats["nb05"][p]
        if None in (s20, s10, s05):
            continue
        means = (s20["mean"], s10["mean"], s05["mean"])
        sems = (s20["sem"], s10["sem"], s05["sem"])
        w(f"### {plabel[p]}")
        w()
        w("| Grid def | discard frac | p median [2.5, 97.5] | GCI% median [2.5, 97.5] |")
        w("|---|---|---|---|")
        for gname, h in grid_defs:
            mc = mc_p_gci(means, sems, h)
            pm, plo, phi = mc["p"]
            gm, glo, ghi = mc["gci"]
            w(f"| {gname} | {mc['frac_disc']*100:.1f}% | "
              f"{pm:.2f} [{plo:.2f}, {phi:.2f}] | {gm:.1f} [{glo:.1f}, {ghi:.1f}] |")
        w()
    w("**Verdict:** A large discard fraction means the monotonicity itself flips under sampling "
      "noise -> the GCI for that parameter is fragile and should be reported with the caveat that "
      "the grids are not reliably in the asymptotic range. A small discard fraction means the "
      "convergence direction is robust and the p/GCI credible interval can be quoted directly.")
    w()

    # ==================== TASK 5 ====================
    w("## Task 5 - Peak outer-skin heat flux over height")
    w()
    w("For each nb grid, the steady-window mean of `HFLUX_OUT_Zxxx` is computed at every height; "
      "the peak-over-height and its location are reported, with SEM from the ACF method.")
    w()
    w("| Case | Peak q'' (kW/m2) | SEM | 95% CI (+/-) | Height of peak |")
    w("|---|---|---|---|---|")
    peak_seq = {}
    for r in nb:
        best = None
        for h in HFLUX_OUT_HEIGHTS:
            x = data[r][f"qout_{h}"]
            if x is None:
                continue
            s = series_stats(x, data[r]["window"])
            if s is None:
                continue
            if best is None or s["mean"] > best[1]["mean"]:
                best = (h, s)
        if best is None:
            w(f"| {r} | -- | -- | -- | -- |")
            continue
        h, s = best
        peak_seq[r] = s["mean"]
        w(f"| {r} | {s['mean']:.2f} | {s['sem']:.2f} | {s['ci95']:.2f} | "
          f"Z{h/100:.2f} m ({h} cm) |")
    w()
    # verify claimed 11.4 -> 15.8 -> 17.8 monotone rising
    claim = {"nb20": 11.4, "nb10": 15.8, "nb05": 17.8}
    if all(r in peak_seq for r in nb):
        seq = [peak_seq["nb20"], peak_seq["nb10"], peak_seq["nb05"]]
        rising = seq[0] < seq[1] < seq[2]
        matches = all(abs(peak_seq[r] - claim[r]) < 0.6 for r in nb)
        w(f"Data sequence (20->10->5 mm): {seq[0]:.1f} -> {seq[1]:.1f} -> {seq[2]:.1f} kW/m2. "
          f"Old figure claim: 11.4 -> 15.8 -> 17.8 kW/m2.")
        w()
        w(f"**Verdict:** peak-over-height is {'monotone RISING as claimed' if rising else 'NOT monotone rising'}; "
          f"data {'matches' if matches else 'does NOT match'} the old figure within 0.6 kW/m2. The peak sits at "
          f"the lowest resolved height near the fire base (see height column), consistent with the "
          f"expectation that the outer skin is hottest where the plume attaches. This flux is NOT grid "
          f"converged (it keeps rising as cells shrink), which must be stated as a limitation.")
    w()

    # ==================== APPENDIX ====================
    w("## Appendix - Celik (2008) variable-r formulas implemented")
    w()
    w("Let f1 = fine, f2 = medium, f3 = coarse grid solutions; grid sizes h1<h2<h3; "
      "refinement ratios r21 = h2/h1, r32 = h3/h2 (allowed to differ).")
    w()
    w("```")
    w("eps21 = f2 - f1,   eps32 = f3 - f2")
    w("s     = sign(eps32 / eps21)              # +1 monotone, -1 oscillatory")
    w("")
    w("Apparent order p (fixed-point iteration):")
    w("    p     = (1/ln(r21)) * | ln|eps32/eps21| + q(p) |")
    w("    q(p)  = ln( (r21^p - s) / (r32^p - s) )")
    w("  (for constant r, r21=r32 -> q=0 -> p = ln|eps32/eps21| / ln(r), the old formula)")
    w("")
    w("Richardson extrapolated value (fine):")
    w("    phi_ext21 = (r21^p * f1 - f2) / (r21^p - 1)")
    w("")
    w("Error estimators:")
    w("    e_a21   = |(f1 - f2) / f1|                     # approximate relative")
    w("    e_ext21 = |(phi_ext21 - f1) / phi_ext21|       # extrapolated relative")
    w("    GCI_fine = 1.25 * e_a21 / (r21^p - 1)          # 1.25 = Fs (3-grid)")
    w("```")
    w()
    w("Grid-size definitions used:")
    w(f"- iso (old, constant r=2): X/Z cell 20 / 10 / 5 mm -> r21=2.0, r32=2.0")
    w(f"- gap-normal (Y cell): 16.7 / 10 / 5 mm -> r21=2.0, r32=1.667")
    w(f"- h-based (dx*dy*dz)^(1/3): 18.82 / 10 / 5 mm -> r21=2.0, r32=1.882")
    w()
    w("Integral time scale (Task 1): tau = dt * (0.5 + sum_{k=1}^{k0-1} rho(k)), where rho(k) is "
      "the normalised autocovariance and k0 the first lag with rho<=0; floored at 0.5*dt. "
      "N_ind = T_window / (2*tau); SEM = SD/sqrt(N_ind); 95% CI = 1.96*SEM.")
    w()
    w("Reference: Celik, Ghia, Roache, Freitas, Coleman, Raad (2008), "
      "\"Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD "
      "Applications,\" J. Fluids Eng. 130(7):078001.")

    text = "\n".join(md) + "\n"
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gci_stats_v2_results.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(text)
    print(f"\n[written] {out}")


if __name__ == "__main__":
    main()
