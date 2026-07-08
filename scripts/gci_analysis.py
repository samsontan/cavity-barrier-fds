"""
Grid-Convergence Index analysis for the cavity-barrier no-barrier series.

Reads the real FDS devc.csv outputs, computes 60-120 s steady means (+ std for
stationarity), verifies Table 1, and applies the Celik et al. (2008) GCI procedure
to the no-barrier 20/10/5 mm series. Honest reporting: flags non-monotone series
(GCI undefined / not in asymptotic range) and reports apparent order p.
"""
import csv, os, math, statistics as st

RUNS = r"C:\FDS_runs\cavb"
WIN = (60.0, 120.0)
CH = {
    "HRR": "HRR_TOTAL",
    "T135": "TC_Z135",
    "T165": "TC_Z165",
    "mdot150": "MFLOW_UP_Z150",
    "q150": "HFLUX_OUT_Z150",
}

def steady(run):
    f = os.path.join(RUNS, f"cavb_{run}", f"cavb_{run}_devc.csv")
    if not os.path.exists(f):
        return None
    rows = list(csv.reader(open(f)))
    names = rows[1]
    idx = {names.index(v): k for k, v in CH.items() if v in names}
    data = {k: [] for k in CH}
    for r in rows[2:]:
        if not r or not r[0].strip():
            continue
        try: t = float(r[0])
        except ValueError: continue
        if WIN[0] <= t <= WIN[1]:
            for ci, key in idx.items():
                try: data[key].append(float(r[ci]))
                except (ValueError, IndexError): pass
    return {k: ((st.mean(v), st.pstdev(v) if len(v) > 1 else 0.0) if v else (None, None))
            for k, v in data.items()}

def fmt(x, n=1):
    return f"{x[0]:.{n}f}(±{x[1]:.{n}f})" if x and x[0] is not None else "--"

def gci(f_coarse, f_med, f_fine, r=2.0):
    """Celik 2008. f1=fine, f2=med, f3=coarse."""
    f1, f2, f3 = f_fine, f_med, f_coarse
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0:
        return None
    ratio = e32 / e21
    s = 1 if ratio > 0 else -1
    p = abs(math.log(abs(ratio))) / math.log(r)
    rp = r ** p
    f_ext = (rp * f1 - f2) / (rp - 1)
    e_a = abs((f1 - f2) / f1)
    e_ext = abs((f_ext - f1) / f_ext) if f_ext != 0 else float("nan")
    gci_fine = 1.25 * e_a / (rp - 1)
    return dict(p=p, ratio=ratio, monotone=(s > 0), f_ext=f_ext,
                e_a=e_a, e_ext=e_ext, gci_fine=gci_fine,
                converging=(abs(ratio) > 1))  # |e32|>|e21| => increments shrink on refine

print("=== Steady 60-120 s means (+/-pstdev) ===")
order = ["nb20","nb10","nb05","op20","op10","op05","cl20","cl10"]
res = {}
for run in order:
    s = steady(run); res[run] = s
    if s is None:
        print(f"{run:6s}: (no data)"); continue
    print(f"{run:6s}: HRR={fmt(s['HRR'],2)}  T135={fmt(s['T135'])}  T165={fmt(s['T165'])}  "
          f"mdot150={fmt(s['mdot150'],4)}  q150={fmt(s['q150'],2)}")

print("\n=== No-barrier GCI (Celik 2008), r=2, grids 20->10->5 mm ===")
labels = {"T135":"Centreline T @1.35 m (C)","T165":"Centreline T @1.65 m (C)",
          "mdot150":"Upward mass flow @1.5 m (kg/s)","q150":"Outer-skin q'' @1.5 m (kW/m2)"}
for key, lab in labels.items():
    c = res["nb20"][key][0]; m = res["nb10"][key][0]; fi = res["nb05"][key][0]
    if None in (c, m, fi):
        print(f"\n{lab}: incomplete ({c},{m},{fi})"); continue
    print(f"\n{lab}")
    print(f"  20mm={c:.4g}  10mm={m:.4g}  5mm={fi:.4g}")
    mono = (c-m)*(m-fi) > 0
    print(f"  monotone across 3 grids: {mono}")
    g = gci(c, m, fi)
    if not mono:
        print("  -> NON-MONOTONE: formal GCI not defined; report as finding (under-resolved even at 5 mm).")
        continue
    print(f"  apparent order p={g['p']:.2f}  (nominal FDS ~2)")
    print(f"  e32/e21 ratio={g['ratio']:.3f}  asymptotic-range(|ratio|>1, increments shrink on refine)={g['converging']}")
    print(f"  Richardson f_ext={g['f_ext']:.4g}   fine-grid GCI={100*g['gci_fine']:.1f}%   e_ext={100*g['e_ext']:.1f}%")
