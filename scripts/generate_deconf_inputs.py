"""
Paper B (FSJ) -- de-confounding campaign (Major Revision, reviewer request).

Four tractable verification runs on top of the REV 2 rig (see generate_inputs.py
for the baseline geometry; all physics, materials, fire, devices and slices are
IDENTICAL to the baseline campaign so results are directly comparable):

  1. cavb_nb10dm  domain-extent      : nb @ dx=10 mm, Z margins doubled
                                       (-0.6..2.7 m vs -0.3..2.1 m).
                                       ACCEPT: mdot@1.5, TC tree within GCI noise
                                       of cavb_nb10 -> draught not boundary-limited.
  2. cavb_nb05if  mesh-interface     : nb @ dx=5 mm, re-decomposed so NO MPI
                                       interface sits at Z=1.3 m (old split) nor
                                       within 50 mm of any TC height.
                                       ACCEPT: metrics within GCI noise of
                                       cavb_nb05 -> interfaces do not contaminate
                                       the convergence series.
  3. cavb_cl10lk  closed leakage     : cl @ dx=10 mm, barrier shortened 10 mm at
                                       each X end (2 edge slots, 5% of cavity
                                       cross-section) -- realistic imperfect seal.
                                       METRIC: how fast the ~2/3 HRR starvation
                                       collapses with leakage.
  4. cavb_op20ds  de-snapped open    : op @ dx=20 mm in X,Z but dy=10 mm (J=5)
                                       so the 20 mm strip resolves EXACTLY
                                       (baseline op20 snaps it to 16.7 mm).
                                       METRIC: op20ds vs op20 vs op10 -- separates
                                       geometry snapping from grid resolution.

Run folders: C:\\FDS_runs\\cavb\\<chid>\\ ; queue file queue_deconf.txt.
Launch with run_queue_deconf.py (detached, pythonw).
"""

import os

RUN_ROOT = r"C:\FDS_runs\cavb"
INPUT_ARCHIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "input")

# Baseline rig constants (MUST match generate_inputs.py REV 2)
X0, X1 = 0.0, 0.4
Y0, Y1 = 0.0, 0.05
ZB0, ZB1 = 0.0, 1.8
BURNER_X0, BURNER_X1 = 0.1, 0.3
BURNER_ZTOP = 0.06
HRR_KW = 6.2
HRRPUA = HRR_KW / ((BURNER_X1 - BURNER_X0) * (Y1 - Y0))   # 620 kW/m2 (kW-based!)
D_STAR_MM = 125.0
BAR_Z0, BAR_Z1 = 1.20, 1.28
OPEN_STRIP_Y1 = 0.02
TC_Z = [0.2, 0.4, 0.6, 0.8, 1.0, 1.15, 1.35, 1.5, 1.65, 1.75]
VEL_Z = [0.4, 0.8, 1.15, 1.35, 1.5, 1.7]
HF_Z = [0.4, 0.8, 1.0, 1.35, 1.5, 1.7]
XC = 0.20

# ---------------------------------------------------------------------------
# Variant definitions
# ---------------------------------------------------------------------------
VARIANTS = {
    "cavb_nb10dm": dict(
        desc="DOMAIN-EXTENT verification: no barrier, dx=10 mm, Z margins doubled",
        dx=10, barrier="nb", t_end=120.0, n_mpi=4, omp=2,
        zd0=-0.6, zd1=2.7,
        meshes=[(40, 5, 90, -0.6, 0.3), (40, 5, 90, 0.3, 1.2),
                (40, 5, 90, 1.2, 2.1), (40, 5, 60, 2.1, 2.7)],
    ),
    "cavb_nb05if": dict(
        desc="MESH-INTERFACE verification: no barrier, dx=5 mm, splits moved off Z=1.3 m",
        dx=5, barrier="nb", t_end=100.0, n_mpi=6, omp=1,
        zd0=-0.3, zd1=2.1,
        # old splits: 0.1/0.5/0.9/1.3/1.7 -- Z=1.3 sat in the plume path.
        # new splits: 0.1/0.5/0.9/1.40/1.75 -- clear of Z=1.3 and >=50 mm from
        # every TC height (nearest: TC 1.35 vs split 1.40, 10 cells away).
        meshes=[(80, 10, 80, -0.3, 0.1), (80, 10, 80, 0.1, 0.5),
                (80, 10, 80, 0.5, 0.9), (80, 10, 100, 0.9, 1.4),
                (80, 10, 70, 1.4, 1.75), (80, 10, 70, 1.75, 2.1)],
    ),
    "cavb_cl10lk": dict(
        desc="CLOSED-STATE LEAKAGE: dx=10 mm, 10 mm edge slot each end (5% leakage area)",
        dx=10, barrier="cl_leak", t_end=120.0, n_mpi=4, omp=2,
        zd0=-0.3, zd1=2.1,
        meshes=[(40, 5, 60, -0.3, 0.3), (40, 5, 60, 0.3, 0.9),
                (40, 5, 60, 0.9, 1.5), (40, 5, 60, 1.5, 2.1)],
        barrier_x0=0.01, barrier_x1=0.39,   # exactly 1 cell open at each end
    ),
    "cavb_op20ds": dict(
        desc="DE-SNAPPED OPEN STATE: dx=20 mm X/Z, dy=10 mm so 20 mm strip is exact",
        dx=20, barrier="op", t_end=120.0, n_mpi=1, omp=4,
        zd0=-0.3, zd1=2.1,
        meshes=[(20, 5, 120, -0.3, 2.1)],   # J=5 -> dy=10 mm (baseline op20: J=3, dy=16.7)
    ),
}


def build_variant(chid, v):
    dx = v["dx"]
    gap_dx = 50.0 / dx
    zd0, zd1 = v["zd0"], v["zd1"]
    t_end = v["t_end"]

    L = []
    a = L.append
    a("! =============================================================================")
    a("! PAPER B (FSJ) -- DE-CONFOUNDING CAMPAIGN (Major Revision)")
    a(f"! {v['desc']}")
    a("! =============================================================================")
    a(f"! Baseline  : identical rig/physics to generate_inputs.py REV 2 campaign")
    a(f"! Grid      : dx = {dx} mm nominal -> cavity-gap/dx = {gap_dx:.1f}")
    a(f"! Fire      : propane line burner 0.2 m x 50 mm, {HRR_KW} kW prescribed (31 kW/m)")
    a(f"! D* check  : D* = 0.125 m for 6.2 kW -> D*/dx = {D_STAR_MM/dx:.1f} at this grid")
    a("! =============================================================================")
    a("")
    a(f"&HEAD CHID='{chid}', TITLE='Cavity barrier de-confounding -- {v['desc'][:60]}' /")
    a("")
    a(f"&TIME T_END={t_end:.1f} /")
    a("")
    a("&DUMP DT_SLCF=2.0, DT_DEVC=0.5, DT_BNDF=15.0, DT_HRR=1.0 /")
    a("")
    a("! ===== MESH =====")
    for (I, J, K, z0, z1) in v["meshes"]:
        a(f"&MESH IJK={I},{J},{K}, XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y1:.2f}, {z0:.2f},{z1:.2f} /")
    a("")
    a("! ===== COMBUSTION =====")
    a("&REAC FUEL='PROPANE', SOOT_YIELD=0.024, CO_YIELD=0.005 /")
    a("")
    a("! ===== MATERIALS =====")
    a("&MATL ID='FIBRE_CEMENT', CONDUCTIVITY=0.25, SPECIFIC_HEAT=1.0, DENSITY=1350. /")
    a("&MATL ID='ALUMINIUM',    CONDUCTIVITY=237., SPECIFIC_HEAT=0.897, DENSITY=2700. /")
    a("&MATL ID='STONE_WOOL',   CONDUCTIVITY=0.04, SPECIFIC_HEAT=0.84, DENSITY=100. /")
    a("&MATL ID='INTUMESCENT',  CONDUCTIVITY=0.20, SPECIFIC_HEAT=1.5,  DENSITY=1000. /")
    a("")
    a("! ===== SURFACES =====")
    a("&SURF ID='INNER_LEAF', MATL_ID='FIBRE_CEMENT', THICKNESS=0.012, BACKING='VOID', COLOR='BEIGE' /")
    a("&SURF ID='OUTER_SKIN', MATL_ID='ALUMINIUM',    THICKNESS=0.003, BACKING='VOID', COLOR='GRAY 60' /")
    a("&SURF ID='BARRIER_CL', MATL_ID='STONE_WOOL',   THICKNESS=0.025, BACKING='INSULATED', COLOR='MAGENTA' /")
    a("&SURF ID='BARRIER_OP', MATL_ID='INTUMESCENT',  THICKNESS=0.010, BACKING='INSULATED', COLOR='MAGENTA' /")
    a("&SURF ID='BURNER',")
    a(f"      HRRPUA={HRRPUA:.1f},")
    a("      RAMP_Q='FIRE_RAMP',")
    a("      COLOR='RED' /")
    a("")
    a("&RAMP ID='FIRE_RAMP', T=0.0,  F=0.0 /")
    a("&RAMP ID='FIRE_RAMP', T=10.0, F=1.0 /")
    a(f"&RAMP ID='FIRE_RAMP', T={t_end:.1f}, F=1.0 /")
    a("")
    a("! ===== BOUNDARIES =====")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {ZB0:.2f},{ZB1:.2f}, SURF_ID='INNER_LEAF' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {zd0:.2f},{ZB0:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {ZB1:.2f},{zd1:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {ZB0:.2f},{ZB1:.2f}, SURF_ID='OUTER_SKIN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {zd0:.2f},{ZB0:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {ZB1:.2f},{zd1:.2f}, SURF_ID='OPEN' /")
    a("&VENT MB='ZMIN', SURF_ID='OPEN' /")
    a("&VENT MB='ZMAX', SURF_ID='OPEN' /")
    a("! XMIN / XMAX default INERT -- sealed cavity edges (test-rig convention)")
    a("")
    a("! ===== FIRE SOURCE =====")
    a(f"&OBST XB={BURNER_X0:.2f},{BURNER_X1:.2f}, {Y0:.2f},{Y1:.2f}, 0.00,{BURNER_ZTOP:.2f}, SURF_ID='INERT' /")
    a(f"&VENT XB={BURNER_X0:.2f},{BURNER_X1:.2f}, {Y0:.2f},{Y1:.2f}, {BURNER_ZTOP:.2f},{BURNER_ZTOP:.2f}, SURF_ID='BURNER' /")
    a("")
    barrier = v["barrier"]
    if barrier == "cl_leak":
        bx0, bx1 = v["barrier_x0"], v["barrier_x1"]
        a("! ===== CAVITY BARRIER -- CLOSED STATE WITH EDGE LEAKAGE =====")
        a(f"! Slots X={X0:.2f}-{bx0:.2f} and {bx1:.2f}-{X1:.2f} (1 cell each end at dx=10 mm)")
        a(f"! Leakage area = 2 x {bx0-X0:.2f} x {Y1-Y0:.2f} = {2*(bx0-X0)*(Y1-Y0)*1e4:.1f} cm2 = "
          f"{100*2*(bx0-X0)/(X1-X0):.0f}% of cavity cross-section")
        a(f"&OBST XB={bx0:.2f},{bx1:.2f}, {Y0:.2f},{Y1:.2f}, {BAR_Z0:.2f},{BAR_Z1:.2f}, SURF_ID='BARRIER_CL', BNDF_OBST=.TRUE. /")
    elif barrier == "op":
        a("! ===== CAVITY BARRIER -- OPEN STATE (strip on inner leaf, 30 mm gap) =====")
        a(f"! dy = 10 mm here -> strip Y=0.00-{OPEN_STRIP_Y1:.2f} resolves EXACTLY (2 cells);")
        a("! baseline op20 (dy=16.7 mm) snapped it to 16.7 mm -- that is the confound under test")
        a(f"&OBST XB={X0:.2f},{X1:.2f}, {Y0:.2f},{OPEN_STRIP_Y1:.2f}, {BAR_Z0:.2f},{BAR_Z1:.2f}, SURF_ID='BARRIER_OP', BNDF_OBST=.TRUE. /")
    else:
        a("! ===== NO CAVITY BARRIER =====")
    a("")
    a("! ===== SLICES =====")
    a(f"&SLCF PBX={XC:.2f}, QUANTITY='TEMPERATURE' /")
    a(f"&SLCF PBX={XC:.2f}, QUANTITY='W-VELOCITY' /")
    a(f"&SLCF PBX={XC:.2f}, QUANTITY='HRRPUV' /")
    a("&SLCF PBY=0.025, QUANTITY='TEMPERATURE' /")
    a("&SLCF PBY=0.025, QUANTITY='W-VELOCITY' /")
    a("&SLCF PBY=0.025, QUANTITY='HRRPUV' /")
    a("&SLCF PBZ=1.50, QUANTITY='TEMPERATURE' /")
    a("")
    a(f"! ===== DEVICES -- centreline thermocouple tree (X={XC:.2f}, Y=0.025) =====")
    for z in TC_Z:
        a(f"&DEVC ID='TC_Z{int(round(z*100)):03d}', XYZ={XC:.2f},0.025,{z:.2f}, QUANTITY='TEMPERATURE' /")
    a("")
    a("! ===== DEVICES -- vertical velocity =====")
    for z in VEL_Z:
        a(f"&DEVC ID='W_Z{int(round(z*100)):03d}', XYZ={XC:.2f},0.025,{z:.2f}, QUANTITY='W-VELOCITY' /")
    a("")
    a("! ===== DEVICES -- outer skin heat flux + wall temperature (IOR=-2) =====")
    for z in HF_Z:
        a(f"&DEVC ID='HFLUX_OUT_Z{int(round(z*100)):03d}', XYZ={XC:.2f},{Y1:.2f},{z:.2f}, IOR=-2, QUANTITY='NET HEAT FLUX' /")
        a(f"&DEVC ID='TWALL_OUT_Z{int(round(z*100)):03d}', XYZ={XC:.2f},{Y1:.2f},{z:.2f}, IOR=-2, QUANTITY='WALL TEMPERATURE' /")
    a("")
    a("! ===== DEVICES -- inner leaf heat flux above barrier (IOR=2) =====")
    for z in [1.35, 1.5, 1.7]:
        a(f"&DEVC ID='HFLUX_IN_Z{int(round(z*100)):03d}', XYZ={XC:.2f},{Y0:.2f},{z:.2f}, IOR=2, QUANTITY='NET HEAT FLUX' /")
    a("")
    a("! ===== DEVICES -- mass flow through horizontal planes =====")
    a(f"&DEVC ID='MFLOW_UP_Z080', XB={X0:.2f},{X1:.2f},{Y0:.2f},{Y1:.2f},0.80,0.80, QUANTITY='MASS FLOW +' /")
    a(f"&DEVC ID='MFLOW_UP_Z150', XB={X0:.2f},{X1:.2f},{Y0:.2f},{Y1:.2f},1.50,1.50, QUANTITY='MASS FLOW +' /")
    a("")
    a("&DEVC ID='HRR_TOTAL', XB={:.2f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}, QUANTITY='HRR' /".format(
        X0, X1, Y0, Y1, zd0, zd1))
    a("")
    a("! ===== BOUNDARY FILES =====")
    a("&BNDF QUANTITY='WALL TEMPERATURE' /")
    a("&BNDF QUANTITY='NET HEAT FLUX' /")
    a("")
    a("&TAIL /")
    a("")
    return "\n".join(L)


def main():
    queue = []
    # cheap -> expensive
    order = ["cavb_op20ds", "cavb_cl10lk", "cavb_nb10dm", "cavb_nb05if"]
    for chid in order:
        v = VARIANTS[chid]
        text = build_variant(chid, v)
        run_dir = os.path.join(RUN_ROOT, chid)
        os.makedirs(run_dir, exist_ok=True)
        for dest in (os.path.join(run_dir, f"{chid}.fds"),
                     os.path.join(INPUT_ARCHIVE, f"{chid}.fds")):
            with open(dest, "w") as f:
                f.write(text)
        queue.append(f"{chid},{v['n_mpi']},{v['omp']}")
        print(f"  wrote {chid}.fds  (n_mpi={v['n_mpi']}, omp={v['omp']})")
    qpath = os.path.join(RUN_ROOT, "queue_deconf.txt")
    with open(qpath, "w") as f:
        f.write("\n".join(queue) + "\n")
    print(f"Queue written: {qpath}")


if __name__ == "__main__":
    main()
