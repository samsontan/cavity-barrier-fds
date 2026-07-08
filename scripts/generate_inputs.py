"""
Paper B (FSJ) -- Cavity barrier CFD campaign: FDS input generator.

REV 2 (2026-06-12): rig descaled to Livkiss laboratory class after measured
throughput on the X1 Carbon (Core Ultra 7 268V) showed the original 3.0 m
one-storey rig would take 2-3 weeks. Same 50 mm gap, same grid series, same
barrier states, same 31 kW/m line intensity -- the science is unchanged and
the validation-class match to Livkiss et al. (Fire Technology 2018 / FSJ 2019)
is now tighter.

Generates 9 FDS input files:

  Grid series (cavity gap 50 mm):
    dx = 20 mm  -> gap/dx = 2.5   (satisfies D*/dx ~ 6 -- the published-study trap)
    dx = 10 mm  -> gap/dx = 5
    dx =  5 mm  -> gap/dx = 10    (Livkiss et al. 2019 threshold)

  Barrier states at each grid:
    nb = no barrier
    cl = closed state (full 50 mm blockage, stone wool, 80 mm tall)
    op = open state   (20 mm intumescent strip on inner leaf, 30 mm air gap)

Rig: representative non-combustible EWS cavity, laboratory scale.
  X = 0.4 m wide, Y = 50 mm gap, Z = 0..1.8 m boards, domain Z = -0.3..2.1 m.
  Fire: propane line burner 0.2 m x 50 mm at base, 6.2 kW prescribed
  (31 kW/m line intensity -- Livkiss-class). Non-combustible boards by design:
  avoids the ABCB report SS6.4.2 material-property critique and isolates the
  barrier aerodynamic/thermal effect. Barrier at Z = 1.20-1.28 m (spandrel zone).

Run folders: C:\\FDS_runs\\cavb\\<chid>\\
"""

import os

RUN_ROOT = r"C:\FDS_runs\cavb"
INPUT_ARCHIVE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "input")

# Domain
X0, X1 = 0.0, 0.4
Y0, Y1 = 0.0, 0.05          # 50 mm cavity gap
ZB0, ZB1 = 0.0, 1.8         # board span
ZD0, ZD1 = -0.3, 2.1        # domain span

# Fire (31 kW/m line intensity, Livkiss-class)
BURNER_X0, BURNER_X1 = 0.1, 0.3
BURNER_ZTOP = 0.06
HRR_KW = 6.2
# FDS HRRPUA is in kW/m2 (NOT W/m2) -- 6.2 kW / 0.01 m2 = 620 kW/m2
HRRPUA = HRR_KW / ((BURNER_X1 - BURNER_X0) * (Y1 - Y0))   # kW/m2
D_STAR_MM = 125.0           # D* = 0.125 m for 6.2 kW

# Barrier (Z = 1.20-1.28 m -- clear of every mesh split plane by >= 2 cells)
BAR_Z0, BAR_Z1 = 1.20, 1.28
OPEN_STRIP_Y1 = 0.02        # open state: strip occupies Y = 0..20 mm on inner leaf

T_END = 120.0               # statistics window 60-120 s
T_END_BY_GRID = {"20": 120.0, "10": 120.0, "05": 100.0}  # 05 trimmed (stats 60-100 s)

GRIDS = {
    # label: (meshes [(I,J,K, z0,z1)], n_mpi, omp_threads)
    "20": ([(20, 3, 120, ZD0, ZD1)], 1, 4),
    "10": ([(40, 5, 60, -0.3, 0.3), (40, 5, 60, 0.3, 0.9),
            (40, 5, 60, 0.9, 1.5), (40, 5, 60, 1.5, 2.1)], 4, 2),
    "05": ([(80, 10, 80, -0.3, 0.1), (80, 10, 80, 0.1, 0.5),
            (80, 10, 80, 0.5, 0.9), (80, 10, 80, 0.9, 1.3),
            (80, 10, 80, 1.3, 1.7), (80, 10, 80, 1.7, 2.1)], 6, 1),
}

BARRIERS = ["nb", "cl", "op"]

# Thermocouple tree heights -- avoids the 1.20-1.28 barrier band
TC_Z = [0.2, 0.4, 0.6, 0.8, 1.0, 1.15, 1.35, 1.5, 1.65, 1.75]
VEL_Z = [0.4, 0.8, 1.15, 1.35, 1.5, 1.7]
HF_Z = [0.4, 0.8, 1.0, 1.35, 1.5, 1.7]   # heat flux / wall temp heights
XC = 0.20                                 # rig centreline


def build_fds(chid, grid_label, barrier, t_end=None):
    if t_end is None:
        t_end = T_END_BY_GRID.get(grid_label, T_END)
    meshes, n_mpi, _omp = GRIDS[grid_label]
    dx = {"20": 20, "10": 10, "05": 5}[grid_label]
    gap_dx = 50.0 / dx
    bar_desc = {"nb": "no cavity barrier",
                "cl": "CLOSED-state barrier (full 50 mm blockage)",
                "op": "OPEN-state barrier (20 mm strip, 30 mm residual gap)"}[barrier]

    L = []
    a = L.append
    a("! =============================================================================")
    a("! PAPER B (FSJ) -- CAVITY BARRIER EFFICACY IN A VENTILATED EWS CAVITY (REV 2)")
    a("! =============================================================================")
    a(f"! Scenario : {bar_desc}")
    a(f"! Grid     : dx = {dx} mm nominal -> cavity-gap/dx = {gap_dx:.1f}")
    a("! Rig      : 0.4 m wide x 50 mm gap x 1.8 m tall non-combustible EWS cavity")
    a("!            inner leaf 12 mm fibre cement, outer skin 3 mm aluminium")
    a(f"! Fire     : propane line burner 0.2 m x 50 mm, {HRR_KW} kW prescribed (31 kW/m)")
    a("! Protocol : ABCB report SS6.4 -- grid convergence on secondary parameters")
    a("!            (cavity gas temperature, wall heat flux, mass flow) because the")
    a("!            fire is prescribed; cavity gap (not D*) is the controlling scale.")
    a(f"! D* check : D* = 0.125 m for 6.2 kW -> D*/dx = {D_STAR_MM/dx:.1f} at this grid")
    a("! =============================================================================")
    a("")
    a(f"&HEAD CHID='{chid}', TITLE='Cavity barrier EWS -- {bar_desc}, dx={dx}mm' /")
    a("")
    a(f"&TIME T_END={t_end:.1f} /")
    a("")
    a("&DUMP DT_SLCF=2.0, DT_DEVC=0.5, DT_BNDF=15.0, DT_HRR=1.0 /")
    a("")
    a("! ===== MESH =====")
    for (I, J, K, z0, z1) in meshes:
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
    a("! Cavity leaves: Y boundary patches over the board span; OPEN above and below")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {ZB0:.2f},{ZB1:.2f}, SURF_ID='INNER_LEAF' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {ZD0:.2f},{ZB0:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y0:.2f}, {ZB1:.2f},{ZD1:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {ZB0:.2f},{ZB1:.2f}, SURF_ID='OUTER_SKIN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {ZD0:.2f},{ZB0:.2f}, SURF_ID='OPEN' /")
    a(f"&VENT XB={X0:.2f},{X1:.2f}, {Y1:.2f},{Y1:.2f}, {ZB1:.2f},{ZD1:.2f}, SURF_ID='OPEN' /")
    a("&VENT MB='ZMIN', SURF_ID='OPEN' /")
    a("&VENT MB='ZMAX', SURF_ID='OPEN' /")
    a("! XMIN / XMAX default INERT -- sealed cavity edges (test-rig convention)")
    a("")
    a("! ===== FIRE SOURCE =====")
    a(f"&OBST XB={BURNER_X0:.2f},{BURNER_X1:.2f}, {Y0:.2f},{Y1:.2f}, 0.00,{BURNER_ZTOP:.2f}, SURF_ID='INERT' /")
    a(f"&VENT XB={BURNER_X0:.2f},{BURNER_X1:.2f}, {Y0:.2f},{Y1:.2f}, {BURNER_ZTOP:.2f},{BURNER_ZTOP:.2f}, SURF_ID='BURNER' /")
    a("")
    if barrier == "cl":
        a("! ===== CAVITY BARRIER -- CLOSED STATE (full blockage) =====")
        a(f"&OBST XB={X0:.2f},{X1:.2f}, {Y0:.2f},{Y1:.2f}, {BAR_Z0:.2f},{BAR_Z1:.2f}, SURF_ID='BARRIER_CL', BNDF_OBST=.TRUE. /")
    elif barrier == "op":
        a("! ===== CAVITY BARRIER -- OPEN STATE (strip on inner leaf, 30 mm gap) =====")
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
        X0, X1, Y0, Y1, ZD0, ZD1))
    a("")
    a("! ===== BOUNDARY FILES =====")
    a("&BNDF QUANTITY='WALL TEMPERATURE' /")
    a("&BNDF QUANTITY='NET HEAT FLUX' /")
    a("")
    a("&TAIL /")
    a("")
    return "\n".join(L), n_mpi


def main():
    nmpi, nomp = {}, {}
    for grid in ["20", "10", "05"]:
        for barrier in BARRIERS:
            chid = f"cavb_{barrier}{grid}"
            text, n_mpi = build_fds(chid, grid, barrier)
            run_dir = os.path.join(RUN_ROOT, chid)
            os.makedirs(run_dir, exist_ok=True)
            for dest in (os.path.join(run_dir, f"{chid}.fds"),
                         os.path.join(INPUT_ARCHIVE, f"{chid}.fds")):
                with open(dest, "w") as f:
                    f.write(text)
            nmpi[chid] = n_mpi
            nomp[chid] = GRIDS[grid][2]
            print(f"  wrote {chid}.fds  (n_mpi={n_mpi}, omp={nomp[chid]})")
    # Queue order: cheap -> expensive; fields: chid,n_mpi,omp
    order = [f"cavb_{b}{g}" for g in ["20", "10", "05"] for b in BARRIERS]
    with open(os.path.join(RUN_ROOT, "queue.txt"), "w") as f:
        for chid in order:
            f.write(f"{chid},{nmpi[chid]},{nomp[chid]}\n")
    print(f"Queue written: {os.path.join(RUN_ROOT, 'queue.txt')}")


if __name__ == "__main__":
    main()
