# cavity-barrier-fds

Simulation inputs, device outputs, and analysis scripts for:

> **A Grid-Resolution Discipline for Cavity-Barrier Fire Simulation: Why the Plume Criterion Fails in Thin Ventilated Façade Cavities**
> Samson B. H. Tan, Teoh Teik Toe, H. M. Iqbal Mahmud, Khalid A. M. Moinuddin
> Submitted to *Fire Safety Journal* (manuscript under preparation/review; citation to be updated on acceptance).

A non-combustible laboratory-scale ventilated cavity (50 mm gap, 6.2 kW propane
line fire) is modelled in FDS 6.10.1 at 3, 5, and 10 gap-normal cells across
three idealised barrier states (none / closed / open), plus a de-confounding
verification campaign. The paper argues the plume-based D*/δx mesh criterion is
unreliable in thin cavities and that grid convergence must be demonstrated on
gap-scale secondary parameters.

## Contents

| Path | What |
|---|---|
| `inputs/` | Complete FDS input files: 9 baseline campaign runs (`cavb_{nb,cl,op}{20,10,05}.fds`) + 5 verification runs (`cavb_op20ds`, `cavb_cl10lk`, `cavb_nb10dm`, `cavb_nb05if`, `cavb_op05if`) |
| `outputs/<chid>/` | Device (`*_devc.csv`) and heat-release (`*_hrr.csv`) time series. Present for 13 of the 14 inputs — see *Which runs have outputs* below |
| `scripts/generate_inputs.py`, `generate_deconf_inputs.py` | Programmatic generators for the input files |
| `scripts/gci_stats_v2.py` | Authoritative statistics: ACF-corrected sampling uncertainty, variable-r Celik GCI, Monte-Carlo p/GCI uncertainty (Tables 1–2) |
| `scripts/gci_analysis.py` | Superseded first-pass GCI script, kept for provenance |
| `scripts/fds_postprocess_cavb.py` | Direct Fortran-binary slice reader + campaign post-processing |
| `scripts/publication_figures.py`, `deconf_jet_figures.py`, `generate_setup_figure.py` | Figure generation (Figures 1–8) |
| `scripts/postprocess_op05if.py` | Outcome classification + steady statistics for the §4.5 disambiguation run |

## Status

All runs are final. `cavb_op05if` (fine-grid open state with relocated MPI
interfaces) completed to `T_END = 100 s` on 9 July 2026 and its outputs are
included.

## Which runs have outputs

Two of the fourteen input files do not correspond to a complete output series,
and we say so rather than let a reader infer otherwise:

- **`cavb_cl05`** (closed state, 5 mm / 10 gap cells) — **never run**. The
  no-barrier series had already shown the secondary parameters to be
  non-asymptotic at this resolution, and the closed-state heat release depends
  on a grid-sensitive extinction sub-model. The input is published for
  completeness; there are no outputs. See §4.3 of the paper.
- **`cavb_op05`** (open state, 5 mm, *original* MPI decomposition) — **diverged**.
  The published CSVs stop at t ≈ 13.5 s and are included as the record of that
  divergence, not as a usable steady-state result. The case was recovered by
  relocating the mesh interfaces clear of the constriction jet; that rerun is
  `cavb_op05if`. See §4.5 and §4.6.5.

The `cavb_op05if` input published here is the as-launched file. The live run
directory additionally carried `&MISC RESTART=.TRUE.` injected by the resume
script after a shutdown; that flag is deliberately **not** published, since a
fresh run must not begin from a restart dump that does not exist.

## Reproducing the paper's numbers

1. Run the input files under FDS 6.10.1 (each case writes `*_devc.csv` / `*_hrr.csv`).
2. `python scripts/gci_stats_v2.py` reproduces every number in Tables 1 and 2
   (steady means, ACF-corrected 95% CIs, two-halves stationarity, variable-r GCI).
3. `python scripts/publication_figures.py` and `python scripts/deconf_jet_figures.py`
   rebuild Figures 1–8 from the raw outputs (paths at the top of each script
   point at the run directory; adjust `RUNS`/`RUN_ROOT` to your local layout).

Simulations were run on the authors' own workstation hardware; no external
funding. Direct `.sf` slice reading is used deliberately (a stale-frame issue
was observed with `fdsreader` on this campaign).

## License

MIT (see `LICENSE`).
