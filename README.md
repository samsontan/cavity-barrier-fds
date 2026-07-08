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
| `outputs/<chid>/` | Device (`*_devc.csv`) and heat-release (`*_hrr.csv`) time series for each completed run |
| `scripts/generate_inputs.py`, `generate_deconf_inputs.py` | Programmatic generators for the input files |
| `scripts/gci_stats_v2.py` | Authoritative statistics: ACF-corrected sampling uncertainty, variable-r Celik GCI, Monte-Carlo p/GCI uncertainty (Tables 1–2) |
| `scripts/gci_analysis.py` | Superseded first-pass GCI script, kept for provenance |
| `scripts/fds_postprocess_cavb.py` | Direct Fortran-binary slice reader + campaign post-processing |
| `scripts/publication_figures.py`, `deconf_jet_figures.py`, `generate_setup_figure.py` | Figure generation (Figures 1–8) |
| `scripts/postprocess_op05if.py` | Outcome classification + steady statistics for the §4.5 disambiguation run |

## Status

The final verification run (`cavb_op05if`, fine-grid open state with relocated
MPI interfaces) is computing at the time of this commit; its output CSVs will
be added on completion. All other runs are final.

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
