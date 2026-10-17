# ecbf-ablation

Code and data for the ablation study of "Anytime-Valid Certification of Robustness Margins for Safe Control via
E-Values" (Trakas, Mirzaeedodangeh, Lindemann, 2026).

## Contents

- `pdf/ablation_report.pdf`: the ablation report, current version (98 pages: implementation details, all sweeps, heatmaps,
  trajectory maps and distributions per setting, the additional S* comparison, and the reference to the arXiv version of the paper).
- `pdf/catalog/`: display catalogs in paper and report style (`display_catalog_*`, `sweeps_{a030,compact,full}_*`),
  built from the paper source as of 12 October 2026; the report above supersedes them where they differ. Their
  captions name files of the paper's LaTeX source (`appendices/`, `sections/`, `figures/`) and the line that includes
  each display; that source is not part of this repository.
- `data/metrics/{final,new_anchors}/metrics.json.gz`: the summary metrics behind every table and figure (N = 200
  trajectories per setting, seeds 970001-970200: means, standard deviations, event counts and uncertainty per cell).
  `final` covers the unified grid, `new_anchors` the added anchor settings.
- `data/`: the display data behind the figures and tables (CSV, JSON, Markdown). Start with `data/README.md`,
  `data/INDEX.md` (one row per display), `data/COLUMNS.csv` (column definitions) and `data/VERIFICATION.md`;
  `data/SHA256SUMS.txt` lists checksums (run `shasum -a 256 -c SHA256SUMS.txt` inside `data/`). Every printed table
  value of the current build is recorded in `code/tools/n200_tables/tables_values.json.gz`; the per-trajectory values
  of the report's sweep box figures are cached in `code/tools/n200_figures/sweepbox_ext_values.json` (as CSV:
  `data/sweeps/sweepbox_ext_values.csv.gz`).
- `code/tools/n200_tables`: table builder, colour switch and checks:
  - `build_tables.py` (all table fragments; `--catalog` also writes the catalog fragments of `catalog_tables.py`),
    `table_colour.py` (coloured and plain variant of each table), `check_tables.py` (re-derives every printed value),
    `compile_tests.py` (compiles every fragment at its target width; test documents in `test/`);
  - `KEY_LINES.md`: the table key, the caption notes and the list of fragments;
  - `starred_setting.json`: the values of the starred row S*;
  - `build_tables_v10_original.py` and `check_tables_v10_original.py`: the earlier table generator and its check,
    which `build_tables.py` extends, kept for reference (only input paths and docstrings were adapted for publication).
- `code/tools/n200_figures`: figure builders and their caches:
  - `extract_sweep_cells.py` (sweep cells from the exports, cache `sweep_cells.json`), `load_setting_assessments.py`
    and `load_default_assessments.py` (per-trajectory values from the saved assessments, caches
    `setting_assessments.json` and `default_assessments.json`);
  - `build_sweeps_box.py` (box-plot sweep of S3 in the paper), `build_sweeps_box_ext.py` (box-plot sweeps of every
    setting in the report), `build_sweeps.py` and `build_sweeps_twin.py` (line-plot sweeps),
    `build_distribution_figures.py` (per-trajectory distributions), `build_theorem_checks_box.py` (theorem checks),
    `build_default_figures.py` (default-setting figures); `n200_style.py` holds the shared styles.
- `code/tools/n200_maps`: report trajectory maps (`build_maps_report.py`, `build_maps_report_pair.py`; the report
  prints the pair maps, `data/maps/maps_report_pair_manifest_*.json`).
- `code/tools/n200_heatmaps` (`build_heatmaps.py`, `run_B4.sh` for the printed heatmaps), `code/tools/n200_catalog`
  (copies of the builders with the catalog variants, `run_catalog_heatmaps.sh`) and `code/catalog/build_catalog.py`:
  heatmap and catalog builders.
- `code/simulation`: the simulation package `ecbf` that generated the trajectories (dynamics and nominal controllers,
  CBF filter, e-value tests and certification, certificate bounds of V-U, V-Delta, V-L and V-C, geometry, risk evaluation
  and assessment) and `run_example.py`, which runs one trajectory. Section 1 of the report gives the implementation
  details this code follows: quadrature, interval refinement, the rectangle condition and the Jacobian bounds.
  `code/simulation/README.md` maps each module to the paper's terms.

The selected settings S1-S5 of the paper and the report carry the labels T1-T5 (asset labels New1-New5) in the data
files and scripts.

## Rebuilding

Requirements: Python 3.12 with numpy and matplotlib (`code/simulation` also needs scipy), a LaTeX installation with
`pdflatex` and `latexmk` (all figure text is set with `text.usetex`), poppler (`pdfinfo`, used to check figure widths)
and bash for the two heatmap scripts (`PYTHON=...` selects their interpreter, default `python3`).

Inputs that are not part of this repository are named by environment variables. Each defaults to the file or folder
named below at the repository root (computed from the script's location under `code/` or `data/`):

| Variable | Default | Content |
|---|---|---|
| `ECBF_N200_NUMBERS` | `N200_numbers/` | summary exports of the N = 200 runs (not included): `{final,new_anchors}/common/metrics.json`, `{final,new_anchors}/common/population.jsonl`, `{final,new_anchors}/temporal/temporal_metrics.json` |
| `ECBF_N200_DATA` | `N200_data/` | trajectory index `N200_INDEX.sqlite` and the saved per-trajectory assessments it points to (not included) |
| `ECBF_N100_CATALOG` | `N100_data/runs.sqlite` | run catalog of the sealed N = 100 campaign (source `n100_prefix_sealed` of the index; not included) |
| `ECBF_CANDIDATES` | `candidate_review.json` | candidate review of the unified grid, read by the table builder (not included) |
| `ECBF_FIGURE_STYLE` | `figure_style/` | figure helper `ieee_figure.py` (strict save: vector PDF at the final width, text-size and line-width checks; not included) |
| `ECBF_MAPS_BUILDER` | `maps_top5/build_maps_top5.py` | trajectory map builder (not included) |
| `ECBF_DISPLAYS` | `data/` | display data folder of this repository (the box-plot builders write their CSV files there; `build_sweeps_box_ext.py` writes `sweeps/sweepbox_ext_values.csv` uncompressed) |
| `ECBF_PREVIOUS_SETUP` | `appendices/f_tab_setup_previous.tex` of the paper source tree | optional: an earlier version of the parameter table (arXiv Table III); when the file exists, `check_tables.py` also checks that its rows are unchanged, otherwise it skips that comparison with a note |

The figure, heatmap and map builders run in place in this repository and write to the folder given by `--out`. The
table and catalog builders need the paper's LaTeX source tree (not included):

1. Copy `code/tools` and `code/catalog` into the paper source tree, next to `main.tex` and `figures/`. The table,
   figure and catalog scripts write their fragments, figures and PDFs relative to that tree. After copying, the
   defaults above resolve relative to the copy, so set every `ECBF_*` variable, with `ECBF_DISPLAYS` pointing to the
   `data/` folder of this repository.
2. Tables: `python3 tools/n200_tables/build_tables.py`, then `python3 tools/n200_tables/check_tables.py` (must report
   0 errors) and `python3 tools/n200_tables/compile_tests.py` (compiles the test documents in both styles).
   `build_tables.py` reads the exports in `ECBF_N200_NUMBERS` and records every printed value in
   `tools/n200_tables/tables_values.json.gz`; `check_tables.py` re-derives each recorded value from the exports.
3. Figures (`--kind paper` or `--kind report`): in `tools/n200_figures`, refresh the caches with
   `extract_sweep_cells.py`, `load_setting_assessments.py` and `load_default_assessments.py`, then run
   `build_sweeps_box.py --omit-zero-rows --out DIR`, `build_sweeps_box_ext.py --kind KIND --out DIR`,
   `build_distribution_figures.py --kind KIND --out DIR`, `build_theorem_checks_box.py --kind KIND --out DIR` and
   `build_default_figures.py --kind KIND --out DIR`; they write vector PDFs at the final printed width and a manifest
   per run. Heatmaps: `tools/n200_heatmaps/run_B4.sh` (or `build_heatmaps.py --help` for single figures).
4. Trajectory maps need the raw trajectory payloads (about 54 GB, not included) and the map builder
   `build_maps_top5.py` (not included, `ECBF_MAPS_BUILDER`); `tools/n200_maps/build_maps_report*.py` call it with the
   report width and text size.
5. The report itself is compiled from the paper source by its `build.py` (not included).

Single trajectories are rerun with `code/simulation/run_example.py` (Python 3.12 with numpy and scipy) for a given
scenario, variant, $(\alpha_0, \alpha, \delta, K)$, budget normalizer $b$, bound factor $\phi$, $L_\mu$ and seed. The
study uses the seeds 970001 to 970200 ($N_{\mathrm{MC}} = 200$ per setting). The campaign tooling that scheduled the full
set of runs is not included.

## Data provenance

The data are the summary values of the N = 200 Monte Carlo trajectories per setting, exported on 7 October 2026 from the
simulation results (3,914 cells in `final`, 720 in `new_anchors`). `data/VERIFICATION.md` compares a sample of 50
CSV values with the exports (`data/tools/verify_sample.py`, exact equality), and `check_tables.py` re-derives every
printed table value from the exports. Raw trajectories are not included; `data/per_trajectory` holds the per-trajectory summaries
used by the displays.

## Paper

The paper links to this repository and cites the ablation report. The arXiv paper reports event frequencies and their
one-sided 95% upper confidence bounds in percent, and mean decision fractions as percentages with the sample SD. The
ablation report gives event counts out of N = 200; the certificate number is a count per trajectory.

No license is granted; contact the authors for reuse.
