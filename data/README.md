# N200 display data

The numbers behind every table, figure, heatmap, sweep plot, distribution plot and trajectory map of the arXiv paper and
of the ablation report (`pdf/ablation_report.pdf`), and behind the display catalogs (`pdf/catalog/`), as flat CSV and
JSON files. N200 is the study with N = 200 Monte Carlo trajectories per configuration (seeds 970001 to 970200); its
summary exports date from 7 October 2026. Nothing was re-simulated: every value is read from the saved summary exports
or from the saved per-trajectory assessments.

Start with `INDEX.md` (one row per display: printed figure or table number, data files, generating script and
command). `COLUMNS.csv` is the column dictionary of every CSV. `VERIFICATION.md` records the check of 50 CSV values
against `metrics.json` (50 of 50 equal). `SHA256SUMS.txt` lists every file (check with `shasum -a 256 -c SHA256SUMS.txt` inside this folder).

The selected settings S1 to S5 of the paper and the report are named T1 to T5 in the data files and New1 to New5 in
figure file names (asset labels); the table data use S1 to S5.

## Contents
| Folder | Files |
|---|---|
| `metrics/` | `{final,new_anchors}/metrics.json.gz`: the summary metrics of the unified grid (`final`) and of the added anchor settings (`new_anchors`): means, SDs, event counts and uncertainty per cell |
| `tables/` | `tables_values.json.gz` (every value of the table fragments with provenance), `tables_values_flat.csv` (same, one row per cell), `tables_values_meta.json` (generator, inputs, rounding rules, fragment lists), `starred_setting.json` (the printed values of the setting S* with their source) |
| `catalog/` | `catalog_values.json.gz` (every value of the catalog table fragments per style), `CATALOG_README.md` (the catalog index: PDFs, order, build) |
| `summary/` | `common_metrics_{final,new_anchors}.csv` (3914 and 720 cells: one row per configuration x method x scenario view; all metric means, SDs, observed counts, trajectory-event counts and one-sided 95% upper bounds, realized-transition events, outcomes, risk-increase checks, assumption flags, assessed populations); `temporal_metrics_{final,new_anchors}.csv` (temporal variance and step RMS) |
| `per_trajectory/` | `per_trajectory_values.csv` (Default, T1-T5 x Corridor, Maze x c=0, c_max, ACI, QT, Unfiltered, V-U, V-Delta, V-L, V-C x 200 seeds: T_mis, m_c/c_max, q+_CBF, q+_unsafe, p-hat^max, h_min, ...), `boxplot_statistics.csv` (quartiles and whiskers as drawn), `source/setting_assessments.json` (the loader cache they come from) |
| `sweeps/` | `sweep_series.csv` (setting, swept parameter alpha/delta/K/b/phi, value, variant or baseline, L_mu, metric, mean, SD, count, upper bound, plotted value); `sweepbox_T3_values.csv` (per-trajectory values behind the arXiv sweep figure `sweepbox_T3_compact.pdf`, 76 boxes x 200 seeds); `sweepbox_ext_values.csv.gz` (per-trajectory values behind the report's sweep figures `sweepbox_<scenario>_<setting>.pdf`, 1024 configurations x 200 seeds, gzip-compressed CSV) |
| `heatmaps/` | `heatmap_grid.csv` (full main grid: scenario, metric, alpha0, alpha, delta, K, variant, mean, SD), `printed_values/*.json` (output of the heatmap builder for the printed heatmaps) |
| `maps/` | the trajectory map manifests: maps_top5 (column and text maps), compact (compact maps), pair (pair maps), report (single maps at the report width) and report pair (the maps the report prints) |
| `tools/` | `build_dump.py`, `make_docs.py`, `verify_sample.py` |

The heatmaps are drawn by the script `code/tools/n200_heatmaps/run_heatmaps.sh`, which writes the paper-style heatmaps to
the folders `out_paper/<name>` and the report-style heatmaps to `out_report/<name>`. The files in
`heatmaps/printed_values/` are named after these folders: `out_<style>__<name>__heatmap_values.json` is the
`heatmap_values.json` of the folder `out_<style>/<name>`. They come from one run of `run_heatmaps.sh` with the scripts in
`code/`, whose heatmaps render identically to the printed ones. The `pdf_sha256` field is the hash of the PDF written in
that run; each PDF records its creation date, so a rerun gives other PDF hashes and otherwise the same files.

## Sources (read only)
- N200 summary exports (not included; the scripts read them from the folder `ECBF_N200_NUMBERS`):
  `final/common/metrics.json` sha256 76d49044...2f76; `final/temporal/temporal_metrics.json` 85ad5851...9e00;
  `new_anchors/common/metrics.json` ea96c428...e5b5e764; `new_anchors/temporal/temporal_metrics.json` 84904d8f...5102.
  `metrics/{final,new_anchors}/metrics.json.gz` hold the two `common/metrics.json` files with every key and value as
  exported, except the free-text descriptions in the field `source_qualification`; their checksums therefore differ
  from those above.
- Per-trajectory values: `per_trajectory/source/setting_assessments.json`, a copy of
  `code/tools/n200_figures/setting_assessments.json`, written by `load_setting_assessments.py`, which checks each saved
  assessment against the trajectory index `N200_INDEX.sqlite` and requires the per-cell means to equal `metrics.json`.
- Tables: `tables/tables_values.json.gz` and `tables/starred_setting.json` are copies of the files of the same name in
  `code/tools/n200_tables/`, written by `build_tables.py` in the table build behind the arXiv paper and the report;
  `tables/tables_values_flat.csv` and `tables/tables_values_meta.json` are derived from that copy by `tools/build_dump.py`.
  Rebuild with `build_tables.py` if the tables change.
- Catalog: `catalog/catalog_values.json.gz` and `catalog/CATALOG_README.md` are copies of
  `code/tools/n200_tables/catalog_values.json.gz` and `code/catalog/README.md`.
- Sweep box plots: `sweeps/sweepbox_T3_values.csv` is written by `code/tools/n200_figures/build_sweeps_box.py`, and
  `sweeps/sweepbox_ext_values.csv.gz` by `tools/build_dump.py` from `code/tools/n200_figures/sweepbox_ext_values.json`,
  the cache of `build_sweeps_box_ext.py` (the same rows as the CSV `sweeps/sweepbox_ext_values.csv` that this script
  writes when it rebuilds its cache). Both builders read the saved per-trajectory assessments through `N200_INDEX.sqlite` and check the cell means,
  event counts and upper bounds against `metrics.json`.
- Maps: the manifests written by the trajectory map builders `build_maps_top5.py` and `build_compact.py` (not included) and by
  `code/tools/n200_maps/build_maps_report.py` and `build_maps_report_pair.py`. In the two report manifests the key
  `v9_builder` names the map builder `build_maps_top5.py`; `reused_generators` in `maps_top5_manifest.json` lists the
  helper scripts of that builder (not included) with their SHA-256.
- Generator hashes: the `generator_sha256` values in the map manifests are those of the scripts as run for the printed
  maps; the scripts in `code/` were edited afterwards, so their hashes differ.

## Conventions
Raw units are kept (fractions, not percent); the plotted percent columns are named `*_pct`, `*_percent` or `plot_value`.
Method names as printed, with the method key of the data files in parentheses where it differs: c=0 (Zero margin),
c_max (Conservative), ACI, QT, Unfiltered, V-U, V-Delta (V-D), V-L (V-Cd), V-C.
Settings, as (alpha0, alpha, delta, K) and paired variant: Default (0.15, 0.10, 0.05, 50) V-C; S1 = T1 (0.30, 0.25, 0.20, 50) V-C;
S2 = T2 (0.30, 0.25, 0.15, 50) V-C; S3 = T3 (0.30, 0.20, 0.20, 50) V-L; S4 = T4 (0.30, 0.25, 0.025, 100) V-Delta;
S5 = T5 (0.30, 0.25, 0.10, 10) V-C. Empty cell = not applicable or undefined (N-A).
Baselines are drawn at the setting's alpha0. Uncertainty: sample SD across trajectories; event bounds are pointwise one-sided 95% (no adjustment for selecting settings).

## Regenerate
Inside `data/`. `make_docs.py` needs only the files of this folder; `build_dump.py` and `verify_sample.py` need the N200
exports in `ECBF_N200_NUMBERS` (see the repository README):
```
python3 tools/build_dump.py      # tables/, catalog/, summary/, per_trajectory/, sweeps/, heatmaps/
python3 tools/make_docs.py       # COLUMNS.csv and INDEX.md
python3 tools/verify_sample.py   # VERIFICATION.md
find . -type f ! -name SHA256SUMS.txt | sed 's|^\./||' | sort | xargs shasum -a 256 > SHA256SUMS.txt
```
`build_dump.py` copies the table, catalog and heatmap files named above from `code/` (the heatmap values only after
`run_heatmaps.sh` has written them) and leaves `maps/` unchanged;
`per_trajectory/source/setting_assessments.json` is copied by hand from `code/tools/n200_figures/`, and
`sweeps/sweepbox_T3_values.csv` is rewritten only by `build_sweeps_box.py --reload` (which writes `S3` in the `setting` column,
where the shipped file has `T3`; the rows are otherwise the same). To redraw a display from the exports use the
command in its `INDEX.md` row (scripts in `code/tools/n200_{figures,tables,heatmaps,maps}`).
Not included: per-trajectory temporal values (they stay in `temporal_metrics.json`), the figure and table PDF/TeX files
themselves, and the raw trajectories.
