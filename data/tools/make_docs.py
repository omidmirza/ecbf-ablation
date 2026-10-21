#!/usr/bin/env python3
"""Writes COLUMNS.csv (column dictionary of every CSV) and INDEX.md (one row per display). Run after build_dump.py.
Command (inside data/):  python3 tools/make_docs.py"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]

CFG = {
    "source": "export folder: final (N200_numbers/final) or new_anchors",
    "view_id": "N200 view id (SHA-256) of the declared metric view in metrics.json",
    "family": "view family: main_grid_200, stress_drift_200 (phi x L_mu), budget_normalizer_200 (b), detailed_200, matched_baselines_200, matched_ACI_200",
    "view_label": "label of the view as recorded in the export ('Stress factor and disturbance drift': the study of the bound factor phi and the distribution-drift bound L_mu)",
    "anchor": "role of a new_anchors view: new_top5_role1..5 = T1..T5 (= S1..S5); empty in final",
    "anchor_id": "anchor setting of the view: current_default = default setting; n200_role1..5_* = S1..S5 (new_anchors); candidate_1..candidate_5 (final) = other candidate settings, not S1..S5, with (alpha0, alpha, delta, K) = (0.30, 0.25, 0.20, 10), (0.30, 0.25, 0.15, 10), (0.30, 0.20, 0.20, 10), (0.30, 0.25, 0.10, 100), (0.30, 0.25, 0.05, 50), each also at K = infinity",
    "role": "cohort role recorded by the export: development for every view (the N200 results are descriptive)", "anchor_alpha0": "alpha0 of the anchor setting (new_anchors)",
    "anchor_alpha": "alpha of the anchor setting", "anchor_delta": "delta of the anchor setting", "anchor_K": "K of the anchor setting",
    "n100_view_id": "view id of the same cell in the N100 export (seeds 970001-970100, the first half of the N200 seeds)", "seed_first": "first seed (970001)",
    "seed_last": "last seed (970200)", "n_seeds": "number of seeds = trajectories (200)",
    "scenario": "corridor or maze", "method": "internal method key (Zero margin=c=0, Conservative=c_max, ACI, QT, Unfiltered, V-U, V-D=V-Delta, V-Cd=V-L, V-C)",
    "method_display": "name printed in the documents", "alpha0": "alpha_0, CBF-violation probability (the baselines run at the same alpha_0)", "alpha": "alpha, margin insufficiency (variants)",
    "delta": "delta, error confidence (variants)", "K": "K (test period in decisions; infinity = never re-test)",
    "L_mu": "L_mu, distribution-drift bound", "phi_bound_scale": "phi (bound factor, configuration key bound_scale)",
    "b_budget_normalizer": "b (budget normalizer, key budget_normalizer; 2.11 default)", "rho_buf": "rho_buf (0)",
    "configuration_id": "SHA-256 id of the configuration", "population_sha256": "SHA-256 of the assessed seed population", "n": "trajectories in the cell (200)",
    "comparison_eligible": "export flag: cell may be compared as a paired view",
}
STAT = {"total_n": "trajectories", "observed_n": "trajectories with a value", "missing_n": "trajectories without a value",
        "missing_reasons": "JSON counts of missing reasons", "applicable": "metric defined for this method", "conditional": "mean is conditional on availability",
        "mean": "sample mean over observed trajectories", "sd": "sample SD (n-1) over observed trajectories", "available_mean": "same as mean (kept from export)",
        "available_sd": "same as sd (kept from export)", "comparison_eligible": "usable in paired comparisons"}
METRIC = {"T_mis": "mission time (s)", "m_c": "mean robustness margin / c_max", "m_Y": "mean barrier loss / c_max", "m_int": "control intervention (stored input coordinates squared)",
          "q_CBF": "fraction of executed decisions with a CBF violation (print x100 for %)", "q_Y": "fraction of barrier-loss exceedance",
          "q_use": "fraction of decisions using a smaller margin", "q_unsafe": "fraction with a safe-set violation (print x100 for %)",
          "q_nocert": "fraction of decisions with no certificate", "q_notrans": "no-deployable-margin fraction q+_nodep (a certified candidate but no deployable margin; report eq. (32))",
          "tau_use": "mean time of first smaller-margin use (s)", "p_max": "mean of the per-trajectory peak estimated exceedance risk",
          "runtime_ms": "controller wall time per decision (ms; descriptive)", "h_min": "mean minimum barrier value (m)",
          "clearance_min": "mean minimum clearance (m)", "n_certificates": "mean number of certificates"}
EVENT = {"applicable_n": "trajectories where the event is defined", "confirmed": "trajectories with a confirmed event",
         "possible": "trajectories with a confirmed or numerically undecided event (the printed count)", "unresolved_trajectories": "undecided trajectories",
         "confirmed_fraction": "confirmed / applicable_n", "possible_fraction": "possible / applicable_n",
         "upper95_possible": "one-sided 95% Clopper-Pearson upper bound for the possible count (fraction)",
         "lower95_confirmed": "one-sided 95% lower bound for the confirmed count (fraction)", "bound": "text of the bound rule",
         "total_n": "trajectories", "comparison_eligible": "usable in paired comparisons"}


def describe_common(c):
    if c in CFG:
        return CFG[c]
    m = re.match(r"statistics\.(\w+)\.(\w+)$", c)
    if m:
        return f"{METRIC.get(m[1], m[1])}: {STAT.get(m[2], m[2])}"
    m = re.match(r"trajectory_events\.(qcert|qdep)\.(\w+)$", c)
    if m:
        return f"{'q_cert (false-certification event, report eq. (27))' if m[1] == 'qcert' else 'q_dep (deployment event: p_t(c_t) > alpha_0 at some decision time, report eq. (27))'} trajectory events: {EVENT.get(m[2], m[2])}"
    m = re.match(r"realized_transition_events\.(\w+)\.(\w+)$", c)
    if m:
        return f"realized-transition event {m[1]}: {m[2]} (trajectories with at least one observed event on executed transitions; fraction, upper95 as fractions)"
    m = re.match(r"deployment_decision_fraction\.(confirmed|possible)\.(\w+)$", c)
    if m:
        return (f"fraction of decisions with p_t(c_t) > alpha_0 (deployment-event decisions; {m[1]}: "
                f"{'lower' if m[1] == 'confirmed' else 'upper'} risk bound above alpha_0): {STAT.get(m[2], m[2])}")
    m = re.match(r"outcomes\.(\w+)\.(\w+)$", c)
    if m:
        return f"termination outcome {m[1]}: {m[2]} (count of trajectories; fraction, lower95/upper95 binomial bounds)"
    m = re.match(r"risk_increase_checks\.(direct|pairwise)\.(\w+)$", c)
    if m:
        return f"risk-increase bound check ({m[1]}): {m[2]} (decisions / trajectories; is_proof=false: numerical check, not a proof)"
    if c.startswith("assumption_observed_flags."):
        return "observed counts of assumption-violation flags: " + c.split(".", 1)[1]
    if c.startswith("numerical_ambiguity."):
        return "counts of numerically unresolved items: " + c.split(".", 1)[1]
    if c.startswith("assumption_assessed_population."):
        return "population size over which the assumption flags were assessed: " + c.split(".", 1)[1]
    return {"termination_counts": "JSON counts of trajectory terminations (goal, contact, infeas)", "source_counts": "JSON counts of the run sources of the 200 seeds",
            "runtime_scope": "text: runtime is descriptive", "scope": "text: scope of the cohort", "deployment_decision_fraction.scope": "text",
            "deployment_decision_fraction.resolved": "number of decisions with resolved deployment interval"}.get(c, c)


def describe_temporal(c):
    if c in CFG:
        return CFG[c]
    m = re.match(r"statistics\.(normalized_margin|estimated_applied_risk)\.(temporal_variance|step_rms)\.(\w+)$", c)
    if m:
        q = "normalized margin" if m[1] == "normalized_margin" else "estimated applied risk"
        mt = "temporal variance" if m[2] == "temporal_variance" else "step RMS (root mean square of consecutive differences)"
        return f"{q}, {mt} over the first 200 decisions, across trajectories: {m[3]} (n, mean, sample SD, min, max)"
    m = re.match(r"statistics\.(normalized_margin|estimated_applied_risk)\.(\w+)$", c)
    if m:
        return f"{m[1]}: {m[2]} (assessed trajectories / missing reasons)"
    return {"statistics.total_n": "trajectories", "statistics.prefix_decisions": "decisions in the fixed prefix (200)",
            "statistics.population": "population note"}.get(c, c)


STARRED = " (starred row S*)"
FILES = {
    "tables/tables_values_flat.csv": {
        "target": "fragment set: paper (arXiv paper and other paper-style fragments) or report (ablation report)",
        "fragment": "table fragment (file name generated_n200/{paper,report}/<fragment>.tex in the paper source tree)",
        "table": "table environment id inside the fragment", "method": "variant / baseline of the row (internal key)", "spec": "metric or quantity printed",
        "printed": "text exactly as printed in the cell (LaTeX)",
        "kind": "cell kind (meanSD, tmeanSD, meanSD_plain, event, event_pct, event_pct_sweep, event_count_compact, count, NA, NA_dash, design_values, a6_total, starred, starred_variant, starred_compact, starred_missing)",
        "fixed": "design_overview rows: fixed parameters printed in the row (LaTeX)",
        "sets": "design_overview rows: JSON of the compared values (alpha per alpha0, delta, K, b, phi, L_mu, methods)",
        "row": "JSON of the row key", "row_label": "row key as text",
        "src": "export the cell comes from: final, new (new_anchors), final_temporal, or the starred-row sources (starred_variant, starred_setting.json)",
        "view_id": "N200 view id", "n100_view_id": "view id of the same cell in the N100 export (seeds 970001-970100, the first half of the N200 seeds)",
        "family": "view family", "configuration_id": "configuration id", "configuration": "JSON configuration",
        "source": "key of the value in the export record (for example statistics.m_c, trajectory_events.qdep)",
        "scale": "display scale (value x scale printed)", "decimals": "printed decimals", "mean_str": "printed mean", "sd_str": "printed SD",
        "conditional": "conditional mean", "boldable": "may be bold", "arrow": "preferred direction", "bold": "printed bold",
        "assessed": "trajectories assessed", "count": "event count printed", "n": "denominator", "upper95": "upper bound used (fraction)",
        "upper_str": "printed bound (%, rounded up)", "check": "check mark printed", "certified": "row is a variant, so the guarantee check applies (false for baselines)",
        "level": "guarantee level used for the check", "counts": "event counts (confirmed, possible)",
        "count_rule": "printed count when confirmed and possible counts differ: range (both), confirmed (lower end) or upper (upper end)",
        "confirmed": "confirmed event count (cells with count_rule)", "pct_str": "printed event frequency (%)",
        "stat_keys": "JSON: statistics fields read from the export (mean and sd, or available_mean and available_sd)",
        "source_file": "source table file of the cell in the paper source tree, not included" + STARRED,
        "source_line": "line in the source table file" + STARRED, "json_spec": "spec key in tables/starred_setting.json" + STARRED,
        "count_str": "printed count", "compact_mean": "mean printed in the compact form" + STARRED,
        "source_text": "text of the cell in the source table" + STARRED, "factor": "factor applied to the source value" + STARRED,
        "source_column": "column in the source table file" + STARRED, "sort_key": "value used to choose the bold entry" + STARRED,
        "digits": "rounded (re-rounded to the printed decimals) or source (digits of the source table)" + STARRED,
        "na": "cell printed as not applicable (N-A)",
        "name": "assumption_checks_arxiv: LaTeX macro of the printed number", "field": "assumption_checks_arxiv: quantity counted",
        "value": "assumption_checks_arxiv: value of the number", "scope": "assumption_checks_arxiv: configurations and views counted"},
    "per_trajectory/per_trajectory_values.csv": {
        "setting": "Default or T1..T5 (= S1..S5 in the paper and the report)", "asset_label": "file label Default, New1..New5", "scenario": "corridor or maze", "method": "method name printed",
        "method_key": "internal key", "is_baseline": "baseline (c=0, c_max, ACI, QT, Unfiltered) at the setting's alpha0", "setting_alpha0": "alpha0 of the setting",
        "setting_alpha": "alpha of the setting", "setting_delta": "delta of the setting", "setting_K": "K of the setting", "setting_paired_variant": "paired variant of the setting",
        "L_mu": "L_mu (0.0005)", "seed": "trajectory seed 970001..970200", "T_mis": "mission time T_mis,k (s)", "m_c": "m_c,k / c_max", "m_Y": "m_Y,k / c_max",
        "q_CBF": "q+_CBF,k fraction of executed decisions", "q_unsafe": "q+_unsafe,k fraction", "q_use": "q+_use,k fraction", "q_nocert": "q+_nocert,k fraction",
        "q_notrans": "q+_nodep,k no-deployable-margin fraction (report eq. (32))", "p_max": "estimated peak exceedance risk p-hat_k^max (empty = N-A)", "h_min": "min_t h (m)",
        "tau_use": "time of first smaller-margin use (s; empty if none)", "N": "executed decisions", "q_CBF_pct": "q+_CBF,k in % (as plotted)", "q_unsafe_pct": "q+_unsafe,k in % (as plotted)"},
    "per_trajectory/boxplot_statistics.csv": {
        "setting": "Default or T1..T5 (= S1..S5)", "asset_label": "file label", "scenario": "scenario", "method": "method", "metric": "T_mis, m_c, q_CBF (%), q_unsafe (%), p_max, h_min",
        "plot_unit": "percent for q rows", "n": "values", "mean": "mean (plot units)", "sd": "sample SD", "min": "minimum", "q1": "first quartile (linear interpolation)",
        "median": "median", "q3": "third quartile", "max": "maximum", "whisker_low": "lowest value >= q1 - 1.5 IQR", "whisker_high": "highest value <= q3 + 1.5 IQR"},
    "sweeps/sweep_series.csv": {
        "setting": "Default or T1..T5 (= S1..S5)", "asset_label": "file label", "alpha0": "setting alpha0", "alpha": "setting alpha", "delta": "setting delta", "K": "setting K",
        "paired_variant": "paired variant of the setting", "scenario": "scenario", "swept_parameter": "alpha, delta, K, b, phi, or baseline (horizontal reference lines)",
        "x_value": "value of the swept parameter (K as 10, 50, 100, 200, infinity)", "x_position": "plotted x position (K: index 0..4)", "variant": "variant or baseline name printed",
        "method_key": "internal key", "L_mu": "L_mu of the cell (phi column: 0.0005 solid, 0.005 dashed)", "is_setting_value": "x equals the setting's own value (vertical dotted line)",
        "source": "final or new_anchors", "view_id": "N200 view id", "metric": "m_c, q_CBF, q_unsafe, q_use, T_mis, p_max, h_min, q_cert, q_dep",
        "kind": "statistic (mean +- SD) or event_count (trajectory count)", "mean": "sample mean (raw units, fractions)", "sd": "sample SD", "count": "q_cert / q_dep: possible (confirmed + undecided) count of 200",
        "count_confirmed": "confirmed count", "upper95": "q_cert / q_dep: one-sided 95% upper bound of the possible count (fraction)", "n_assessed": "trajectories assessed",
        "n_total": "trajectories (200)", "applicable": "metric defined", "plot_value": "value as plotted (q_CBF and q_unsafe in %)", "plot_unit": "unit of plot_value"},
    "sweeps/sweepbox_T3_values.csv": {
        "setting": "T3 (= S3; build_sweeps_box.py writes S3 when it rewrites this file)", "asset_label": "file label New3", "alpha0": "setting alpha0", "alpha": "setting alpha", "delta": "setting delta", "K": "setting K",
        "scenario": "corridor or maze", "swept_parameter": "alpha, delta, K or b", "x_value": "value of the swept parameter (K as 10, 50, 100, 200, infinity)",
        "x_position": "plotted x position", "is_setting_value": "x equals the setting's own value (grey column band)", "variant": "V-L or V-C",
        "method_key": "internal key (V-Cd = V-L)", "source": "final or new_anchors", "family": "view family of the cell", "view_id": "N200 view id",
        "configuration_id": "SHA-256 id of the configuration", "seed": "trajectory seed 970001..970200",
        "index_source": "run source of the trajectory in the trajectory index", "assessment_sha256": "SHA-256 of the saved per-trajectory assessment",
        "m_c_over_cmax": "m_c,k / c_max", "q_CBF_fraction": "q+_CBF,k (fraction)", "q_CBF_percent": "q+_CBF,k (%, as plotted)", "T_mis_s": "T_mis,k (s)",
        "qdep_applicable": "q_dep event defined for this trajectory", "qdep_confirmed": "confirmed q_dep event", "qdep_possible": "confirmed or numerically undecided q_dep event",
        "cell_mean_m_c": "cell mean of m_c,k / c_max (equals metrics.json)", "cell_mean_q_CBF_fraction": "cell mean of q+_CBF,k (equals metrics.json)",
        "cell_mean_T_mis_s": "cell mean of T_mis,k (equals metrics.json)", "cell_qdep_possible_count": "possible q_dep count of the cell",
        "cell_qdep_confirmed_count": "confirmed q_dep count of the cell", "cell_qdep_upper95": "one-sided 95% Clopper-Pearson upper bound of the possible count (fraction)"},
    "sweeps/sweepbox_ext_values.csv.gz": {
        "configuration_id": "SHA-256 id of the configuration", "alpha0": "alpha0", "alpha": "alpha", "delta": "delta", "K": "K (infinity = never re-test)",
        "budget_normalizer": "b (budget normalizer)", "bound_scale": "phi (bound factor)", "L_mu": "L_mu (drift constant)", "benchmark": "scenario: corridor or maze",
        "method": "internal method key (V-U, V-D = V-Delta, V-Cd = V-L, V-C)", "seed": "trajectory seed 970001..970200",
        "index_source": "run source of the trajectory in the trajectory index", "assessment_sha256": "SHA-256 of the saved per-trajectory assessment",
        "m_c_over_cmax": "m_c,k / c_max", "q_CBF_fraction": "q+_CBF,k (fraction)", "q_unsafe_fraction": "q+_unsafe,k (fraction)", "T_mis_s": "T_mis,k (s)",
        "qcert_applicable": "q_cert event defined for this trajectory", "qcert_confirmed": "confirmed false certification", "qcert_possible": "confirmed or numerically undecided false certification",
        "qdep_applicable": "q_dep event defined for this trajectory", "qdep_confirmed": "confirmed q_dep event", "qdep_possible": "confirmed or numerically undecided q_dep event"},
    "heatmaps/heatmap_grid.csv": {
        "scenario": "scenario", "metric": "metric key (see summary dictionary)", "alpha0": "alpha0", "alpha": "alpha", "delta": "delta", "K": "K", "variant": "variant printed",
        "method_key": "internal key", "mean": "sample mean (raw units; heatmaps print mean x display scale)", "sd": "sample SD", "n_assessed": "trajectories with a value",
        "n_total": "trajectories (200)", "applicable": "metric defined", "conditional": "conditional mean", "view_id": "N200 view id"},
}


def header(f):
    opener = gzip.open if f.endswith(".gz") else open
    with opener(OUT / f, "rt", newline="") as fh:
        return next(csv.reader(fh))


def main():
    rows = []
    for f, d in FILES.items():
        for c in header(f):
            rows.append((f, c, d.get(c, "?")))
    for f, fn in (("summary/common_metrics_final.csv", describe_common), ("summary/common_metrics_new_anchors.csv", describe_common),
                  ("summary/temporal_metrics_final.csv", describe_temporal), ("summary/temporal_metrics_new_anchors.csv", describe_temporal)):
        rows += [(f, c, fn(c)) for c in header(f)]
    unknown = [r for r in rows if r[2] == "?"]
    if unknown:
        raise SystemExit(f"undocumented columns: {unknown[:10]}")
    with open(OUT / "COLUMNS.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["file", "column", "description"])
        w.writerows(rows)
    print("COLUMNS.csv", len(rows))
    index()


# ----------------------------------------------------------------------------------------------- INDEX.md
SRC = "code"        # the scripts as kept in this repository (code/tools/...)
PY = "python3"
SC = ("corridor", "maze")
KEYS = ("Default", "New1", "New2", "New3", "New4", "New5")
LAB = {"Default": "default", "New1": "T1", "New2": "T2", "New3": "T3", "New4": "T4", "New5": "T5"}
SNAME = {"Default": "default", "New1": "S1", "New2": "S2", "New3": "S3", "New4": "S4", "New5": "S5"}
SET = {"Default": "(0.15, 0.10, 0.05, 50) V-C", "New1": "(0.30, 0.25, 0.20, 50) V-C", "New2": "(0.30, 0.25, 0.15, 50) V-C",
       "New3": "(0.30, 0.20, 0.20, 50) V-L", "New4": "(0.30, 0.25, 0.025, 100) V-Delta", "New5": "(0.30, 0.25, 0.10, 10) V-C"}
FIG = f"cd {SRC}/tools/n200_figures && {PY}"
HEAT = f"cd {SRC}/tools/n200_heatmaps && ./run_heatmaps.sh"
MAPS_EXT = "map builder `build_maps_top5.py` (not included, `ECBF_MAPS_BUILDER`; needs the raw trajectory payloads)"
# printed figure numbers of the ablation report (pdf/ablation_report.pdf)
R_MAP = {"corridor": {"Default": 1, "New3": 2, "New1": 3, "New2": 4, "New4": 5, "New5": 6},
         "maze": {"Default": 7, "New3": 8, "New1": 9, "New2": 10, "New4": 11, "New5": 12}}
R_DIST = {"Default": 13, "New3": 14, "New1": 15, "New2": 16, "New4": 17, "New5": 18}
R_SWEEP = {"corridor": {"Default": 19, "New3": 21, "New1": 23, "New2": 25, "New4": 27, "New5": 29},
           "maze": {"Default": 20, "New3": 22, "New1": 24, "New2": 26, "New4": 28, "New5": 30}}
R_CHECKS = {"corridor": 42, "maze": 43, "Default": 44}


def setting_filter(k):
    return "Default" if k == "Default" else LAB[k]


def index():
    L = []
    add = L.append
    add("# Index of displays\n")
    add("One row per figure and table of the arXiv paper and of the ablation report (`pdf/ablation_report.pdf`), with the data "
        "files in this folder and the script that draws it. Paths in the data columns are relative to this folder (`data/`); "
        "`code/tools/...` paths are relative to the repository root; `appendices/...` and `figures/...` paths are relative to "
        "the paper source tree (not included). The column Printed as gives the figure or table number in the arXiv paper "
        "and in the report. N200 is the study with N = 200 trajectories per configuration.\n")
    add("Settings: " + "; ".join(f"{SNAME[k]} {v}" for k, v in SET.items()) + ", as (alpha0, alpha, delta, K) and paired variant. "
        "S1..S5 of the paper and the report are named T1..T5 in the data files and New1..New5 in figure file names "
        "(the table data use S1..S5).\n")
    add("The CSV files were written by `tools/build_dump.py` from the N200 summary exports, and `COLUMNS.csv` and this index by "
        "`tools/make_docs.py`. A display is redrawn from the exports only by its script in `code/tools/...` (no simulation). "
        "The figure and heatmap commands below run inside this repository and need the figure helper `ieee_figure.py` "
        "(`ECBF_FIGURE_STYLE`, not included). With it, the distribution and sweep figures run from the caches in "
        "`code/tools/n200_figures`; the theorem-check figures and the heatmaps also need the exports (`ECBF_N200_NUMBERS`), "
        "and refreshing the caches needs the saved assessments (`ECBF_N200_DATA`). The table builder runs in the paper source "
        "tree (repository README, Rebuilding). The scripts write their PDFs to the `--out` folder given.\n")

    add("\n## Figures\n")
    add("| Display (file) | Printed as | Data file(s) here | Script and command |")
    add("|---|---|---|---|")
    # trajectory maps (report: pair maps at the report width)
    for s in SC:
        for k in KEYS:
            pr = f"report Fig. {R_MAP[s][k]}"
            if s == "corridor" and k == "New3":
                pr = "arXiv Fig. 4, " + pr
            add(f"| `map_{s}_{k}_pair.pdf` ({s.title()}, {SNAME[k]}) | {pr} | "
                f"`maps/maps_report_pair_manifest_Default_New1_New2_New3_New4_New5.json` (per-figure values and hashes at the report width); "
                f"seeds and per-seed records: `maps/pair_manifest.json` | "
                f"`{SRC}/tools/n200_maps/build_maps_report_pair.py` (calls the {MAPS_EXT}) |")
    add(f"| `map_maze_New3_pair_zoom.pdf` (Maze, S3, enlarged view) | arXiv Fig. 5 | none (trajectory map; no manifest in this folder) | {MAPS_EXT} |")
    add(f"| `map_maze_New3_compact_main_abcd.pdf` (Maze, S3 with V-L: paths, margin, loss and risk of one trajectory) | "
        f"arXiv Fig. 3 (Sec. VI) | none (time series of one trajectory; raw trajectories not included) | {MAPS_EXT} |")
    # distributions
    for k in KEYS:
        pr = ("arXiv Fig. 6, " if k == "New3" else "") + f"report Fig. {R_DIST[k]}"
        add(f"| `distributions_{k}_half.pdf` (Corridor and Maze, {SNAME[k]}) | {pr} | "
            f"`per_trajectory/per_trajectory_values.csv` (setting={setting_filter(k)}), `per_trajectory/boxplot_statistics.csv` (setting={setting_filter(k)}) | "
            f"`{SRC}/tools/n200_figures/build_distribution_figures.py`: `{FIG} build_distribution_figures.py --kind report --out out_report "
            f"--form half --labels {k}` (arXiv: `--kind paper --out out_paper`); cache `setting_assessments.json` from `load_setting_assessments.py` |")
    # heatmaps, in the order of the report (Figs. 31-41)
    heat = [("corridor", "q_CBF", "0p15", "report/a"), ("maze", "q_CBF", "0p15", "report/a"),
            ("corridor", "m_c", "0p15", "report/a"), ("maze", "m_c", "0p15", "report/a"),
            ("corridor", "m_c", "0p3", "report/a"), ("corridor", "q_CBF", "0p3", "report/a"),
            ("corridor", "q_unsafe", "0p3", "report/b"), ("maze", "m_c", "0p3", "report/a"),
            ("maze", "q_CBF", "0p3", "report/a"), ("maze", "q_unsafe", "0p3", "report/b"),
            ("corridor", "q_unsafe", "0p25", "report/c")]
    for i, (s, m, a0, folder) in enumerate(heat):
        pv = "out_" + folder.replace("/", "__") + "__heatmap_values.json"
        add(f"| `heatmap_{s}_{m}_a0_{a0}_d_all.pdf` | report Fig. {31 + i} | `heatmaps/heatmap_grid.csv` (scenario={s}, metric={m}, "
            f"alpha0={a0.replace('p', '.')}; all alpha, delta, K, variants); printed values: `heatmaps/printed_values/{pv}` | "
            f"`{SRC}/tools/n200_heatmaps/build_heatmaps.py`: `{HEAT}` (output folder `out_{folder}`) |")
    add(f"| `heatmap_maze_m_c_a0_0p3_d_0p05-0p1.pdf` (Maze, alpha0 = 0.30, delta = 0.05 and 0.10) | arXiv Fig. 7 | "
        f"`heatmaps/heatmap_grid.csv` (scenario=maze, metric=m_c, alpha0=0.3, delta in 0.05, 0.1); printed values: "
        f"`heatmaps/printed_values/out_paper__arxiv_maze__heatmap_values.json` | `{HEAT}` (output folder `out_paper/arxiv_maze`) |")
    # sweeps
    for s in SC:
        for k in ("Default", "New3", "New1", "New2", "New4", "New5"):
            add(f"| `sweepbox_{s}_{k}.pdf` ({s.title()}, {SNAME[k]}: one-parameter sweeps) | report Fig. {R_SWEEP[s][k]} | "
                f"`sweeps/sweepbox_ext_values.csv.gz` (rows of the sweeps around the setting: filter on the configuration columns and "
                f"benchmark={s}); QT and ACI: `sweeps/sweep_series.csv` (setting={setting_filter(k)}, scenario={s}, swept_parameter=baseline) | "
                f"`{SRC}/tools/n200_figures/build_sweeps_box_ext.py`: `{FIG} build_sweeps_box_ext.py --kind report --out out_report "
                f"--labels {k} --scens {s}`; cache `sweepbox_ext_values.json` |")
    add(f"| `sweepbox_T3_compact.pdf` (Corridor and Maze, S3, V-L and V-C: sweeps of alpha, delta, K, b) | arXiv Fig. 8 | "
        f"`sweeps/sweepbox_T3_values.csv` (per-trajectory values, 76 boxes x 200 seeds); QT and ACI: `sweeps/sweep_series.csv` (setting=T3, "
        f"swept_parameter=baseline) | `{SRC}/tools/n200_figures/build_sweeps_box.py`: `{FIG} build_sweeps_box.py --omit-zero-rows --out out_paper` "
        f"(`--reload` re-reads the assessments and repeats the checks); cache `sweepbox_T3_values.json` |")
    # theorem checks
    for s in SC:
        add(f"| `theorem_checks_{s}_New3.pdf` (S3 with V-L: event frequencies against the guarantee levels) | report Fig. {R_CHECKS[s]} | "
            f"`summary/common_metrics_new_anchors.csv` (family=detailed_200, anchor=new_top5_role3, scenario={s}: `trajectory_events.*`, "
            f"`risk_increase_checks.*`) | `{SRC}/tools/n200_figures/build_theorem_checks_box.py`: `{FIG} build_theorem_checks_box.py --kind report --out out_report` |")
    add(f"| `theorem_checks_corridor_Default.pdf` (default setting with V-C) | report Fig. {R_CHECKS['Default']} | "
        f"`summary/common_metrics_final.csv` (family=detailed_200, anchor_id=current_default, scenario=corridor: `trajectory_events.*`) | "
        f"`{FIG} build_theorem_checks_box.py --kind report --out out_report --setting Default` |")

    add("\nOther outputs of the scripts, not printed in the arXiv paper or the report: `sweepplot_<scen>_<label>.pdf` and "
        "`sweepplot_T3_compact.pdf` (`build_sweeps.py` after `extract_sweep_cells.py`; data `sweeps/sweep_series.csv`), "
        "`sweepplot_T3_compact_twin.pdf` (`build_sweeps_twin.py`), `distributions_<scen>_<label>.pdf` "
        "(`build_distribution_figures.py --form separate`), the remaining heatmaps of `run_heatmaps.sh` (output folders "
        "`out_paper/all`, `out_paper/a015`, `out_paper/arxiv_maze` and `out_paper/arxiv_corridor`; their printed values are in "
        "`heatmaps/printed_values/`), `map_<scen>_<label>_text.pdf` (`build_maps_report.py`; "
        "`maps/maps_report_manifest_Default_New1_New2_New3_New4_New5.json`), the column and compact maps of the map builder "
        "(`maps/maps_top5_manifest.json`, `maps/compact_manifest.json`), and the default-setting figures of "
        "`build_default_figures.py` (data `per_trajectory/per_trajectory_values.csv`, setting=Default).\n")
    add("`run_heatmaps.sh` writes the paper-style heatmaps to `out_paper/<name>` and the report-style heatmaps to "
        "`out_report/<name>`; the files in `heatmaps/printed_values/` are named after these folders "
        "(`out_<style>__<name>__heatmap_values.json`).\n")

    add("\n## Tables\n")
    add("Every value of the table fragments: `tables/tables_values_flat.csv` (also `tables/tables_values.json.gz`); filter on the "
        "columns `fragment` and `target` (paper or report). The fragment files are `appendices/generated_n200/{paper,report}/<fragment>.tex` "
        "of the paper source tree. Generator for all of them: `python3 tools/n200_tables/build_tables.py`, then `check_tables.py`, "
        "run in the paper source tree (repository README, Rebuilding). Summary data behind each: `summary/*.csv` "
        "(family named below).\n")
    add("| Fragment | Printed as | Data (filter) | Underlying summary rows |")
    add("|---|---|---|---|")
    T = [("main_maze_arxiv_v8", "arXiv Table I (Sec. VI)", "fragment=main_maze_arxiv_v8, target=paper",
          "common_metrics_final.csv: main_grid_200 at S1..S5 (four variants) and matched_baselines_200 at alpha0=0.30, Maze"),
         ("main_corridor_arxiv_v8", "arXiv Table IV", "fragment=main_corridor_arxiv_v8, target=paper", "same, Corridor"),
         ("stress_representative_arxiv_percent", "arXiv Table V", "fragment=stress_representative_arxiv_percent, target=paper",
          "common_metrics_new_anchors.csv: stress_drift_200 at S3 with V-L (phi x L_mu)"),
         ("assumption_checks_arxiv", "numbers in the text of arXiv App. F-D", 
          "fragment=assumption_checks_arxiv (columns name, field, value, scope)", "main_grid_200, detailed_200, budget_normalizer_200 and stress_drift_200 views (counts of configurations, trajectories, checked decisions and unresolved events; column scope)"),
         ("table_key", "key box of report Section 2.1.3 (no caption)", "none (text)", "-"),
         ("design_overview", "report Table 3", "fragment=design_overview, target=report (columns fixed, sets)", "compared values per study (common_metrics_*.csv column family)"),
         ("default_{scen}", "report Tables 4 (Corridor) and 5 (Maze)", "fragment=default_corridor / default_maze, target=report",
          "common_metrics_final.csv: main_grid_200 at the default setting (four variants) and matched_baselines_200 at alpha0=0.15"),
         ("topfive_full_{scen}", "report Tables 6 (Corridor) and 7 (Maze)", "fragment=topfive_full_*, target=report",
          "common_metrics_final.csv: main_grid_200 at S1..S5 and matched_baselines_200 at alpha0=0.30"),
         ("star_maze", "report Table 8", "fragment=star_maze, target=report; tables/starred_setting.json",
          "S* = (0.25, 0.225, 0.15, 50), Maze, four variants: values from the N = 200 table of the paper source (columns source_file, source_line); not on the main grid"),
         ("delta_{scen}", "report Tables 9 (Corridor) and 10 (Maze)", "fragment=delta_*, target=report", "main_grid_200, delta swept at each setting"),
         ("period_{scen}", "report Tables 11 (Corridor) and 12 (Maze)", "fragment=period_*, target=report", "main_grid_200, K swept"),
         ("b_{scen}", "report Tables 13 (Corridor) and 14 (Maze)", "fragment=b_*, target=report", "budget_normalizer_200 (final: default setting; new_anchors: S1..S5)"),
         ("phi_{scen}", "report Tables 15 (Corridor) and 16 (Maze)", "fragment=phi_*, target=report", "stress_drift_200 (phi x L_mu)"),
         ("stress_{scen}_cases", "report Tables 17 (Corridor) and 18 (Maze)", "fragment=stress_*_cases, target=report", "common_metrics_new_anchors.csv: stress_drift_200 at S1..S5 (event counts)"),
         ("checks_{scen}_L0005, _L005", "report Tables 19, 20 (Corridor) and 21, 22 (Maze)", "fragment=checks_*_L0005 / _L005, target=report",
          "detailed_200 and matched_baselines_200 (risk-increase checks, assumption flags)"),
         ("variation_{scen}", "report Tables 23 (Corridor) and 24 (Maze)", "fragment=variation_*, target=report", "temporal_metrics_*.csv (temporal variance, step RMS over the first 200 decisions)")]
    for t in T:
        add("| " + " | ".join(f"`{x}`" if i == 0 else x for i, x in enumerate(t)) + " |")
    add("\nThe other fragments in the same files are written by the table builder but are not printed in the arXiv paper or the report: "
        "`arxiv_design_overview`, `checks_{scen}_arxiv`, `checks_{scen}_L0005_nodef`, `detailed_events_{scen}`, `main_{scen}_extended_v8`, "
        "`n30_record_{scen}`, `sec6_combined`, "
        "`sec6_maze_arxiv`, `sec6_short`, `secondary_{scen}_arxiv`, `stress_{scen}_arxiv`, `stress_{scen}_arxiv_percent`, "
        "`table_key_arxiv`, `topfive_full_{scen}_arxiv_percent`, `topfive_full_{scen}_onepage`, `variation_{scen}_arxiv`, and the "
        "`target=paper` versions of the report fragments.\n")
    add("\n## Data files\n")
    add("| File | Content |")
    add("|---|---|")
    for f, c in (("metrics/final/metrics.json.gz, metrics/new_anchors/metrics.json.gz", "summary metrics per view and cell (the exports' `common/metrics.json`; see README.md, Sources)"),
                 ("tables/tables_values.json.gz, tables/tables_values_flat.csv", "every value of the table fragments with provenance (one row per cell in the CSV)"),
                 ("tables/tables_values_meta.json", "generator, inputs and their SHA-256, rounding and count rules, fragment lists"),
                 ("tables/starred_setting.json", "setting S* (V-C, alpha=0.225): every printed value with source file and line"),
                 ("catalog/catalog_values.json.gz, catalog/CATALOG_README.md", "values of the catalog table fragments per style and the catalog index (PDFs, order, build)"),
                 ("summary/common_metrics_final.csv, common_metrics_new_anchors.csv", "all metric means, SDs, counts, upper bounds, check flags, assessed populations per view and cell (3914 and 720 rows)"),
                 ("summary/temporal_metrics_final.csv, temporal_metrics_new_anchors.csv", "temporal variance and step RMS summaries per view and cell"),
                 ("per_trajectory/per_trajectory_values.csv", "per-trajectory values behind the distribution figures, 6 settings x 2 scenarios x 9 methods x 200 seeds"),
                 ("per_trajectory/boxplot_statistics.csv", "box statistics per panel as drawn in the distribution figures"),
                 ("per_trajectory/source/setting_assessments.json", "the per-trajectory loader cache the two files above come from"),
                 ("sweeps/sweep_series.csv", "series of the sweep plots (means, SDs, counts, bounds), baselines included"),
                 ("sweeps/sweepbox_T3_values.csv", "per-trajectory m_c,k/c_max, q+_CBF,k, T_mis,k and q_dep flags behind `sweepbox_T3_compact.pdf` (S3; V-L, V-C; "
                  "swept alpha, delta, K, b; 200 seeds per box; cell means, q_dep counts and 95% upper bound equal metrics.json)"),
                 ("sweeps/sweepbox_ext_values.csv.gz", "per-trajectory m_c,k/c_max, q+_CBF,k, q+_unsafe,k, T_mis,k and q_cert/q_dep flags behind `sweepbox_<scen>_<label>.pdf` "
                  "(1024 configurations x 200 seeds; gzip-compressed form of the CSV that `build_sweeps_box_ext.py` writes)"),
                 ("heatmaps/heatmap_grid.csv", "full main grid, every metric"),
                 ("heatmaps/printed_values/*.json", "heatmap builder output (printed numbers per figure)"),
                 ("maps/*.json", "map manifests (settings, seeds, per-seed records, payload hashes)"),
                 ("COLUMNS.csv", "column dictionary of every CSV"),
                 ("VERIFICATION.md", "check of 50 CSV values against the original metrics.json")):
        add(f"| `{f}` | {c} |")
    (OUT / "INDEX.md").write_text("\n".join(L) + "\n")
    print("INDEX.md", len(L), "lines")


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    main()
