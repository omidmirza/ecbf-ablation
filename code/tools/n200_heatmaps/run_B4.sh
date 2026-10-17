#!/usr/bin/env bash
# Heatmaps printed in the paper (out_B4_paper) and the report (out_B4_report); report width 505.89 pt at 11 pt.
set -e
PY=${PYTHON:-python3}   # interpreter (Python 3.12 with numpy and matplotlib); override with PYTHON=...
cd "$(dirname "$0")"
D="0.025 0.05 0.1 0.15 0.2"
$PY build_heatmaps.py --out out_B4_paper/all --scenarios corridor maze --alpha0 0.30 --deltas $D --deltas-per-figure 5 --metrics m_c q_CBF q_unsafe --tag all
$PY build_heatmaps.py --out out_B4_paper/arxiv_maze --scenarios maze --alpha0 0.30 --deltas 0.05 0.1 --deltas-per-figure 2 --metrics m_c q_unsafe
$PY build_heatmaps.py --out out_B4_paper/arxiv_corridor --scenarios corridor --alpha0 0.30 --deltas 0.15 0.2 --deltas-per-figure 2 --metrics q_CBF q_unsafe
$PY build_heatmaps.py --out out_B4_report/a --width-texpt 505.89 --text-pt 11 --rotate-delta --bold-narrow --cell-h 22 --sd-offset 5.3 --scenarios corridor maze --alpha0 0.15 0.30 --deltas $D --deltas-per-figure 5 --metrics m_c q_CBF --tag all
$PY build_heatmaps.py --out out_B4_report/b --width-texpt 505.89 --text-pt 11 --rotate-delta --bold-narrow --cell-h 22 --sd-offset 5.3 --scenarios corridor maze --alpha0 0.30 --deltas $D --deltas-per-figure 5 --metrics q_unsafe --tag all
$PY build_heatmaps.py --out out_B4_report/c --width-texpt 505.89 --text-pt 11 --rotate-delta --bold-narrow --cell-h 22 --sd-offset 5.3 --scenarios corridor --alpha0 0.25 --deltas $D --deltas-per-figure 5 --metrics q_unsafe --tag all
echo HEATMAPS DONE
