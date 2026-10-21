#!/usr/bin/env bash
# Catalog heatmaps: scenario x metric (m_c, q_CBF, q_unsafe) x alpha0 {0.15,0.20,0.25,0.30}, paper (505.89 pt, 10 pt) and report (505.89 pt, 11 pt);
# d_all = five delta slices per figure; two-delta (arXiv style) = consecutive delta pairs.
# The jobs run in parallel; the script exits with status 1 if any job fails.
set -e
PY=${PYTHON:-python3}   # interpreter (Python 3.12 with numpy and matplotlib); override with PYTHON=...
cd "$(dirname "$0")"
O=../../figures/n200/catalog   # paper source tree layout: tools/n200_catalog/ next to figures/
A="0.15 0.20 0.25 0.30"; M="m_c q_CBF q_unsafe"; D="0.025 0.05 0.1 0.15 0.2"
REP=(--width-texpt 505.89 --text-pt 11 --rotate-delta --bold-narrow --cell-h 22 --sd-offset 5.3)
pids=()
$PY build_heatmaps.py --out $O/paper/heatmaps/d_all --scenarios corridor maze --alpha0 $A --deltas $D --deltas-per-figure 5 --metrics $M --tag all &
pids+=($!)
$PY build_heatmaps.py --out $O/report/heatmaps/d_all "${REP[@]}" --scenarios corridor maze --alpha0 $A --deltas $D --deltas-per-figure 5 --metrics $M --tag all &
pids+=($!)
for P in "0.025 0.05" "0.05 0.1" "0.1 0.15" "0.15 0.2"; do
  t=${P// /-}; t=${t//./p}
  $PY build_heatmaps.py --out $O/paper/heatmaps/two_delta_$t --scenarios corridor maze --alpha0 $A --deltas $P --deltas-per-figure 2 --metrics $M &
  pids+=($!)
  $PY build_heatmaps.py --out $O/report/heatmaps/two_delta_$t "${REP[@]}" --scenarios corridor maze --alpha0 $A --deltas $P --deltas-per-figure 2 --metrics $M &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
if [ "$fail" -ne 0 ]; then
  echo "HEATMAPS FAILED" >&2
  exit 1
fi
echo HEATMAPS DONE
