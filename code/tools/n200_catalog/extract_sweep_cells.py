#!/usr/bin/env python3
"""Compact per-cell extract of the N200 exports for the sweep plots. Read only.

Inputs (folder ECBF_N200_NUMBERS, not part of this repository): N200_numbers/final/common/metrics.json (main_grid_200, stress_drift_200,
budget_normalizer_200, detailed_200, matched_baselines_200, matched_ACI_200) and .../new_anchors/common/metrics.json
(stress_drift_200, budget_normalizer_200, detailed_200 around the top-five settings T1-T5).
Writes sweep_cells.json: one record per (source file, view, cell) with the family, the anchor of the view, the
configuration and, per metric, mean / sd / applicable, plus the qcert and qdep trajectory-event counts.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/<dir>/<this file>)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports
SRC = {"final": NUM / "final/common/metrics.json",
       "new_anchors": NUM / "new_anchors/common/metrics.json"}
OUT = Path(__file__).resolve().parent / "sweep_cells.json"
KEYS = ("T_mis", "m_c", "q_CBF", "q_unsafe", "q_use", "p_max", "h_min")


def recorded_path(p):
    """A path as recorded in manifests and caches: '<default folder>/<rest>' under an input folder named by an ECBF_*
    variable, else relative to the repository root, else the file name (no absolute local paths)."""
    p = Path(p).resolve()
    for var, name in (("ECBF_N200_NUMBERS", "N200_numbers"), ("ECBF_N200_DATA", "N200_data"),
                      ("ECBF_FIGURE_STYLE", "figure_style"), ("ECBF_DISPLAYS", "data")):
        if os.environ.get(var):
            try:
                return f"{name}/{p.relative_to(Path(os.environ[var]).resolve()).as_posix()}"
            except ValueError:
                pass
    try:
        return p.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return p.name


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    rows = []
    for name, path in SRC.items():
        d = json.loads(path.read_text())
        for v in d["views"]:
            df = v["definition"]
            for c in v["cells"]:
                st = c["statistics"]
                ev = c.get("trajectory_events") or {}
                rows.append(dict(
                    source=name, view_id=v.get("view_id"), family=df["family"], anchor=df.get("anchor"),
                    anchor_setting=df.get("anchor_setting"), configuration=c["configuration"], n=c["n"],
                    stats={k: dict(mean=st[k].get("available_mean", st[k].get("mean")), sd=st[k].get("sd"),
                                   applicable=st[k].get("applicable")) for k in KEYS if k in st},
                    events={f: {x: ev[f].get(x) for x in ("confirmed", "possible", "applicable_n", "upper95_possible")}
                            for f in ("qcert", "qdep") if ev.get(f)}))
        print(name, len(rows), flush=True)
    OUT.write_text(json.dumps(dict(generator=recorded_path(Path(__file__)),
                                   inputs={recorded_path(p): sha_file(p) for p in SRC.values()}, cells=rows)) + "\n")
    print("->", OUT, len(rows))


if __name__ == "__main__":
    main()
