#!/usr/bin/env python3
"""Report trajectory maps of kind 'pair' (map_<scenario>_<label>_pair.pdf) at the 11 pt report body.
Copy of build_maps_report.py with kind 'appendix' -> kind 'pair' (the paper's pair maps are set at 10 pt, below the
report's 11 pt body).

As in build_maps_report.py (report-width trajectory maps):

Runs the trajectory map builder build_maps_top5.py unchanged (loaded from its file, path ECBF_MAPS_BUILDER; the
builder and the raw trajectory payloads are not part of this repository) with four settings changed in this process
only:
  - text width figs.FULL = 505.89 pt (report \\textwidth: geometry letterpaper, margin 0.75in; 11 pt);
  - body text size figs.TXT = 11 and a floor of 11 on every font size (FontProperties.set_size), so no glyph is
    smaller than the 11 pt report body (the generators also hard-code 10 in places);
  - the figure helper's registered template widths include the report width and its minimum text size is 10.95 TeX pt;
  - output folder = --out (default ./out_report), preview scratch = MAPS_PREVIEW (never the builder's folder).
Only the kind 'pair' (320 pt tall, two maps per page) is drawn here.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_maps/<this file>)
MAPS_BUILDER = Path(os.environ.get("ECBF_MAPS_BUILDER", REPO / "maps_top5" / "build_maps_top5.py"))
FIGURE_STYLE = Path(os.environ.get("ECBF_FIGURE_STYLE", REPO / "figure_style"))   # ieee_figure.py (not included)
sys.path.insert(0, str(FIGURE_STYLE))
REPORT_TEXTWIDTH_TEXPT = 505.89   # report geometry letterpaper, margin 0.75in
REPORT_TXT = 11


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


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", nargs="+", default=["Default", "New1", "New2", "New3", "New4", "New5"])
    ap.add_argument("--out", type=Path, default=HERE / "out_report_pair")
    ap.add_argument("--zoom-extra", type=float, default=2.0, help="pt added to the scaled zoom band")
    ap.add_argument("--preview", type=Path, default=Path(tempfile.gettempdir()) / "ecbf_maps_top5_preview")   # payload cache
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    os.environ["MAPS_PREVIEW"] = str(a.preview)
    sys.argv = [str(MAPS_BUILDER)]
    spec = importlib.util.spec_from_file_location("maps_top5_v9", MAPS_BUILDER)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["maps_top5_v9"] = mod
    spec.loader.exec_module(mod)

    import matplotlib.font_manager as fm
    import ieee_figure as fs
    width_in = REPORT_TEXTWIDTH_TEXPT / 72.0
    figs = mod.base.figs
    figs.FULL = width_in
    figs.TXT = REPORT_TXT
    if not any(abs(width_in - w) < 1e-6 for w in fs.TEMPLATE_WIDTHS_IN):
        fs.TEMPLATE_WIDTHS_IN = tuple(fs.TEMPLATE_WIDTHS_IN) + (width_in,)
    fs.MIN_TEXT_PT = 10.95
    if hasattr(figs, "MIN_TEXT_PT"):
        figs.MIN_TEXT_PT = 10.95
    orig_set_size = fm.FontProperties.set_size

    def set_size(self, size):
        orig_set_size(self, size)
        if self.get_size_in_points() < REPORT_TXT:
            orig_set_size(self, REPORT_TXT)
    fm.FontProperties.set_size = set_size
    import matplotlib.pyplot as plt
    for key in ("font.size", "xtick.labelsize", "ytick.labelsize", "legend.fontsize", "axes.labelsize",
                "axes.titlesize", "legend.title_fontsize"):
        v = plt.rcParams[key]
        if isinstance(v, (int, float)) and v < REPORT_TXT:
            plt.rcParams[key] = REPORT_TXT

    # the 11 pt legend of panel (a) is taller: raise the corridor/planar zoom band by the text ratio (+ 2 pt)
    for k in list(mod.base.ex.ZOOM_H):
        mod.base.ex.ZOOM_H[k] = round(mod.base.ex.ZOOM_H[k] * REPORT_TXT / 10.0 + a.zoom_extra, 1)
    mod.HERE = a.out
    mod.load_all(a.labels)
    info = mod.draw(a.labels, ["corridor", "maze"], kinds=("pair",))
    manifest = dict(generator=recorded_path(Path(__file__)), generator_sha256=sha(Path(__file__).resolve()),
                    v9_builder=recorded_path(MAPS_BUILDER), v9_builder_sha256=sha(MAPS_BUILDER), report_textwidth_texpt=REPORT_TEXTWIDTH_TEXPT,
                    text_pt=REPORT_TXT, zoom_h_pt=dict(mod.base.ex.ZOOM_H), figures=info)
    (a.out / f"maps_report_pair_manifest_{'_'.join(a.labels)}.json").write_text(json.dumps(manifest, indent=1, default=str) + "\n")
    print("DONE", len(info))


if __name__ == "__main__":
    main()
