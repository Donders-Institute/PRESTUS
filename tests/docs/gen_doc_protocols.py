#!/usr/bin/env python3
"""Write the example-protocol figures of doc_multitransducer_async.md from the Nexus drawings.

Each example protocol is built as a Nexus config (the same edits the Sonication sequence
card makes) and drawn by the Nexus itself: the card's step strip (SequenceStrip) or the
timing diagram with its firing-schedule lanes (TimingDiagram). Both are painted widgets,
so they render to text SVG (QSvgGenerator) without raster images. Fonts are set to the
docs font stack. Needs the Nexus environment (PyQt6 with QtSvg).

Usage (from the PRESTUS root):
    python tests/docs/gen_doc_protocols.py      # (re)write doc/assets/figures/fig_protocol-*.svg
Re-run after changing the sequence card, the timing diagram, or the examples below.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "nexus"))
OUT = ROOT / "doc" / "assets" / "figures"
FONT = "Roboto, 'Helvetica Neue', Arial, 'DejaVu Sans', sans-serif"   # stylesheets / figures
WIDTH = 700          # timing diagram: drawn narrow so its labels stay legible
STRIP_WIDTH = 760    # step strip: wide enough for the step labels

from PyQt6.QtCore import QRect  # noqa: E402
from PyQt6.QtGui import QPainter  # noqa: E402
from PyQt6.QtSvg import QSvgGenerator  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

app = QApplication.instance() or QApplication([])

from prestus_nexus import sequence as sq  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.schema import Schema  # noqa: E402
from prestus_nexus.ui.export import PaintedExporter  # noqa: E402
from prestus_nexus.ui.sequence_card import SequenceCard  # noqa: E402
from prestus_nexus.ui.timing_view import TimingPanel  # noqa: E402

LEFT = ([-52, -8, 58], [-29, -6, 7])          # tutorial_ernie.yaml stance [MNI mm]
RIGHT = ([52, -8, 58], [29, -6, 7])           # its mirror image


def base(schema: Schema) -> ConfigModel:
    """20 ms pulses at 5 Hz for 20 s, 60 s steady state (tutorial_ernie.yaml)."""
    m = ConfigModel(schema)
    for k, v in {"timing.pd": 0.02, "timing.pri": 0.2, "timing.ptd": 20.0, "timing.ptri": 20.0,
                 "timing.ptrd": 20.0, "timing.pt_timestep": 0.02, "timing.post_pt_timestep": 1.0,
                 "timing.post_ptri_dur": 60.0, "modules.run_heating_sims": 1}.items():
        m.set(k, v)
    m.set("transducer.serial", "CTX_500_026")
    m.set("transducer.trans_pos", LEFT[0])
    m.set("transducer.focus_pos", LEFT[1])
    return m


def place(m: ConfigModel, ti: int, stance) -> None:
    m.set("transducer.trans_pos", stance[0], ti)
    m.set("transducer.focus_pos", stance[1], ti)


def sequence(schema):
    """Left, right, left again (twice as long), 30 s pauses."""
    m = base(schema)
    sq.add_step(m, 30.0)
    place(m, 1, RIGHT)
    sq.add_step(m, 30.0)                      # T3 keeps T1's stance: the field is reused
    sq.set_timing(m, 2, {"ptrd": 40.0})
    return m


def together(schema):
    """Left and right fire together, then left alone."""
    m = sequence(schema)
    sq.join_previous(m, 1)
    return m


def mixed(schema):
    """Left alone; two transducers together on a deep target; left again."""
    m = base(schema)
    sq.add_step(m, 30.0)
    place(m, 1, ([-72, -22, 8], [-12, -18, 6]))
    sq.add_step(m, 0.0)
    place(m, 2, ([-50, -60, 50], [-12, -18, 6]))
    sq.join_previous(m, 2)
    sq.add_step(m, 30.0)                      # T4 keeps T1's stance
    sq.set_timing(m, 2, {"ptrd": 40.0})
    m.set("simulation.uncertainty", True)
    return m


def interleaved(schema):
    """Two transducers in the same 1 s trains, T2's pulses 100 ms after T1's."""
    m = base(schema)
    for k, v in {"timing.ptd": 1.0, "timing.ptri": 2.0}.items():
        m.set(k, v)
    m.add_transducer()
    m.set("simulation.transducer_coupling", "async")
    m.set("timing.async_schedule", "custom")
    place(m, 1, RIGHT)
    m.set("transducer.timing_offset_s", 0.0, 0)
    m.set("transducer.timing_offset_s", 0.1, 1)
    return m


def strip_svg(m: ConfigModel, path: Path) -> None:
    card = SequenceCard(m)
    card.refresh()
    strip = card.strip
    strip.setFixedSize(STRIP_WIDTH, strip.H)
    gen = QSvgGenerator()
    gen.setFileName(str(path))
    gen.setSize(strip.size())
    gen.setViewBox(QRect(0, 0, strip.width(), strip.height()))
    gen.setTitle(path.stem)
    qp = QPainter(gen)
    strip.render(qp)
    qp.end()


def timing_svg(m: ConfigModel, path: Path) -> None:
    panel = TimingPanel(m)
    panel.diagram.resize(WIDTH, panel.diagram.sizeHint().height())
    PaintedExporter(panel.diagram).save(path, "svg", None)


def tidy(path: Path) -> None:
    svg = path.read_text(encoding="utf-8")
    if "<image" in svg:
        raise SystemExit(f"{path.name}: contains a raster image; draw it as vector instead")
    svg = re.sub(r'font-family="[^"]*"', f'font-family="{FONT}"', svg)
    # QSvgGenerator gives the size in mm at 90 dpi: use the drawn pixels, as on screen
    box = re.search(r'viewBox="0 0 (\d+) (\d+)"', svg)
    if box:
        svg = re.sub(r'<svg width="[^"]*" height="[^"]*"',
                     f'<svg width="{box.group(1)}" height="{box.group(2)}"', svg, count=1)
    path.write_text(svg, encoding="utf-8")


FIGURES = {
    "fig_protocol-sequence_diagram.svg": (sequence, strip_svg),
    "fig_protocol-together_diagram.svg": (together, strip_svg),
    "fig_protocol-mixed_diagram.svg": (mixed, strip_svg),
    "fig_protocol-interleaved_diagram.svg": (interleaved, timing_svg),
}


def main() -> None:
    schema = Schema(ROOT / "config" / "config_default.yaml")
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (build, draw) in FIGURES.items():
        m = build(schema)
        draw(m, OUT / name)
        tidy(OUT / name)
        print(f"wrote {name}: {sq.describe(m) or 'no sequence'}")


if __name__ == "__main__":
    main()
