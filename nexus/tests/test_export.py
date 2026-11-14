"""Figure export (ui/export.py): the download button on matplotlib and hand-painted figures,
light-theme export with the screen theme restored, and dark-mode painting."""

import xml.etree.ElementTree as ET

import pytest

pytest.importorskip("matplotlib")
from matplotlib.colors import to_hex  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from PyQt6.QtGui import QImage  # noqa: E402

from prestus_nexus.ui import export, mpl_theme, style  # noqa: E402


@pytest.fixture
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def dark():
    old = style._dark
    style._dark = True
    yield
    style._dark = old


@pytest.fixture(autouse=True)
def _session():
    keep = dict(export._SESSION)
    yield
    export._SESSION.update(keep)


def _canvas(qapp, **kw):
    fig = Figure(figsize=(3, 2))
    ax = fig.add_subplot()
    ax.plot([0, 1], [0, 1], color=style.SERIES[0], label="data")
    ax.set_xlabel("depth [mm]")
    ax.legend()
    canvas = mpl_theme.ThemedCanvas(fig, export_name="test_view", **kw)
    canvas.resize(300, 200)
    canvas.draw()
    return canvas, ax


def test_canvas_has_button_pinned_top_right(qapp):
    canvas, _ax = _canvas(qapp)
    btn = canvas.export_button
    assert isinstance(btn, export.ExportButton)
    assert btn.toolTip() == "Export figure"
    canvas.show()
    qapp.processEvents()
    canvas.resize(400, 250)
    qapp.processEvents()
    assert btn.geometry().right() <= canvas.width() and btn.x() > canvas.width() - 40
    btn._fill_menu()
    texts = [a.text() for a in btn.menu().actions() if a.text()]
    assert texts[:5] == ["SVG (vector)…", "PDF (vector)…", "PNG, 300 dpi…", "PNG, 600 dpi…", "Copy image"]
    assert "Keep current theme" in texts
    canvas.close()
    assert mpl_theme.ThemedCanvas(Figure(), export=False).export_button is None


def test_mpl_export_svg_pdf_png(qapp, tmp_path):
    canvas, _ax = _canvas(qapp)
    btn = canvas.export_button
    svg = btn.export(tmp_path / "f.svg")
    ET.parse(svg)                                                    # well-formed SVG
    pdf = btn.export(tmp_path / "f.pdf")
    assert pdf.read_bytes().startswith(b"%PDF")
    p300 = btn.export(tmp_path / "a.png", "png", 300)
    p600 = btn.export(tmp_path / "b.png", "png", 600)
    i300, i600 = QImage(str(p300)), QImage(str(p600))
    assert not i300.isNull() and abs(i600.width() / i300.width() - 2) < 0.05
    assert export._SESSION["dir"] == str(tmp_path)                    # folder remembered


def test_export_light_then_restores_dark(qapp, tmp_path, dark):
    canvas, ax = _canvas(qapp)
    assert to_hex(ax.patch.get_facecolor()) == style.C_DARK["surface"]
    seen = {}
    real = canvas.figure.savefig

    def spy(*a, **k):
        seen["bg"] = to_hex(ax.patch.get_facecolor())
        seen["label"] = to_hex(ax.xaxis.label.get_color())
        return real(*a, **k)
    canvas.figure.savefig = spy
    canvas.export_button.export(tmp_path / "f.png", "png", 300)
    assert seen == {"bg": style.C["surface"], "label": style.C["ink2"]}
    assert to_hex(ax.patch.get_facecolor()) == style.C_DARK["surface"]       # theme restored
    assert to_hex(ax.xaxis.label.get_color()) == style.C_DARK["ink2"]
    export.set_keep_theme(True)
    canvas.export_button.export(tmp_path / "g.png", "png", 300)
    assert seen["bg"] == style.C_DARK["surface"]


def test_palette_colours_follow_theme(qapp):
    """A colour a plot took from the active palette follows the theme by its role."""
    fig = Figure()
    ax = fig.add_subplot()
    ln, = ax.plot([0, 1], [0, 1], color=style.C_DARK["ink3"])
    ax.text(0, 0, "note", color=style.C_DARK["amber"])
    mpl_theme.apply(fig, style.C)
    assert to_hex(ln.get_color()) == style.C["ink3"]
    assert to_hex(ax.texts[0].get_color()) == style.C["amber"]


def test_filename_hint():
    class M:
        def get(self, key):
            return {"subject_id": 7, "io.output_affix": "lowdose"}.get(key)
    assert export.subject_stem(M(), "beam")() == "sub-007_beam_desc-lowdose"
    assert export.subject_stem(lambda: M(), "beam")() == "sub-007_beam_desc-lowdose"


class _Model:
    n_transducers = 1

    def __init__(self):
        self.v = {"transducer.type": "annular", "transducer.annular.elem_od_mm": [20.0, 40.0],
                  "transducer.annular.elem_id_mm": [0.0, 22.0], "transducer.annular.curv_radius_mm": 60.0,
                  "transducer.annular.dist_geom_ep_mm": 55.0, "transducer.focal_distance_ep": 40.0,
                  "subject_id": 1}

    def get(self, key, ti=None):
        return self.v.get(key)

    def device(self, ti):
        return None


def test_transducer_painted_export(qapp, tmp_path):
    from prestus_nexus.ui.transducer_view import TransducerSchematic
    w = TransducerSchematic(_Model())
    w.resize(400, 240)
    btn = w.export_button
    assert btn.default_name("svg") == "sub-001_transducer.svg"
    ET.parse(btn.export(tmp_path / "t.svg"))
    assert btn.export(tmp_path / "t.pdf").read_bytes().startswith(b"%PDF")
    small = QImage(str(btn.export(tmp_path / "t3.png", "png", 3)))
    big = QImage(str(btn.export(tmp_path / "t6.png", "png", 6)))
    assert (small.width(), small.height()) == (1200, 720)
    assert (big.width(), big.height()) == (2400, 1440)
    assert not btn.isHidden()


def _bg(img: QImage) -> str:
    return img.pixelColor(2, img.height() - 2).name()


def test_transducer_dark_on_screen_light_in_export(qapp, dark):
    from prestus_nexus.ui.transducer_view import TransducerSchematic
    w = TransducerSchematic(_Model())
    w.resize(400, 240)
    screen = w.grab().toImage()
    assert _bg(screen) == style.C_DARK["surface"]
    assert _bg(export.PaintedExporter(w).raster(1.0)) == style.C["surface"]
    export.set_keep_theme(True)
    assert _bg(export.PaintedExporter(w).raster(1.0)) == style.C_DARK["surface"]


def test_timing_and_r2_strip_exportable(qapp, tmp_path):
    from prestus_nexus.ui.calibration_plots import R2Strip
    from prestus_nexus.ui.timing_view import TimingDiagram
    strip = R2Strip()
    strip.resize(500, 70)
    strip.set_data([30.0, 40.0], [0.95, 0.7])
    ET.parse(strip.export_button.export(tmp_path / "r2.svg"))
    t = TimingDiagram()
    t.resize(500, 200)
    assert QImage(str(t.export_button.export(tmp_path / "t.png", "png", 3))).width() == 1500


def test_coupling_empty_state(qapp):
    from prestus_nexus.ui.coupling_view import CouplingPrediction

    class M:
        n_transducers = 2

        def effective_transducers(self):
            return [{"target_isppa_wcm2": 10, "freq_hz": 5e5}, {"target_isppa_wcm2": 10, "freq_hz": 5e5}]

        def get(self, key, ti=None):
            return None
    w = CouplingPrediction(M(), lambda: None)
    hits = []
    w.open_placement_requested.connect(lambda: hits.append(1))
    w.toggle.setChecked(True)
    assert w.stack.is_empty()
    assert w.stack.empty.button.text() == "Open placement viewer"
    w.stack.empty.button.click()
    assert hits == [1]
