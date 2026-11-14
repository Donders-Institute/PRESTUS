"""ITRUSST intensity notation (I<sub>SPPA</sub>) in rich text, path fields that show their whole
value, and hover text naming the transducer under the cursor (slices and 3D)."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import QPoint, Qt  # noqa: E402
from PyQt6.QtGui import QHelpEvent  # noqa: E402
from PyQt6.QtCore import QEvent  # noqa: E402
from PyQt6.QtWidgets import QApplication, QLabel, QToolTip  # noqa: E402

from prestus_nexus import labels  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def window(app, schema, prestus_root, tmp_path_factory):
    from prestus_nexus.ui.main_window import MainWindow
    w = MainWindow(ConfigModel(schema), Session(tmp_path_factory.mktemp("s") / "session.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def test_intensity_notation_helpers():
    assert labels.rich("Target free-water ISPPA") == "Target free-water I<sub>SPPA</sub>"
    assert labels.rich("ISPTA, ISPPA.") == "I<sub>SPTA</sub>, I<sub>SPPA</sub>."
    assert labels.rich("multi-ISPPA sweep") == "multi-ISPPA sweep"            # a pipeline name
    assert labels.rich("desc-multiIsppa_report.html") == "desc-multiIsppa_report.html"
    assert labels.mathtext("Target ISPPA [W/cm²]") == r"Target I$_{\mathrm{SPPA}}$ [W/cm²]"


def test_plain_label_stays_plain(schema):
    node = schema.get("transducer.target_isppa_wcm2")
    assert labels.label(node) == "Target free-water ISPPA"        # combos, Markdown, file names
    assert labels.rich_label(node) == "Target free-water I<sub>SPPA</sub>"


def test_field_tooltip_uses_subscripts(schema):
    from prestus_nexus.ui.tooltip import html_for
    node = schema.get("transducer.target_isppa_wcm2")
    tip = html_for(node, labels.label(node), "30")
    assert "I<sub>SPPA</sub>" in tip and "&lt;sub" not in tip and "ISPPA</div>" not in tip


def test_form_labels_render_subscripts_not_raw_html(window):
    from prestus_nexus.ui.section import SectionPage
    rows = [r for page in window.findChildren(SectionPage) for r in page.rows]
    row = next(r for r in rows if r.node.path == "transducer.target_isppa_wcm2")
    assert "I<sub>SPPA</sub>" in row.label.text()
    # every label holding the subscript markup renders it as rich text (never shown raw)
    for lbl in window.findChildren(QLabel):
        if "<sub>" in lbl.text():
            assert lbl.textFormat() != Qt.TextFormat.PlainText, lbl.text()
            assert lbl.textFormat() == Qt.TextFormat.RichText or Qt.mightBeRichText(lbl.text())


def test_path_field_shows_start_and_full_value_on_hover(app, schema):
    from prestus_nexus.ui.fields import TextEditor, full_value_tip
    value = "/opt/simnibs/lib:/usr/lib64:/opt/matlab/R2023b/runtime/glnxa64:/opt/cuda/7.2.0/lib64"
    assert full_value_tip(value).split("<br>") == value.split(":")
    assert full_value_tip("/one/folder") == "/one/folder"
    ed = TextEditor(schema.get("hpc.ld_library_path"))
    ed.resize(160, 30)
    ed.show()
    ed.set_value(value, "")
    assert ed.edit.cursorPosition() == 0                 # the start, not the clipped tail
    ev = QHelpEvent(QEvent.Type.ToolTip, QPoint(5, 5), ed.edit.mapToGlobal(QPoint(5, 5)))
    assert ed.eventFilter(ed.edit, ev) is True
    assert QToolTip.isVisible() and "/opt/cuda/7.2.0/lib64" in QToolTip.text()
    QToolTip.hideText()
    ed.edit.setToolTip("Cannot parse")                   # an error keeps its own tooltip
    assert ed.eventFilter(ed.edit, ev) is False
    ed.close()


def test_marker_tips_pick_radius(app):
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
    from matplotlib.figure import Figure
    from prestus_nexus.ui.tooltip import MarkerTips
    fig = Figure(figsize=(3, 3))
    canvas = FigureCanvasQTAgg(fig)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    tips = MarkerTips(canvas)
    tips.add(ax, 50, 50, "Transducer 1: CTX-500")
    tips.add(ax, np.nan, 10, "never")                    # off-slice/NaN points are skipped
    x, y = ax.transData.transform((50, 50))
    assert tips.at(ax, x + 3, y - 3) == "Transducer 1: CTX-500"
    assert tips.at(ax, x + 40, y) == ""
    tips.clear(ax)
    assert tips.at(ax, x, y) == ""


def test_ortho_view_keeps_marker_texts(app):
    pytest.importorskip("matplotlib")
    from prestus_nexus.ui.ortho_view import OrthoViewer
    v = OrthoViewer()
    v.set_markers([(np.zeros(3), "s", "#f00", "Transducer 1: X"), (np.ones(3), "x", "#f00")])
    assert [m[1] for m in v.markers] == ["s", "x"]
    assert v.marker_texts == ["Transducer 1: X", ""]


def test_3d_scene_picks_transducers_only(app):
    from PyQt6.QtCore import QByteArray, QUrl
    from PyQt6.QtQml import QQmlComponent, QQmlEngine
    from PyQt6.QtWidgets import QWidget
    view3d = pytest.importorskip("prestus_nexus.ui.view3d")
    if not view3d.HAVE_Q3D:
        pytest.skip("no Qt Quick 3D")
    m = view3d._tipped((np.zeros((3, 3)), np.zeros((1, 3), int), "#fff", 1.0), "Transducer 1: X")
    v, f, c, a = m                                       # still unpacks as a 4-tuple
    assert m.tip == "Transducer 1: X" and a == 1.0
    engine = QQmlEngine()
    host = QWidget()
    engine.rootContext().setContextProperty("tipSink", view3d._TipSink(host))
    engine.rootContext().setContextProperty("meshes", [])
    comp = QQmlComponent(engine)
    comp.setData(QByteArray(view3d.QML.encode()), QUrl("file:///placement_3d.qml"))
    assert comp.errors() == [] and comp.create() is not None
