"""Figures follow the Nexus theme: default backgrounds and ink recoloured, chosen colours kept."""

import pytest

pytest.importorskip("matplotlib")
from matplotlib.colors import to_hex  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from prestus_nexus.ui import mpl_theme, style  # noqa: E402


def _fig():
    fig = Figure()
    ax = fig.add_subplot()
    ax.plot([0, 1], [0, 1], color="#c62828", label="data")
    ax.set_xlabel("depth")
    ax.set_title("warning", color="#ef6c00")
    ax.legend()
    return fig, ax


def test_dark_and_back_to_light():
    fig, ax = _fig()
    mpl_theme.apply(fig, style.C_DARK)
    d = style.C_DARK
    assert to_hex(fig.patch.get_facecolor()) == d["surface"]
    assert to_hex(ax.patch.get_facecolor()) == d["surface"]
    assert to_hex(ax.xaxis.label.get_color()) == d["ink2"]
    assert to_hex(ax.spines["left"].get_edgecolor()) == d["ink3"]
    assert to_hex(ax.title.get_color()) == "#ef6c00"                 # chosen colour kept
    assert to_hex(ax.lines[0].get_color()) == "#c62828"
    assert to_hex(ax.get_legend().get_frame().get_facecolor()) == d["surface"]
    fig.canvas.draw()                                                # ticks made while drawing
    assert to_hex(ax.get_xticklabels()[0].get_color()) == d["ink2"]
    mpl_theme.apply(fig, style.C)                                    # switching back recolours
    assert to_hex(fig.patch.get_facecolor()) == style.C["surface"]
    assert to_hex(ax.xaxis.label.get_color()) == style.C["ink2"]


def test_dark_image_background_kept():
    fig = Figure()
    ax = fig.add_subplot()
    ax.set_facecolor("black")                                        # e.g. an MRI slice view
    mpl_theme.apply(fig, style.C_DARK)
    assert to_hex(ax.patch.get_facecolor()) == "#000000"


def test_canvas_passes_wheel_to_page_unless_it_scrolls():
    """A plot without a scroll_event handler leaves the wheel to the scroll area around it."""
    from PyQt6.QtCore import QPoint, QPointF, Qt
    from PyQt6.QtGui import QWheelEvent
    from PyQt6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])

    def wheel():
        return QWheelEvent(QPointF(10, 10), QPointF(10, 10), QPoint(0, 0), QPoint(0, -120),
                           Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                           Qt.ScrollPhase.NoScrollPhase, False)

    canvas = mpl_theme.ThemedCanvas(_fig()[0], export=False)
    ev = wheel()
    canvas.wheelEvent(ev)
    assert not ev.isAccepted()                       # propagates: the page scrolls
    canvas.mpl_connect("scroll_event", lambda e: None)
    ev = wheel()
    canvas.wheelEvent(ev)
    assert ev.isAccepted()                           # slice / zoom views keep the wheel


@pytest.mark.parametrize("facecolor", [None, "surface"])     # matplotlib's white, or the card colour set
@pytest.mark.parametrize("in_card", [False, True])
def test_figure_background_follows_its_container_and_exports_white(facecolor, in_card):
    """On the page a figure takes the page colour (no white block); in a card it stays card
    white; a download is white either way."""
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication, QFrame, QVBoxLayout
    app = QApplication.instance() or QApplication([])  # noqa: F841
    host = QFrame()
    if in_card:
        host.setObjectName("Card")
    fig = Figure(facecolor=style.active()["surface"]) if facecolor else Figure()
    fig.add_subplot(111).plot([0, 1], [0, 1])
    canvas = mpl_theme.ThemedCanvas(fig)
    QVBoxLayout(host).addWidget(canvas)
    canvas.draw()
    expected = style.active()["surface"] if in_card else style.active()["bg"]
    assert to_hex(fig.patch.get_facecolor()) == expected.lower()
    if canvas.export_button is not None:
        from prestus_nexus.ui.export import MplExporter
        img = MplExporter(canvas).image()
        assert img.pixelColor(1, 1).name() == "#ffffff"           # downloads: white
