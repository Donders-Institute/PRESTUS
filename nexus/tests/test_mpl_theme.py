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
