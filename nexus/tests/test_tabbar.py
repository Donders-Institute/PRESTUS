"""Every tab row in the Nexus is an EmphasisTabBar: light dividers between neighbouring tabs,
and a tab stretched to share the bar keeps its icon beside its label."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtWidgets import QApplication, QTabBar, QTabWidget  # noqa: E402

from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui import icons  # noqa: E402
from prestus_nexus.ui.tabbar import EmphasisTabBar, tab_widget  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _all_tab_bars(app, window) -> dict:
    """Each tab and sub-tab in turn, so pages built on first show exist."""
    seen = {}
    for i in range(window.tabs.count()):
        window.tabs.setCurrentIndex(i)
        app.processEvents()
        for sub in window.tabs.widget(i).findChildren(QTabWidget):
            for j in range(sub.count()):
                sub.setCurrentIndex(j)
                app.processEvents()
        for bar in window.findChildren(QTabBar):
            seen[id(bar)] = bar
    return seen


def test_every_tab_bar_is_the_nexus_tab_bar(app, schema, prestus_root, tmp_path):
    from prestus_nexus.ui.main_window import MainWindow
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "s.yaml"), prestus_root)
    w.resize(1000, 800)
    w.show()
    bars = _all_tab_bars(app, w)
    assert len(bars) >= 5, "the walk over tabs found too few tab bars to mean anything"
    plain = [type(b.parentWidget()).__name__ for b in bars.values() if not isinstance(b, EmphasisTabBar)]
    assert not plain, f"tab bars without the Nexus dividers: {plain}"
    w.dirty = False
    w.close()


@pytest.mark.parametrize("with_icons", [True, False])
def test_stretched_tab_keeps_icon_beside_label(app, with_icons):
    tabs = tab_widget()                                       # expanding, like the sub-tabs
    for name in ("Settings", "Calibration", "Equipment"):
        page = QTabWidget()
        i = tabs.addTab(page, name)
        if with_icons:
            tabs.setTabIcon(i, icons.icon("advanced"))
    tabs.resize(1200, 300)
    tabs.show()
    app.processEvents()
    bar = tabs.tabBar()
    for i in range(bar.count()):
        tab, label = bar.tabRect(i), bar.label_rect(i)
        assert tab.width() > 300                              # stretched well beyond its contents
        if with_icons:
            # the icon + label group is narrow and centred, not spread over the whole tab
            assert label.width() < tab.width() / 2
            assert abs(label.center().x() - tab.center().x()) <= 1
        else:
            assert label == tab                               # negative control: text-only tabs as before
    xs = bar.separators()
    assert len(xs) == bar.count() - 1
    assert xs == [bar.tabRect(i).right() + 1 for i in range(bar.count() - 1)]
    bar.setTabVisible(1, False)
    assert len(bar.separators()) == bar.count() - 2           # hidden tabs get no divider
    tabs.close()
