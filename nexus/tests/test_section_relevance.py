"""Only the settings the requested setup reads are shown; every card and sub-tab has a name,
an icon and (cards) a line saying what it decides."""

import pytest
from PyQt6.QtWidgets import QApplication

from prestus_nexus.model import ConfigModel
from prestus_nexus.session import Session
from prestus_nexus.ui import icons
from prestus_nexus.ui.main_window import SUBTAB_ICONS, MainWindow
from prestus_nexus.ui.section import SectionPage


def test_water_hides_head_settings(schema):
    m = ConfigModel(schema)
    head = ["path.anat", "path.seg", "path.t1_pattern", "path.t2_pattern", "startup.simnibs_bin_path",
            "io.save_MNI", "io.deface_plots", "analysis.target_mni_mm",
            "modules.run_targeting_analysis", "modules.run_posthoc_water_sims"]
    assert all(m.is_visible(k) for k in head)                         # layered by default
    m.set("simulation.medium", "water")
    assert not any(m.is_visible(k) for k in head)
    m.set("placement.mode", "mni")                                    # MNI placement reads the segmentation
    assert m.is_visible("path.anat") and m.is_visible("startup.simnibs_bin_path")
    assert not m.is_visible("path.t2_pattern")


def test_scheduler_fields_follow_the_computer(schema):
    m = ConfigModel(schema)
    m.set("platform", "auto")
    assert not m.is_visible("hpc.partition")                         # no scheduler known
    m.context["__scheduler"] = True
    assert m.is_visible("hpc.partition") and m.is_visible("hpc.posthoc_water_timelimit")
    m.set("modules.run_posthoc_water_sims", 0)
    assert not m.is_visible("hpc.posthoc_water_timelimit")
    m.context["__scheduler"] = False
    m.set("platform", "slurm")
    assert m.is_visible("hpc.partition")
    m.set("platform", "matlab")
    assert not m.is_visible("hpc.partition")


def test_section_names_icons_subtitles(schema):
    names = {}
    for key, spec in schema.sections.items():
        assert spec.get("label") and spec.get("sub"), key
        assert spec.get("icon") in icons.names(), key
        names.setdefault(spec["label"], []).append(key)
    assert not [n for n, ks in names.items() if len(ks) > 1]          # one name per card
    assert schema.get("jwave").label == schema.section_label("jwave")  # no second name for a group
    assert all(v in icons.names() for v in SUBTAB_ICONS.values())


@pytest.fixture
def window(schema, prestus_root, tmp_path):
    app = QApplication.instance() or QApplication([])  # noqa: F841
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "s.yaml"), prestus_root)
    yield w
    w.dirty = False
    w.close()


def test_cards_and_subtabs_have_icons(window):
    cards = [p for p in window.pages.values() if isinstance(p, SectionPage)]
    assert cards and all(p.icon_name in icons.names() and p.icon.pixmap() is not None for p in cards)
    for i in range(window.tabs.count()):
        sub = getattr(window.tabs.widget(i), "subtabs", None)
        if sub is not None:
            assert all(not sub.tabIcon(j).isNull() for j in range(sub.count())), window.tabs.tabText(i)


def test_derived_keys_hidden_until_set(window):
    io = window.pages["io"]
    shown = lambda: [r.node.path for r in io.rows if io._row_shown(r)]   # noqa: E731
    assert "io.adopted_heatmap" not in shown()
    io.set_show_advanced(True)
    assert "io.adopted_heatmap" not in shown()                        # empty: hidden even with advanced
    window.model.set("io.adopted_heatmap", "/x/end.nii.gz")
    assert "io.adopted_heatmap" in shown()
