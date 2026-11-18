"""Import trajectory of a Localite XML: which config values choose the series (Qt-free part)."""

from prestus_nexus import localite as lc
from prestus_nexus.model import ConfigModel

SERIES = [
    {"index": 1, "description": "", "is_pose": False, "tracked": True},      # planned point
    {"index": 2, "description": "L_M1", "is_pose": True, "tracked": True, "n_markers": 80, "pos_dev_mm": 1.25},
    {"index": 3, "description": "", "is_pose": True, "tracked": False},     # at the origin
    {"index": 4, "description": "dup", "is_pose": True, "tracked": True, "n_markers": 1},
    {"index": 5, "description": "dup", "is_pose": True, "tracked": True, "n_markers": 80},
]


def test_usable_series_drops_points_and_untracked():
    assert [s["index"] for s in lc.usable_series(SERIES)] == [2, 4, 5]


def test_series_label():
    assert lc.series_label(SERIES[1]) == "Series 2: L_M1 (80 markers, mean deviation 1.2 mm)"
    assert lc.series_label(SERIES[3]) == "Series 4: dup (1 marker)"


def test_unique_description_chooses_by_target(schema):
    m = ConfigModel(schema)
    m.set("placement.localite.position", 4)
    ch = dict(lc.file_import_settings(m, "/x.xml", SERIES[1], SERIES))
    assert ch["placement.mode"] == "localite" and ch["placement.localite.file"] == "/x.xml"
    assert ch["placement.localite.target"] == "L_M1"
    assert ch["placement.localite.position"] is None                   # reset: target decides


def test_shared_description_chooses_by_number(schema):
    m = ConfigModel(schema)
    m.set("placement.localite.target", "old")
    ch = dict(lc.file_import_settings(m, "/x.xml", SERIES[4], SERIES))
    assert ch["placement.localite.position"] == 5
    assert ch["placement.localite.target"] is None


def test_several_transducers_set_their_own_entry(schema):
    m = ConfigModel(schema)
    m.add_transducer()
    ch = dict(lc.file_import_settings(m, "/x.xml", SERIES[1], SERIES, ti=1))
    assert ch["placement.localite.target"] == ["", "L_M1"]
    m.set("placement.localite.target", ["R_M1", "L_M1"])
    ch = dict(lc.file_import_settings(m, "/x.xml", SERIES[4], SERIES, ti=0))
    assert ch["placement.localite.position"] == [5, 1]
    assert ch["placement.localite.target"] == ["", "L_M1"]             # by number for transducer 1
