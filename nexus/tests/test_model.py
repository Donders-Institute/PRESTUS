"""ConfigModel mirrors load_parameters' merge (mergestruct + restore_struct_arrays)."""

from prestus_nexus import yamlio
from prestus_nexus.model import ConfigModel, mergestruct, values_equal


def test_mergestruct_semantics():
    a = {"x": 1, "g": {"p": 1, "q": 2}, "e": [], "s": "a"}
    b = {"g": {"q": 3}, "e": {"new": 1}, "s": ""}
    assert mergestruct(a, b) == {"x": 1, "g": {"p": 1, "q": 3}, "e": {"new": 1}, "s": ""}


def test_values_equal_follows_matlab_types():
    assert values_equal(1, 1.0)
    assert not values_equal(True, 1)
    assert not values_equal("", [])
    assert values_equal([1, 2.0], [1.0, 2])


def test_setting_baseline_removes_override(schema):
    m = ConfigModel(schema)
    m.set("grid.resolution_mm", 1.0)
    assert m.override_dict() == {"grid": {"resolution_mm": 1.0}}
    m.set("grid.resolution_mm", 0.5)
    assert m.override_dict() == {}


def test_later_transducers_inherit_transducer_one(schema):
    m = ConfigModel(schema)
    m.set("transducer.type", "annular")
    m.set("transducer.annular.elem_n", 4)
    m.add_transducer()
    m.set("transducer.trans_pos", [1, 2, 3], 1)
    assert m.get("transducer.type", 1) == "annular"          # inherited
    m.set("transducer.type", "annular", 1)                     # equals inherited -> no override
    assert m.tx[1] == {"trans_pos": [1, 2, 3]}
    tx = m.effective()["transducer"]
    assert tx[1]["annular"]["elem_n"] == 4
    out = m.override_dict()["transducer"]
    assert isinstance(out, list) and len(out) == 2


def test_override_order_follows_default(schema):
    m = ConfigModel(schema)
    m.set("timing.pd", 0.02)
    m.set("subject_id", 5)
    m.set("path.anat", "/data")
    assert list(m.override_dict()) == ["subject_id", "path", "timing"]


def test_unknown_and_deprecated_keys(schema):
    m = ConfigModel(schema)
    m.load_dict({"grid": {"resolutoin_mm": 1},
                 "transducer": {"annular": {"elem_amp": 1}},
                 "calibration": {"opt_method": "FEXminimize"}})
    found = {(i.level, i.path): i.message for i in m.unknown_keys()}
    assert "resolution_mm" in found[("warning", "grid.resolutoin_mm")]
    assert "elem_p_amp" in found[("warning", "transducer.annular.elem_amp")]
    assert ("info", "calibration.opt_method") in found


def test_loaded_keys_are_kept_until_pruned(schema):
    m = ConfigModel(schema)
    m.load_dict({"grid": {"resolution_mm": 0.5, "min_ppw": 8}})
    assert m.override_dict() == {"grid": {"resolution_mm": 0.5, "min_ppw": 8}}
    assert m.prune_defaults() == 1
    assert m.override_dict() == {"grid": {"min_ppw": 8}}


def test_save_and_reload(schema, tmp_path):
    m = ConfigModel(schema)
    m.set("path.t1_pattern", "T1_%03d.nii.gz")
    m.set("hpc.timelimit", "12:00:00")
    m.set("transducer.focus_pos", [10.5, -3, 2])
    target = m.save(tmp_path / "config_x.yaml", header="test")
    again = ConfigModel.from_file(schema, target)
    assert again.override_dict() == m.override_dict()
    assert yamlio.load_file(target)["hpc"]["timelimit"] == "12:00:00"


def test_transducer_list_has_identical_top_level_keys(schema):
    # yaml.loadFile concatenates list elements into a struct array, which
    # fails unless every element has the same fields (see test_matlab_roundtrip)
    m = ConfigModel(schema)
    m.set("transducer.name", "left")
    m.set("transducer.annular.elem_n", 4)
    m.add_transducer()
    m.set("transducer.trans_pos", [1, 2, 3], 1)
    out = m.override_dict()["transducer"]
    assert set(out[0]) == set(out[1])
    assert out[1]["name"] == "left"          # inherited value, written explicitly
    assert out[1]["annular"] == {}           # group: merged from transducer 1
    assert out[0]["trans_pos"] == []         # set on transducer 2 only: the default on 1
    again = ConfigModel(schema)
    again.load_dict({"transducer": out})
    assert values_equal(again.effective(), m.effective())


def test_keys_of_later_transducers_reach_every_element(schema):
    """Every key any element sets, in every direction: plain values and groups, set on
    transducer 2 only, on 3 only, or on 2 and 3 (a sequence's timing_offset_s and timing)."""
    m = ConfigModel(schema)
    for _ in range(2):
        m.add_transducer()
    m.set("transducer.timing_offset_s", 80.0, 1)
    m.set("transducer.timing_offset_s", 160.0, 2)
    m.set("transducer.timing.ptrd", 4.0, 1)
    m.set("transducer.focal_distance_ep", 50, 2)
    out = m.override_dict()["transducer"]
    assert set(out[0]) == set(out[1]) == set(out[2])
    assert out[0]["timing_offset_s"] == [] and out[0]["timing"] == {} and out[2]["timing"] == {}
    assert out[1]["focal_distance_ep"] == []
    again = ConfigModel(schema)
    again.load_dict({"transducer": out})
    # an absent key and its written empty default mean the same (unset) to PRESTUS
    for ti in range(3):
        for key in ("transducer.timing_offset_s", "transducer.timing.ptrd",
                    "transducer.focal_distance_ep"):
            assert values_equal(again.get(key, ti), m.get(key, ti))


def test_transducer_name_title_and_tooltip(schema):
    m = ConfigModel(schema)
    m.set("transducer.serial", "")
    assert m.transducer_title(0) == "Transducer 1"
    assert "no device set" in m.transducer_tooltip(0)
    m.set("transducer.serial", "CTX-TEST")
    m.add_transducer()
    assert m.transducer_title(1) == "Transducer 2: CTX-TEST"   # inherited from transducer 1
    m.set("transducer.serial", "NS-TEST", 1)
    assert m.transducer_name(1) == "NS-TEST"
    assert m.transducer_tooltip(1).startswith("Transducer 2: NS-TEST")
