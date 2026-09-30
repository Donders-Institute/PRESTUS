"""ConfigModel mirrors load_parameters' merge (mergestruct + restore_struct_arrays)."""

from prestus_planner import yamlio
from prestus_planner.model import ConfigModel, mergestruct, values_equal


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
    assert set(out[0]) - {"trans_pos"} <= set(out[1])
    assert out[1]["name"] == "left"          # inherited value, written explicitly
    assert out[1]["annular"] == {}           # group: merged from transducer 1
    again = ConfigModel(schema)
    again.load_dict({"transducer": out})
    assert values_equal(again.effective(), m.effective())
