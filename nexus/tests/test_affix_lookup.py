"""Run outputs found by the affix they carry: as load_parameters writes it ('_desc-...'), else
io.output_affix as written (outputs of scripts that set the affix after load_parameters)."""

import os
import time

import pytest

from prestus_nexus import outputs
from prestus_nexus.study import affix_path, sanitize_affix


def _make(tmp_path):
    return lambda a: tmp_path / f"sub-009_layered_medium_stats{a}.json"


@pytest.mark.parametrize("present,expected", [
    (["_desc-kosciessa_k-plan"], "_desc-kosciessa_k-plan"),     # current runs
    (["_kosciessa_k-plan"], "_kosciessa_k-plan"),               # affix set after load_parameters
    (["_desc-kosciessa_k-plan", "_kosciessa_k-plan"], "_desc-kosciessa_k-plan"),
    ([], "_desc-kosciessa_k-plan"),                             # neither: the current form
])
def test_affix_path_prefers_the_load_parameters_form(tmp_path, present, expected):
    for a in present:
        _make(tmp_path)(a).write_text("{}")
    assert affix_path(_make(tmp_path), "_kosciessa_k-plan") == _make(tmp_path)(expected)


def test_affix_path_does_not_touch_the_disk_when_forms_agree(tmp_path, monkeypatch):
    from pathlib import Path
    calls = []
    monkeypatch.setattr(Path, "exists", lambda self: calls.append(self) or False)
    assert affix_path(_make(tmp_path), "_desc-a") == _make(tmp_path)("_desc-a")
    assert affix_path(_make(tmp_path), "") == _make(tmp_path)("")
    assert calls == []


def test_placement_record_and_profile_paths_follow_the_affix(tmp_path, monkeypatch, schema):
    from prestus_nexus.model import ConfigModel
    m = ConfigModel(schema)
    m.set("subject_id", 9)
    m.set("io.output_affix", "_kosciessa_k-plan")
    monkeypatch.setattr(outputs, "dir_output", lambda model: tmp_path)
    (tmp_path / "img").mkdir()
    (tmp_path / "sub-009_placement_kosciessa_k-plan.json").write_text("{}")
    (tmp_path / "img" / "sub-009_layered_axial_profile_kosciessa_k-plan.csv").write_text("")
    assert outputs.placement_record_path(m).name == "sub-009_placement_kosciessa_k-plan.json"
    prof, target = outputs.axial_profile_paths(m)
    assert prof.name == "sub-009_layered_axial_profile_kosciessa_k-plan.csv"
    assert target.name == "sub-009_layered_desc-kosciessa_k-plan_targeting.csv"   # missing: current form


def test_skull_stats_take_this_configs_file_not_the_newest(tmp_path):
    pytest.importorskip("PyQt6")
    from prestus_nexus.ui.pct_page import stats_candidates, stats_title
    own = tmp_path / "sub-009_layered_medium_stats_kosciessa_k-plan.json"
    other = tmp_path / "sub-009_layered_medium_stats_wiesinger_mc.json"
    own.write_text("{}")
    other.write_text("{}")
    t = time.time()
    os.utime(own, (t - 100, t - 100))
    os.utime(other, (t, t))                              # negative control: newer, another run's
    real, _ = stats_candidates(tmp_path, 9, "layered", "_kosciessa_k-plan")
    assert real == own
    own.unlink()
    real, _ = stats_candidates(tmp_path, 9, "layered", "_kosciessa_k-plan")
    assert real == other                                 # nothing of its own: the newest, labelled
    title = stats_title({"n_voxels": 1606731.0, "skull_model": "pct"}, other, own=False)
    assert "n=1606731 skull voxels" in title and "from another run" in title
    assert "another run" not in stats_title({"n_voxels": 5}, own, own=True)


def test_sanitize_affix_twin_of_load_parameters():
    assert sanitize_affix("_kosciessa_k-plan") == "_desc-kosciessa_k-plan"
    assert sanitize_affix("_desc-x") == "_desc-x" and sanitize_affix("") == ""
