"""Project copies of config_default.yaml: layering onto the toolbox default (as
load_parameters does), stale-copy detection, and refreshing the copy."""

from pathlib import Path

from prestus_nexus import checks, yamlio
from prestus_nexus.model import ConfigModel
from prestus_nexus.project_defaults import carry_pinned_values, replace_project_default
from prestus_nexus.schema import Schema


def _project(tmp_path: Path) -> Path:
    project = tmp_path / "config_default.yaml"
    project.write_text("simulation:\n  medium: water\ngrid:\n  resolution_mm: 0.75\n", encoding="utf-8")
    return project


def test_toolbox_default_is_not_layered(prestus_root):
    s = Schema(prestus_root / "config" / "config_default.yaml",
               toolbox_path=prestus_root / "config" / "config_default.yaml")
    assert s.project_path is None
    assert s.stale_keys == [] and s.pinned == {}


def test_project_copy_is_merged_onto_toolbox(tmp_path, prestus_root):
    toolbox = prestus_root / "config" / "config_default.yaml"
    s = Schema(_project(tmp_path), toolbox_path=toolbox)
    assert s.project_path == tmp_path / "config_default.yaml"
    assert s.defaults["grid"]["resolution_mm"] == 0.75                 # project pin wins
    assert "cache_signature_check" in s.defaults["io"]                 # toolbox key kept
    assert "io.cache_signature_check" in s.stale_keys
    assert s.pinned["grid.resolution_mm"][0] == 0.75
    assert "io.cache_signature_check" in s                             # known to the schema


def test_stale_copy_is_an_issue(tmp_path, prestus_root):
    s = Schema(_project(tmp_path), toolbox_path=prestus_root / "config" / "config_default.yaml")
    msgs = [i.message for i in checks.defaults_issues(ConfigModel(s))]
    assert len(msgs) == 1 and "Update project defaults" in msgs[0]


def test_refresh_keeps_effective_values(tmp_path, prestus_root):
    toolbox = prestus_root / "config" / "config_default.yaml"
    project = _project(tmp_path)
    model = ConfigModel(Schema(project, toolbox_path=toolbox))
    model.set("subject_id", 5)
    before = model.effective()

    written = carry_pinned_values(model, model.schema.pinned)
    assert "grid.resolution_mm" in written
    study = tmp_path / "config_study.yaml"
    model.save(study)
    backup = replace_project_default(project, toolbox)
    assert backup.is_file() and "resolution_mm: 0.75" in backup.read_text(encoding="utf-8")
    assert project.read_text(encoding="utf-8") == toolbox.read_text(encoding="utf-8")

    reloaded = ConfigModel.from_file(Schema(project, toolbox_path=toolbox), study)
    assert reloaded.schema.stale_keys == [] and reloaded.schema.pinned == {}
    after = reloaded.effective()
    assert after["grid"]["resolution_mm"] == before["grid"]["resolution_mm"] == 0.75
    assert after["simulation"]["medium"] == before["simulation"]["medium"] == "water"
    assert yamlio.load_file(study)["grid"]["resolution_mm"] == 0.75


def test_carry_does_not_overwrite_study_values(tmp_path, prestus_root):
    model = ConfigModel(Schema(_project(tmp_path),
                               toolbox_path=prestus_root / "config" / "config_default.yaml"))
    model.overrides.setdefault("grid", {})["resolution_mm"] = 1.0
    assert "grid.resolution_mm" not in carry_pinned_values(model, model.schema.pinned)
    assert model.overrides["grid"]["resolution_mm"] == 1.0
