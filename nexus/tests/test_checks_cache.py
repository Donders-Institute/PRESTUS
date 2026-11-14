"""checks.subject_issues reuses its file probes while the settings it read are unchanged."""

from prestus_nexus import checks
from prestus_nexus.model import ConfigModel


def _layered(schema, tmp_path):
    anat = tmp_path / "anat"
    anat.mkdir()
    m = ConfigModel(schema)
    for key, value in {"subject_id": 3, "path.anat": str(anat), "path.sim": str(tmp_path / "sim"),
                       "simulation.medium": "layered"}.items():
        m.set(key, value)
    return m, anat


def test_unrelated_edit_does_no_file_io(schema, tmp_path, monkeypatch):
    m, _ = _layered(schema, tmp_path)
    first = checks.subject_issues(m)
    calls = []
    real = checks._subject_issues
    monkeypatch.setattr(checks, "_subject_issues", lambda model: calls.append(1) or real(model))
    m.set("transducer.n_elements", 7)               # not read by subject_issues
    assert checks.subject_issues(m) == first and not calls
    m.set("path.t1_pattern", "sub-%1$03d_other.nii.gz")   # read by subject_issues
    checks.subject_issues(m)
    assert calls == [1]


def test_clear_status_cache_forces_a_new_probe(schema, tmp_path):
    m, anat = _layered(schema, tmp_path)
    assert ("path.t1_pattern", "error") in {(i.path, i.level) for i in checks.subject_issues(m)}
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    checks.clear_status_cache()
    assert ("path.t1_pattern", "info") in {(i.path, i.level) for i in checks.subject_issues(m)}


def test_cache_expires(schema, tmp_path, monkeypatch):
    m, anat = _layered(schema, tmp_path)
    checks.subject_issues(m)
    (anat / "sub-003_T1w.nii.gz").write_bytes(b"")
    monkeypatch.setattr(checks, "_STATUS_TTL_S", 0.0)
    assert ("path.t1_pattern", "info") in {(i.path, i.level) for i in checks.subject_issues(m)}
