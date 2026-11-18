import shutil
from pathlib import Path

import pytest

from prestus_nexus import checks
from prestus_nexus import simnibs as sn
from prestus_nexus.model import ConfigModel

INI = "# standard settings\n[preprocess]\ndenoise = false\n\n[samseg]\nthreads = 8\n"


def _conda_env(home: Path, name: str, version: str) -> Path:
    """A conda environment with SimNIBS: bin/charm is a console script of the environment."""
    env = home / ".conda" / "envs" / name
    (env / "bin").mkdir(parents=True)
    (env / "bin" / "charm").write_text("#!/usr/bin/env python\nfrom simnibs.cli.charm import main\n")
    site = env / "lib" / "python3.11" / "site-packages"
    (site / "simnibs").mkdir(parents=True)
    (site / "simnibs" / "charm.ini").write_text(INI)
    (site / f"simnibs-{version}.dist-info").mkdir()
    return env / "bin"


def _installer(home: Path, version: str) -> Path:
    """The SimNIBS installer layout: bin/charm wraps simnibs_env's python and charm.py."""
    root = home / "Applications" / f"SimNIBS-{version}"
    site = root / "simnibs_env" / "lib" / "python3.9" / "site-packages"
    (site / "simnibs" / "cli").mkdir(parents=True)
    (site / f"simnibs-{version}.dist-info").mkdir()
    (root / "bin").mkdir(parents=True)
    (root / "bin" / "charm").write_text(
        f'#! /bin/bash -e\n"{root}/simnibs_env/bin/python" -E -u "{site}/simnibs/cli/charm.py" "$@"\n')
    return root / "bin"


def _m2m(seg: Path, sid: int, version: str, settings: str = INI) -> Path:
    m2m = seg / f"m2m_sub-{sid:03d}"
    m2m.mkdir(parents=True)
    (m2m / "final_tissues.nii.gz").write_bytes(b"")
    (m2m / "charm_log.html").write_text(f"INFO: simnibs version {version}\n<br>more\n")
    (m2m / "settings.ini").write_text(settings)
    (m2m / "charm_report.html").write_text("<html></html>")
    return m2m


def test_install_info_reads_version_without_running(tmp_path):
    conda = sn.install_info(_conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"))
    assert conda.version == "4.6.0" and conda.supported
    assert conda.default_settings.name == "charm.ini"
    inst = sn.install_info(_installer(tmp_path, "4.0.0"))
    assert inst.version == "4.0.0" and inst.supported
    assert sn.install_info(tmp_path) is None                   # no charm


def test_simnibs_3_is_not_supported(tmp_path):
    old = sn.install_info(_conda_env(tmp_path, "simnibs_v3", "3.2.6"))
    assert old.supported is False
    assert sn.is_supported(None) is None


def test_find_installs_lists_several_newest_first(tmp_path, monkeypatch):
    monkeypatch.setattr(sn.shutil, "which", lambda _n: None)
    _conda_env(tmp_path, "simnibs_v4.5.0", "4.5.0")
    _conda_env(tmp_path, "PRESTUS_env_4.6.0", "4.6.0")
    _installer(tmp_path, "4.0.0")
    found = [i for i in sn.find_installs(tmp_path) if str(i.bin).startswith(str(tmp_path))]
    assert [i.version for i in found] == ["4.6.0", "4.5.0", "4.0.0"]


def test_settings_compare_ignores_comments_and_spacing(tmp_path):
    a, b, c = tmp_path / "a.ini", tmp_path / "b.ini", tmp_path / "c.ini"
    a.write_text(INI)
    b.write_text("[preprocess]\ndenoise=false   # off\n[samseg]\n  threads =  8\n")
    c.write_text(INI.replace("threads = 8", "threads = 4"))
    assert sn.settings_equal(a, b)
    assert not sn.settings_equal(a, c)


def test_differences_version_and_settings(tmp_path):
    inst = sn.install_info(_conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"))
    same = sn.segmentation_state(_m2m(tmp_path / "seg", 1, "4.6.0"))
    assert sn.differences(same, inst, None) == []
    other = sn.segmentation_state(_m2m(tmp_path / "seg", 2, "4.0.0", INI.replace("false", "true")))
    diffs = sn.differences(other, inst, None)
    assert any("4.0.0" in d for d in diffs) and any("settings" in d for d in diffs)
    custom = tmp_path / "custom.ini"
    custom.write_text(INI.replace("false", "true"))
    assert not any("settings" in d for d in sn.differences(other, inst, str(custom)))


def _layered(schema, tmp_path, bin_path, m2m_version=None, **extra):
    anat = tmp_path / "anat"
    anat.mkdir(exist_ok=True)
    (anat / "sub-001_T1w.nii.gz").write_bytes(b"")
    if m2m_version:
        _m2m(anat, 1, m2m_version)
    m = ConfigModel(schema)
    for k, v in {"subject_id": 1, "path.anat": str(anat), "path.sim": str(tmp_path / "sim"),
                 "simulation.medium": "layered", "startup.simnibs_bin_path": str(bin_path), **extra}.items():
        m.set(k, v)
    return m


def _issues(m, path, level):
    return [i.message for i in checks.subject_issues(m) if i.path == path and i.level == level]


def test_checks_reject_simnibs_3_and_missing_charm(schema, tmp_path):
    m = _layered(schema, tmp_path, _conda_env(tmp_path, "simnibs_v3", "3.2.6"))
    assert any("SimNIBS 4" in s for s in _issues(m, "startup.simnibs_bin_path", "error"))
    m.set("startup.simnibs_bin_path", str(tmp_path))
    assert any("No charm" in s for s in _issues(m, "startup.simnibs_bin_path", "error"))


def test_checks_warn_about_a_differing_segmentation(schema, tmp_path):
    m = _layered(schema, tmp_path, _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"), m2m_version="4.0.0")
    assert any("SimNIBS 4.0.0" in s for s in _issues(m, "path.seg", "warning"))
    m.set("io.overwrite_simnibs", 1)                           # re-run: the warning is about deleting it
    assert not _issues(m, "path.seg", "warning")
    assert _issues(m, "io.overwrite_simnibs", "warning")


def test_checks_settings_file(schema, tmp_path):
    m = _layered(schema, tmp_path, _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"), m2m_version="4.6.0")
    m.set("segmentation.settings_file", str(tmp_path / "missing.ini"))
    assert _issues(m, "segmentation.settings_file", "error")
    inside = tmp_path / "anat" / "m2m_sub-001" / "my.ini"
    inside.write_text(INI)
    m.set("segmentation.settings_file", str(inside))
    assert any("outside" in s for s in _issues(m, "segmentation.settings_file", "error"))


# ── installs on the cluster (scan over SSH) ──────────────────────────────

def _run_scan(home: Path) -> str:
    import subprocess
    from prestus_nexus import remote
    env = {"HOME": str(home), "PATH": "/usr/bin:/bin"}
    return subprocess.run(["bash", "-c", remote.SIMNIBS_SCAN], env=env, capture_output=True,
                          text=True, timeout=60).stdout


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_remote_scan_script_finds_installs(tmp_path):
    from prestus_nexus import remote
    b46 = _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0")
    (b46 / "pydeface").write_text("")
    _installer(tmp_path, "4.0.0")
    (tmp_path / "SimNIBS_v4.0.0").mkdir()      # postinstall link folder: same install via a symlink
    (tmp_path / "SimNIBS_v4.0.0" / "bin").symlink_to(tmp_path / "Applications" / "SimNIBS-4.0.0" / "bin")
    (tmp_path / "SimNIBS_copy").mkdir()
    (tmp_path / "SimNIBS_copy" / "bin").symlink_to(tmp_path / "Applications" / "SimNIBS-4.0.0" / "bin")
    found = remote.parse_simnibs_scan(_run_scan(tmp_path), "login.cluster")
    found = {str(i.bin): i for i in found if str(i.bin).startswith(str(tmp_path))}
    conda = found[str(b46)]
    assert conda.version == "4.6.0" and conda.host == "login.cluster" and conda.pydeface
    # the scan's charm.ini matches the local normalisation
    assert list(conda.settings_lines) == sn.normalised_settings(
        tmp_path / ".conda/envs/simnibs_v4.6.0/lib/python3.11/site-packages/simnibs/charm.ini")
    # the installer layout (charm wraps simnibs_env); not listed twice through the symlink
    vers = sorted(i.version for i in found.values())
    assert vers == ["4.0.0", "4.6.0"]


def test_remote_scan_unfinished_is_none():
    from prestus_nexus import remote
    assert remote.parse_simnibs_scan("ssh: connect to host x: Connection refused", "x") is None


def test_lookup_uses_installs_of_the_submission_host(tmp_path, monkeypatch):
    inst = sn.Install(sn.PurePosixPath("/home/u/.conda/envs/s46/bin"), "4.6.0", None, host="hpc",
                      settings_lines=("[samseg]", "threads=8"))
    sn.set_remote_installs("hpc", [inst])
    monkeypatch.setattr(sn, "remote_host", lambda: "hpc")
    assert sn.lookup("/home/u/.conda/envs/s46/bin/") == inst
    assert sn.Install.from_dict(inst.to_dict()) == inst
    monkeypatch.setattr(sn, "remote_host", lambda: None)            # local runs: not this install
    assert sn.lookup("/home/u/.conda/envs/s46/bin") is None
    sn.set_remote_installs("hpc", [])


def test_checks_accept_a_scanned_cluster_install(schema, tmp_path, monkeypatch):
    remote_bin = "/home/u/.conda/envs/simnibs_v4.6.0/bin"
    m = _layered(schema, tmp_path, remote_bin, m2m_version="4.0.0")
    monkeypatch.setattr(sn, "remote_host", lambda: "hpc")
    assert any("scan the cluster" in s for s in _issues(m, "startup.simnibs_bin_path", "warning"))   # not scanned yet
    sn.set_remote_installs("hpc", [sn.Install(sn.PurePosixPath(remote_bin), "4.6.0", None, host="hpc",
                                              pydeface=False)])
    try:
        assert not _issues(m, "startup.simnibs_bin_path", "warning")
        assert not _issues(m, "startup.simnibs_bin_path", "error")
        assert any("4.0.0" in s for s in _issues(m, "path.seg", "warning"))     # compared with the cluster's 4.6.0
        m.set("segmentation.deface", 1)
        assert any("pydeface" in s for s in _issues(m, "segmentation.deface", "warning"))
    finally:
        sn.set_remote_installs("hpc", [])
    monkeypatch.setattr(sn, "remote_host", lambda: None)            # local runs: the folder must exist here
    assert _issues(m, "startup.simnibs_bin_path", "error")


# ── settings that follow the SimNIBS environment ─────────────────────────

def _with_python(bin_dir: Path) -> Path:
    """Give a conda env fixture a real interpreter on charm's '#!' line and bundled libraries."""
    env = bin_dir.parent
    py = env / "bin" / "python3.11"
    py.write_text("")
    (bin_dir / "charm").write_text(f"#!{py}\nfrom simnibs.cli.charm import main\n")
    (env / "lib" / "python3.11" / "site-packages" / "simnibs" / "external" / "lib" / "linux").mkdir(parents=True)
    return py


def test_install_python_and_libraries(tmp_path):
    b = _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0")
    py = _with_python(b)
    inst = sn.install_info(b)
    assert inst.python == str(py)
    assert inst.ld_library_path.endswith("simnibs/external/lib/linux")
    assert Path(str(inst.env_root)) == b.parent
    # installer: the wrapper's quoted python
    ib = _installer(tmp_path, "4.0.0")
    ipy = ib.parent / "simnibs_env" / "bin" / "python"
    ipy.parent.mkdir(parents=True, exist_ok=True)
    ipy.write_text("")
    assert sn.install_info(ib).python == str(ipy)


def test_owns_and_linked_mismatches(tmp_path):
    a = sn.install_info(_conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"))
    other = _conda_env(tmp_path, "simnibs_v4.0.0", "4.0.0")
    inst = sn.install_info(_installer(tmp_path, "4.5.0"))
    assert a.owns(str(tmp_path / ".conda/envs/simnibs_v4.6.0/bin")) is True
    assert a.owns(str(other)) is False
    assert inst.owns(str(tmp_path / "Applications" / "SimNIBS-4.5.0")) is True     # install root (PlanTUS)
    values = {"hpc.ld_library_path": str(other.parent / "lib/python3.11/site-packages/simnibs/external/lib/linux"),
              "placement.plantus.env_path": str(other)}
    assert [k for k, _ in sn.linked_mismatches(values.get, a)] == list(sn.LINKED_KEYS)
    values["hpc.ld_library_path"] = "/opt/gcc/7.2.0/lib64"           # a system folder: not SimNIBS's
    assert [k for k, _ in sn.linked_mismatches(values.get, a)] == ["placement.plantus.env_path"]


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
def test_remote_scan_reports_python_and_libraries(tmp_path):
    from prestus_nexus import remote
    b = _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0")
    py = _with_python(b)
    found = [i for i in remote.parse_simnibs_scan(_run_scan(tmp_path), "hpc") if str(i.bin) == str(b)]
    assert found[0].python == str(py)
    assert found[0].ld_library_path.endswith("simnibs/external/lib/linux")


def test_checks_plantus_env_follows_simnibs(schema, tmp_path):
    b = _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0")
    _with_python(b)
    m = _layered(schema, tmp_path, b, m2m_version="4.6.0")
    m.set("placement.mode", "plantus")
    issues = [i for i in checks.all_issues(m) if i.path == "placement.plantus.env_path"]
    assert not [i for i in issues if i.level == "error"]
    m.set("placement.plantus.env_path", str(_conda_env(tmp_path, "simnibs_v4.0.0", "4.0.0")))
    assert any("another SimNIBS install" in i.message for i in checks.all_issues(m)
               if i.path == "placement.plantus.env_path")


def test_panel_relinks_on_environment_change(schema, tmp_path):
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])  # noqa: F841  (kept alive for the widget)
    from prestus_nexus.ui.simnibs_panel import SimnibsPanel
    a, b = _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"), _conda_env(tmp_path, "simnibs_v4.5.0", "4.5.0")
    m = _layered(schema, tmp_path, a)
    own = str(a.parent / "lib/python3.11/site-packages/simnibs/external/lib/linux")
    m.set("hpc.ld_library_path", own)                     # pointed into the old install
    m.set("placement.plantus.env_path", "/somewhere/custom")  # not the old install's: kept
    p = SimnibsPanel(m)
    p._choose(str(b))
    assert m.get("startup.simnibs_bin_path") == str(b)
    assert m.get("hpc.ld_library_path") in ("", None)
    assert m.get("placement.plantus.env_path") == "/somewhere/custom"


def test_segmentation_card_has_no_rerun_switch(schema, tmp_path):
    """Segment again is the button (segmentation_panel); io.overwrite_simnibs, set elsewhere,
    is flagged on the card instead of being a second switch."""
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication, QCheckBox
    app = QApplication.instance() or QApplication([])  # noqa: F841
    from prestus_nexus.ui.simnibs_panel import SimnibsPanel
    m = _layered(schema, tmp_path, _conda_env(tmp_path, "simnibs_v4.6.0", "4.6.0"), m2m_version="4.6.0")
    p = SimnibsPanel(m)
    assert not [c for c in p.findChildren(QCheckBox) if "segment" in c.text().lower()]
    p.refresh()
    off = p.state.text()
    m.set("io.overwrite_simnibs", 1)
    p.refresh()
    assert "segments again" in p.state.text() and "segments again" not in off
