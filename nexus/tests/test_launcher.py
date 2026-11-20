"""Launcher script, MATLAB commands and the Run page."""

import os
import stat
import subprocess
from pathlib import Path

import pytest

from prestus_nexus import launcher
from prestus_nexus.model import ConfigModel


def test_launcher_script(tmp_path, prestus_root):
    study = tmp_path / "it's" / "config_pilot.yaml"
    study.parent.mkdir()
    study.write_text("subject_id: 1\n")
    path = launcher.write_launcher(study, prestus_root)
    assert path.name == "run_pilot.m"
    text = path.read_text()
    assert "load_parameters('config_pilot.yaml', '" in text
    assert "it''s" in text                           # MATLAB quote escaping
    assert "parameters.platform = 'matlab';" in text
    assert "prestus_pipeline_start(parameters);" in text
    assert "platform" not in launcher.launcher_script(study, prestus_root, local=False).split("parameters =")[1]


def test_plantus_script(tmp_path, prestus_root):
    study = tmp_path / "config_pilot.yaml"
    study.write_text("subject_id: 1\n")
    path = launcher.write_plantus_script(study, prestus_root, interactive=True)
    assert path.name == "plantus_pilot.m"
    text = path.read_text()
    assert "'placement.plantus.interactive', 1;" in text
    assert "'simulation.interactive', 0;" in text and "'placement.reuse', 'never'" in text
    assert "parameters.modules.run_transducer_placement = 1;" in text
    assert "'placement.plantus.interactive', 0;" in launcher.plantus_script(study, prestus_root, False)


def test_commands_wrap_sentinels(tmp_path, prestus_root):
    argv = launcher.validate_command("/x/matlab", tmp_path / "config_a.yaml", prestus_root)
    assert argv[:2] == ["/x/matlab", "-batch"]
    assert launcher.OK_SENTINEL in argv[2] and launcher.FAIL_SENTINEL in argv[2]
    run = launcher.run_command("/x/matlab", tmp_path / "run_a.m")
    assert run[1] == "-sd" and "run('" in run[-1]


MATLAB = launcher.find_matlab()


@pytest.mark.matlab
@pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")
@pytest.mark.parametrize("good", [True, False])
def test_validate_in_matlab(schema, prestus_root, tmp_path, good):
    m = ConfigModel(schema)
    for k, v in {"subject_id": 2, "path.sim": str(tmp_path), "simulation.medium": "water",
                 "transducer.type": "annular", "transducer.freq_hz": 5e5,
                 "transducer.focal_distance_ep": 60,
                 "transducer.annular.elem_n": 2, "transducer.annular.elem_id_mm": [0, 30],
                 "transducer.annular.elem_od_mm": [29, 60] if good else [29],   # wrong length
                 "transducer.annular.elem_p_amp": 1e5, "transducer.annular.elem_phase_deg": [0, 0],
                 "transducer.annular.curv_radius_mm": 63.2}.items():
        m.set(k, v)
    study = m.save(tmp_path / "config_v.yaml")
    proc = subprocess.run(launcher.validate_command(MATLAB, study, prestus_root),
                          capture_output=True, text=True, timeout=600, cwd=tmp_path)
    out = proc.stdout + proc.stderr
    if good:
        assert proc.returncode == 0 and launcher.OK_SENTINEL in out, out[-2000:]
        assert "Loaded: subject 2, 1 transducer(s), medium water" in out
    else:
        assert proc.returncode != 0 and launcher.FAIL_SENTINEL in out
        assert "elem_od_mm length must match elem_n" in out


# ── Run page with a fake MATLAB ─────────────────────────────────────────

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture
def fake_matlab(tmp_path):
    def make(ok: bool) -> str:
        script = tmp_path / ("matlab_ok" if ok else "matlab_fail")
        body = (f"echo 'Loading...'\necho 'Warning: something odd'\necho '{launcher.OK_SENTINEL}'\nexit 0\n"
                if ok else f"echo 'Error using load_parameters'\necho '{launcher.FAIL_SENTINEL}: boom' >&2\nexit 1\n")
        script.write_text("#!/bin/sh\n" + body)
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        return str(script)
    return make


@pytest.mark.posix
@pytest.mark.parametrize("ok", [True, False])
def test_run_page_streams_and_reports(schema, prestus_root, tmp_path, fake_matlab, ok):
    pytest.importorskip("PyQt6")
    from PyQt6.QtCore import QEventLoop, QTimer
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841 (keep a reference)
    m = ConfigModel(schema)
    m.set("path.sim", str(tmp_path))
    m.save(tmp_path / "config_r.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    page = window._page(RUN)
    page.matlab.setText(fake_matlab(ok))
    results = []
    loop = QEventLoop()
    page.finished.connect(lambda task, success: (results.append((task, success)), loop.quit()))
    QTimer.singleShot(20000, loop.quit)
    page.validate()
    loop.exec()
    assert results == [("Validate", ok)]
    log = page.log.toPlainText()
    assert ("Warning: something odd" in log) if ok else ("boom" in log)
    assert page.proc is None and page.b_validate.isEnabled()


def test_matlab_gui_commands(tmp_path, prestus_root):
    study = tmp_path / "config_a.yaml"
    argv = launcher.gui_command("/x/matlab", prestus_root, launcher.neuronav_gui_call(study))
    assert argv[:3] == ["/x/matlab", "-nosplash", "-r"]
    assert f"neuronav_gui('{study.resolve()}')" in argv[3]
    assert launcher.calibration_gui_call("/lib/it's") == \
        "calibration_gui(struct('OutputPath', '/lib/it''s'))"
    assert launcher.calibration_gui_call(None) == "calibration_gui()"


@pytest.mark.matlab
@pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")
def test_matlab_gui_entry_points_exist(prestus_root):
    # the GUIs themselves need a display; check the functions resolve and take our inputs
    body = (f"addpath(genpath(fullfile('{prestus_root}', 'functions'))); "
            "assert(~isempty(which('neuronav_gui')) && nargin('neuronav_gui') == 1); "
            "assert(~isempty(which('calibration_gui')) && nargin('calibration_gui') == 1)")
    proc = subprocess.run([MATLAB, "-batch", body], capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]


@pytest.mark.matlab
@pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")
def test_plantus_script_loads_in_matlab(schema, prestus_root, tmp_path):
    """The placement-only script's overrides load and leave only placement on."""
    m = ConfigModel(schema)
    for k, v in {"subject_id": 2, "path.sim": str(tmp_path), "simulation.medium": "water",
                 "simulation.interactive": 1, "placement.mode": "plantus"}.items():
        m.set(k, v)
    study = m.save(tmp_path / "config_p.yaml")
    script = launcher.plantus_script(study, prestus_root, interactive=True).replace(
        "prestus_pipeline_start(parameters);",
        "mods = parameters.modules; on = {}; for f = fieldnames(mods)', "
        "if isnumeric(mods.(f{1})) && isscalar(mods.(f{1})) && mods.(f{1}) == 1, on{end+1} = f{1}; end, end\n"
        "fprintf('ON=%s PT=%d SIM=%d REUSE=%s\\n', strjoin(on, ','), "
        "parameters.placement.plantus.interactive, parameters.simulation.interactive, "
        "parameters.placement.reuse);")
    path = tmp_path / "plantus_check.m"
    path.write_text(script)
    proc = subprocess.run(launcher.run_command(MATLAB, path), capture_output=True, text=True,
                          timeout=600, cwd=tmp_path)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0 and launcher.OK_SENTINEL in out, out[-2000:]
    assert "ON=run_transducer_placement PT=1 SIM=0 REUSE=never" in out, out[-2000:]



def test_plantus_submit_script_keeps_platform(tmp_path, prestus_root):
    study = tmp_path / "config_pilot.yaml"
    study.write_text("subject_id: 1\n")
    local = launcher.plantus_script(study, prestus_root, False)
    job = launcher.plantus_script(study, prestus_root, False, local=False)
    assert "parameters.platform = 'matlab'" in local and "parameters.platform" not in job
    assert launcher.PLACEMENT_ONLY in job
    assert launcher.matlab_value(0) == "0" and launcher.matlab_value("never") == "'never'"
    assert launcher.matlab_value([]) == "[]" and launcher.matlab_value(True) == "true"


def test_kwave_drive_submit_script(tmp_path):
    """Submitting the k-Wave drive check calls kwave_drive_submit with the script, fit and config."""
    from prestus_nexus import launcher as L
    script = tmp_path / "kwave_drive_combo.m"
    script.write_text("% drive check\n")
    fit = tmp_path / "combo_fit.mat"
    cfg = tmp_path / "study.yaml"
    path = L.write_kwave_drive_submit(script, fit, cfg, tmp_path / "PRESTUS")
    assert path == script.with_name("submit_kwave_drive_combo.m")
    text = path.read_text()
    assert f"kwave_drive_submit('{script.resolve()}', '{fit.resolve()}', '{cfg}')" in text
    assert "addpath(genpath(fullfile(" in text
    body = L.kwave_drive_submit_body("/p/s.m", "/p/f.mat", None, "/p/PRESTUS")
    assert body.endswith("kwave_drive_submit('/p/s.m', '/p/f.mat', '')")
