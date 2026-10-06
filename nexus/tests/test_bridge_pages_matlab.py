"""End to end with MATLAB: the Localite and Calibration pages drive the real nexus_bridge
(synthetic Localite study; CTX_250_001 hydrophone profiles from local/profiles)."""

import glob
import os
import shutil
import subprocess
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6")

from PyQt6.QtCore import Qt  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from prestus_nexus import bridge as br  # noqa: E402
from prestus_nexus import yamlio  # noqa: E402
from prestus_nexus.equipment import EquipmentLibrary  # noqa: E402
from prestus_nexus.model import ConfigModel  # noqa: E402
from prestus_nexus.session import Session  # noqa: E402
from prestus_nexus.ui.main_window import CALIBRATION, LOCALITE, MainWindow  # noqa: E402


def _find_matlab():
    cand = [os.environ.get("PRESTUS_MATLAB"), shutil.which("matlab")]
    cand += sorted(glob.glob("/Applications/MATLAB_R20*.app/bin/matlab"), reverse=True)
    return next((c for c in cand if c and Path(c).exists()), None)


MATLAB = _find_matlab()
pytestmark = [pytest.mark.matlab, pytest.mark.skipif(MATLAB is None, reason="MATLAB not available")]


@pytest.fixture(scope="module")
def app():
    a = QApplication.instance() or QApplication([])
    yield a
    br.shutdown_shared()


@pytest.fixture(scope="module")
def study(prestus_root, tmp_path_factory):
    root = tmp_path_factory.mktemp("study")
    body = ("addpath(fullfile('{p}', 'tests', 'fixtures')); "
            "rz = @(d) [cosd(d) -sind(d) 0; sind(d) cosd(d) 0; 0 0 1]; "
            "make_synthetic_localite_export('{r}/raw', 'sub-001', 'ses-01', "
            "{{[rz(20), [40; 90; 140]; 0 0 0 1], [rz(160), [140; 95; 135]; 0 0 0 1]}}, [1; 1; 0], "
            "{{'LEFT PUL', 'RIGHT PUL'}})").format(p=prestus_root, r=root)
    subprocess.run([MATLAB, "-batch", body], check=True, capture_output=True, timeout=300)
    return root


def _wait(app, cond, timeout=300.0):
    t0 = time.monotonic()
    while not cond() and time.monotonic() - t0 < timeout:
        app.processEvents()
        time.sleep(0.02)
    return cond()


def test_localite_and_calibration_pages_with_matlab(app, schema, prestus_root, study, monkeypatch):
    monkeypatch.setenv("PRESTUS_MATLAB", MATLAB)
    m = ConfigModel(schema, EquipmentLibrary(prestus_root))
    m.set("path.localite_raw", str(study / "raw"))
    m.set("path.localite_post", str(study / "post"))
    m.set("transducer.serial", "CTX_250_001")
    m.set("transducer.combo.ds_serial", "SC_105_010")
    m.set("transducer.library_path", str(study / "lib"))
    m.save(study / "config_study.yaml")
    w = MainWindow(m, Session(study / "session.yaml"), prestus_root)
    try:
        loc = w.show_special(LOCALITE)
        loc.detect()
        assert _wait(app, lambda: loc.write_btn.isEnabled() or "✗" in loc.status.text()), loc.status.text()
        assert loc.table.rowCount() == 2, loc.status.text()
        assert loc.table.item(0, 3).text() == "LEFT PUL"
        loc.write()
        assert _wait(app, lambda: loc.use_btn.isEnabled() or "✗" in loc.status.text()), loc.status.text()
        assert (study / "post" / "sub-001" / "ses-01" / "InstrumentMarker_sub-001_ses-01.json").is_file()

        prof = prestus_root / "local" / "profiles" / "CTX_250_001"
        if not prof.is_dir():
            pytest.skip("no hydrophone profiles in local/profiles (calibration part)")
        cal = w.show_special(CALIBRATION)
        assert _wait(app, lambda: cal.measured.count() > 0 and not cal._busy), cal.log.toPlainText()
        for i in range(cal.measured.count()):
            cal.measured.item(i).setCheckState(Qt.CheckState.Checked if i == 6 else Qt.CheckState.Unchecked)
        cal.nstarts.setValue(2)
        cal.fit()
        assert _wait(app, lambda: cal.has_fit or "✗" in cal.log.toPlainText()), cal.log.toPlainText()
        cal.save()
        entry = study / "lib" / "CTX_250_001_SC_105_010.yaml"
        assert _wait(app, lambda: entry.is_file()), cal.log.toPlainText()
        meta = yamlio.load_file(entry)["meta"]
        assert meta["calibration_type"] == "user_empirical" and meta["fit"]["tool"] == "nexus"
        # reopening the entry restores its fit settings
        cal.opened = None
        cal.nstarts.setValue(9)
        cal._open()
        assert _wait(app, lambda: not cal._busy and cal.nstarts.value() == 2), cal.log.toPlainText()
    finally:
        w.dirty = False
        w.close()
