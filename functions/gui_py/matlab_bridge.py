"""MATLAB bridge: run PRESTUS pipeline from Python.

Priority order:
  1. matlab.engine  (bidirectional, live stdout)
  2. subprocess     (matlab -batch, suitable for HPC dispatch)

The bridge runs the pipeline in a background QThread and emits
log lines and a finished/error signal back to the GUI thread.
"""

from __future__ import annotations
import os
import sys
import tempfile
import subprocess
import threading
from pathlib import Path

from PyQt6.QtCore import QThread, pyqtSignal, QObject


class MatlabWorker(QThread):
    """Runs prestus_pipeline_start in a background thread.

    Signals
    -------
    log_line : str      — one line of MATLAB output
    finished : bool     — True = success, False = error
    error    : str      — error message when finished(False)
    """

    log_line = pyqtSignal(str)
    finished = pyqtSignal(bool)
    error    = pyqtSignal(str)

    def __init__(self, config_path: str, prestus_root: str, parent=None):
        super().__init__(parent)
        self._config_path  = config_path
        self._prestus_root = prestus_root
        self._use_engine   = _matlab_engine_available()

    def run(self):
        try:
            if self._use_engine:
                self._run_engine()
            else:
                self._run_subprocess()
            self.finished.emit(True)
        except Exception as exc:
            self.error.emit(str(exc))
            self.finished.emit(False)

    # ── matlab.engine path ────────────────────────────────────────────────────

    def _run_engine(self):
        import matlab.engine
        self.log_line.emit("▶ Starting MATLAB engine…")
        eng = matlab.engine.start_matlab()
        try:
            eng.addpath(self._prestus_root, nargout=0)
            self.log_line.emit("▶ Running load_parameters…")
            params = eng.load_parameters(self._config_path, nargout=1)
            self.log_line.emit("▶ Running prestus_pipeline_start…")

            # Capture stdout by redirecting via a StringIO buffer
            import io
            out_buf = io.StringIO()
            eng.prestus_pipeline_start(params, nargout=0, stdout=out_buf, stderr=out_buf)
            for line in out_buf.getvalue().splitlines():
                self.log_line.emit(line)

            # Also tail log file if it exists
            log_file = None
            try:
                log_file = eng.eval("parameters.io.log_file", nargout=1)
            except Exception:
                pass
            if log_file:
                self._tail_log(log_file)
        finally:
            eng.quit()

    # ── subprocess path ───────────────────────────────────────────────────────

    def _run_subprocess(self):
        matlab_bin = _find_matlab_bin()
        if matlab_bin is None:
            raise RuntimeError(
                "MATLAB not found on PATH and matlab.engine is not available.\n"
                "Set startup.matlab_bin_path in your config or install the "
                "MATLAB Engine API for Python."
            )

        script = (
            f"addpath(genpath('{self._prestus_root}')); "
            f"parameters = load_parameters('{self._config_path}'); "
            f"prestus_pipeline_start(parameters);"
        )
        cmd = [matlab_bin, "-batch", script]
        self.log_line.emit(f"▶ Subprocess: {' '.join(cmd[:2])} -batch <script>")

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        for line in proc.stdout:
            self.log_line.emit(line.rstrip())
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError(f"MATLAB exited with code {proc.returncode}")

    def _tail_log(self, log_file: str):
        if not os.path.isfile(log_file):
            return
        with open(log_file) as f:
            for line in f:
                line = line.rstrip()
                if line:
                    self.log_line.emit(line)


class CalibrationWorker(QThread):
    """Runs calibration_pipeline_start in a background thread."""

    log_line = pyqtSignal(str)
    finished = pyqtSignal(bool)
    error    = pyqtSignal(str)

    def __init__(self, config_path: str, prestus_root: str, parent=None):
        super().__init__(parent)
        self._config_path  = config_path
        self._prestus_root = prestus_root

    def run(self):
        try:
            matlab_bin = _find_matlab_bin()
            if matlab_bin is None:
                raise RuntimeError("MATLAB not found on PATH.")
            script = (
                f"addpath(genpath('{self._prestus_root}')); "
                f"params = load_parameters('{self._config_path}'); "
                f"calibration_pipeline_start(params);"
            )
            cmd = [matlab_bin, "-batch", script]
            self.log_line.emit(f"▶ Running calibration via MATLAB…")
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            for line in proc.stdout:
                self.log_line.emit(line.rstrip())
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"MATLAB exited with code {proc.returncode}")
            self.finished.emit(True)
        except Exception as exc:
            self.error.emit(str(exc))
            self.finished.emit(False)


class IngestWorker(QThread):
    """Runs neuronav_ingest_markers in a background thread."""

    log_line = pyqtSignal(str)
    finished = pyqtSignal(bool, str)   # success, xml_path
    error    = pyqtSignal(str)

    def __init__(self, prestus_root: str, sub_id: str, ses_id: str,
                 raw_path: str, out_path: str, parent=None):
        super().__init__(parent)
        self._prestus_root = prestus_root
        self._sub_id   = sub_id
        self._ses_id   = ses_id
        self._raw_path = raw_path
        self._out_path = out_path

    def run(self):
        try:
            matlab_bin = _find_matlab_bin()
            if matlab_bin is None:
                raise RuntimeError("MATLAB not found on PATH.")
            script = (
                f"addpath(genpath('{self._prestus_root}')); "
                f"params.path.localite_raw = '{self._raw_path}'; "
                f"params.path.localite_post = '{self._out_path}'; "
                f"params.io.overwrite_files = 'always'; "
                f"tmap.series_index = 1; tmap.transducer_id = 1; tmap.target_name = 'target'; "
                f"[xml, ~] = neuronav_ingest_markers(params, '{self._sub_id}', '{self._ses_id}', tmap); "
                f"disp(['XML_PATH:' xml]);"
            )
            proc = subprocess.Popen(
                [matlab_bin, "-batch", script],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1
            )
            xml_path = ""
            for line in proc.stdout:
                line = line.rstrip()
                self.log_line.emit(line)
                if line.startswith("XML_PATH:"):
                    xml_path = line[9:]
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"MATLAB exited with code {proc.returncode}")
            self.finished.emit(True, xml_path)
        except Exception as exc:
            self.error.emit(str(exc))
            self.finished.emit(False, "")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _matlab_engine_available() -> bool:
    try:
        import matlab.engine
        return True
    except ImportError:
        return False


def _find_matlab_bin() -> str | None:
    import shutil
    candidate = shutil.which("matlab")
    if candidate:
        return candidate
    for search in (
        "/Applications/MATLAB_R2024b.app/bin/matlab",
        "/Applications/MATLAB_R2024a.app/bin/matlab",
        "/Applications/MATLAB_R2023b.app/bin/matlab",
        "/usr/local/bin/matlab",
    ):
        if os.path.isfile(search):
            return search
    return None
