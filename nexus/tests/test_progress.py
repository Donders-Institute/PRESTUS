"""Run progress from PRESTUS logs (progress.py), job logs (jobs.read_progress), and its
display on the Run page and in the jobs view."""

import os
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from prestus_nexus import jobs
from prestus_nexus import progress as P

T0 = 1_000_000.0


def _solver(name, stage="acoustic", t=T0):
    return f"▶ prestus_pipeline [#1] 💿/x\n▶ {stage} [#1] 💿/x\n[progress] {stage} 0/1 @{t} solver {name}\n"


# One log shape per way a solver reports progress: (solver, body, now, fraction, eta, estimated)
SOLVER_LOGS = {
    "kwave_cpp": ("cpp_gpu",
                  "| Progress |  Elapsed time  |  Time to go  |  Est. finish time  |\n"
                  "|    45%   |      12.345s  |     15.000s  |  19/11/26 10:11:12 |\n",
                  T0 + 20, 0.45, 15.0, False),
    "kwave_mlx": ("cpp_mlx", "|     5%   |       1.000s  |     19.000s  |  19/11/26 10:11:12 |\n",
                  T0 + 2, 0.05, 19.0, False),
    "kwave_matlab": ("matlab_gpu",
                     "  precomputation completed in 1min 0s\n  starting time loop...\n"
                     "  estimated simulation time 10min 0s...\n",
                     T0 + 60 + 300, 0.5, 300.0, True),
    "kwave_elastic": ("kwave_elastic",
                      "  precomputation completed in 2.5s\n  estimated simulation time 2hours 0min 0s...\n",
                      T0 + 2.5 + 3600, 0.5, 3600.0, True),
    "jwave_time_domain": ("jwave_pseudo_spectral",
                          "[run_jwave_sim]   step  250/1000\n[run_jwave_sim]   estimated remaining time: 30.0 s\n",
                          T0 + 10, 0.25, 30.0, False),
    "thermal_steps": ("kWaveDiffusion",
                      f"[progress] thermal 10/100 @{T0 + 10} t=1.00s Tmax=37.50C\n"
                      f"[progress] thermal 20/100 @{T0 + 20} t=2.00s Tmax=38.10C\n",
                      T0 + 25, 0.2, 75.0, False),
    "babelbrain_acoustic": ("babelbrain",
                            f"[progress] acoustic 101/1000 @{T0 + 10} FDTD run 1\n"
                            f"[progress] acoustic 501/1000 @{T0 + 50} FDTD run 1\n",
                            T0 + 50, 0.501, 49.9, False),
    "babelbrain_thermal": ("babelbrain",
                           f"[progress] thermal 10/100 @{T0 + 10} Tnow=37.20C\n"
                           f"[progress] thermal 20/100 @{T0 + 20} Tnow=38.10C\n",
                           T0 + 20, 0.2, 80.0, False),
}


@pytest.mark.parametrize("kind", sorted(SOLVER_LOGS))
def test_each_solver_kind_gives_fraction_and_time_left(kind):
    name, body, now, frac, eta, estimated = SOLVER_LOGS[kind]
    stage = "thermal" if "thermal" in kind else "acoustic"
    p = P.RunProgress.from_text(_solver(name, stage) + body)
    est = p.estimate(now)
    assert est.fraction == pytest.approx(frac, abs=1e-6)
    assert est.eta_s == pytest.approx(eta, abs=1e-6)
    assert est.estimated is estimated
    assert P.SOLVER_NAMES[name] in est.text


@pytest.mark.parametrize("name", sorted(P.SILENT_SOLVERS))
def test_silent_solvers_show_elapsed_time_only(name):
    est = P.RunProgress.from_text(_solver(name)).estimate(T0 + 90)
    assert est.fraction is None and est.eta_s is None
    assert "1 min 30 s elapsed" in est.text and "reports no progress" in est.text


def test_feed_in_pieces_equals_whole_text():
    name, body, now, *_ = SOLVER_LOGS["kwave_cpp"]
    text = _solver(name) + body
    whole = P.RunProgress.from_text(text)
    pieces = P.RunProgress()
    for i in range(0, len(text), 7):              # chunks split lines anywhere
        pieces.feed(text[i:i + 7])
    pieces.flush()
    assert pieces.estimate(now) == whole.estimate(now)


def test_stages_and_solver_end_with_their_stage():
    text = _solver("cpp_cpu") + "|    50%   |   1.0s  |   1.0s  | x |\n⏱ acoustic   12.0s | ΔRAM\n▶ acoustic_analysis [#1] x\n"
    p = P.RunProgress.from_text(text)
    assert p.stage == "acoustic_analysis" and p.solver is None
    assert "acoustic" in p.done_stages
    assert p.estimate(T0).fraction is None        # busy: no solver in this stage


def test_solver_outside_a_stage_is_cleared_by_the_next_stage():
    text = ("▶ prestus_pipeline [#1] x\n⏱ source  1.0s\n[progress] acoustic 0/1 @1 solver cpp_gpu\n"
            "|   100%   |   1.0s  |   0.0s  | x |\n▶ acoustic [#1] x\n")
    p = P.RunProgress.from_text(text)
    assert p.stage == "acoustic" and p.solver is None


def test_finished_failed_and_live_frames():
    p = P.RunProgress.from_text(_solver("matlab_cpu") + "[live] acoustic /out/live/acoustic_sub-001.png\n")
    assert p.live == {"acoustic": "/out/live/acoustic_sub-001.png"}
    done = P.RunProgress.from_text("⏱ prestus_pipeline 3.0s\nPipeline finished successfully\n")
    assert done.estimate(T0).fraction == 1.0
    failed = P.RunProgress.from_text(_solver("cpp_gpu") + "[PRESTUS Nexus] FAILED: out of memory\n")
    assert failed.failed and "Failed" in failed.estimate(T0).text


def test_no_progress_lines_means_nothing_to_show():
    p = P.RunProgress.from_text("Loading parameters\nsome output\n")
    assert not p.has_news() and not p.live


@pytest.mark.parametrize("text,seconds", [("12.5s", 12.5), ("1min 23.4s", 83.4),
                                          ("2hours 3min 4s", 7384), ("1 days, 2 hours, and 3 min", 93780),
                                          ("not a time", None)])
def test_scale_time(text, seconds):
    assert P.scale_time_seconds(text) == (pytest.approx(seconds) if seconds is not None else None)


def test_relevant_lines_keep_what_the_parser_reads():
    for kind, (name, body, now, *_rest) in SOLVER_LOGS.items():
        stage = "thermal" if "thermal" in kind else "acoustic"
        text = _solver(name, stage) + "noise line\n" + body + "[live] acoustic /a.png\n"
        kept = P.relevant_lines(text)
        assert "noise line" not in kept
        assert P.RunProgress.from_text(kept).estimate(now) == P.RunProgress.from_text(text).estimate(now)


# ── job logs ────────────────────────────────────────────────────────────

def _job_with_log(tmp_path, body, state="RUNNING"):
    (tmp_path / "sub-001_slurm_output_77.log").write_text(body, encoding="utf-8")
    return jobs.Job("77", "PS", state=state, work_dir=str(tmp_path))


def test_read_progress_local(tmp_path):
    name, body, now, frac, *_ = SOLVER_LOGS["kwave_cpp"]
    job = _job_with_log(tmp_path, "noise\n" * 1000 + _solver(name) + body)
    jobs.read_progress([job])
    assert "noise" not in job.progress_log
    assert jobs.job_progress(job).estimate(now).fraction == pytest.approx(frac)


def test_read_progress_skips_pending_and_reads_a_finished_job_once(tmp_path):
    job = _job_with_log(tmp_path, _solver("cpp_cpu"), state="PENDING")
    jobs.read_progress([job])
    assert job.progress_log == ""
    job.state = "COMPLETED"
    jobs.read_progress([job])
    assert job.progress_log.endswith("#end")
    (tmp_path / "sub-001_slurm_output_77.log").write_text("changed later", encoding="utf-8")
    before = job.progress_log
    jobs.read_progress([job])
    assert job.progress_log == before
    assert jobs.job_progress(job).finished


def test_failed_job_progress_is_failed(tmp_path):
    job = _job_with_log(tmp_path, _solver("cpp_cpu"), state="FAILED")
    jobs.read_progress([job])
    assert jobs.job_progress(job).failed


def test_progress_over_ssh_one_call_for_all_jobs(monkeypatch):
    calls = []

    class Remote:
        def to_local_or_none(self, path):
            return None

        def ssh(self, script):
            calls.append(script)
            return ["ssh", "host", script]
    out = (f"{jobs.PROGRESS_MARK}1\n" + _solver("cpp_gpu") + "|  10%  |  1.0s  |  9.0s  | x |\n"
           f"{jobs.PROGRESS_MARK}2\n")
    monkeypatch.setattr(jobs, "_run", lambda argv, timeout=15: out)
    js = [jobs.Job("1", state="RUNNING", work_dir="/project/a"),
          jobs.Job("2", state="RUNNING", work_dir="/project/b dir")]
    jobs.read_progress(js, Remote())
    assert len(calls) == 1 and "'/project/b dir'" in calls[0] and "grep -a -E" in calls[0]
    assert jobs.job_progress(js[0]).estimate(T0).fraction == pytest.approx(0.10)
    assert js[1].progress_log == ""


def test_poll_reads_progress_for_slurm_only(monkeypatch):
    seen = []
    monkeypatch.setattr(jobs, "query", lambda js, kind=None, remote=None: js)
    monkeypatch.setattr(jobs, "read_progress", lambda js, remote=None: seen.append(len(js)))
    jobs.poll([jobs.Job("1")], "slurm")
    jobs.poll([jobs.Job("1")], "qsub")
    assert seen == [1]


# ── display ─────────────────────────────────────────────────────────────

@pytest.fixture
def page(schema, prestus_root, tmp_path):
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])
    m = ConfigModel(schema)
    m.save(tmp_path / "config_p.yaml")
    window = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)
    p = window._page(RUN)
    p.app = app
    yield p
    p.jobs_timer.stop()
    p.progress_timer.stop()


def test_run_page_shows_local_progress_and_live_frame(page, tmp_path):
    from PyQt6.QtGui import QColor, QImage
    frame = tmp_path / "acoustic_sub-001.png"
    img = QImage(40, 20, QImage.Format.Format_RGB32)
    img.fill(QColor("red"))
    img.save(str(frame))
    now = time.time()
    page.task = "Run"
    page.progress = P.RunProgress()
    page.progress.feed(_solver("cpp_gpu", t=now) + "|    30%   |   3.0s  |   7.0s  | x |\n"
                       f"[live] acoustic {frame}\n", now)
    page._show_progress()
    panel = page.progress_panel
    assert not panel.isHidden() and panel.bar.value() == 300
    assert panel.title.text() == "Acoustic simulation" and panel.percent.text() == "30 %"
    assert "k-Wave C++, CUDA" in panel.detail.text() and not panel.b_live.isHidden()
    side = page.window.run_side                         # the Run section's sub-tabs
    live = page.live_page
    i = side.indexOf(live)
    assert i >= 0 and side.isTabVisible(i) and side.currentWidget() is live   # switched to, once
    assert live.view.tabs.count() == 1 and live.view.tabs.tabText(0) == "Pressure · live"
    assert live.view.image.pixmap() is not None and not live.view.image.pixmap().isNull()
    assert live.card.title.text() == "Acoustic simulation" and live.card.b_live.isHidden()
    side.setCurrentWidget(page)                         # the user goes back to Launch
    page._show_progress()
    assert side.currentWidget() is page                 # not switched again in this run
    panel.b_live.click()                                # the card's button: to the Live view tab
    assert side.currentWidget() is live


def test_run_without_frames_shows_no_live_view(page):
    page.task = "Run"
    page.progress = P.RunProgress.from_text(_solver("cpp_gpu", t=time.time()))
    page._show_progress()
    side = page.window.run_side
    assert page.progress_panel.b_live.isHidden()
    assert not side.isTabVisible(side.indexOf(page.live_page))     # hidden from the start


def test_run_page_hides_panel_without_progress_and_for_submit(page):
    page.task = "Validate"
    page.progress = P.RunProgress.from_text("Loaded: subject 1\n")
    page._show_progress()
    assert page.progress_panel.isHidden()
    page.task = "Submit"
    page.progress = P.RunProgress.from_text(_solver("cpp_gpu"))
    page._show_progress()
    assert page.progress_panel.isHidden()


def test_jobs_view_progress_column_and_completion_signal(page, tmp_path, monkeypatch):
    from prestus_nexus.ui.jobs_view import PROGRESS_COL
    name, body, *_ = SOLVER_LOGS["thermal_steps"]
    now = time.time()
    job = jobs.Job("5", "PS", state="RUNNING", work_dir=str(tmp_path),
                   progress_log=P.relevant_lines(_solver(name, "thermal", t=now - 30) +
                                                 body.replace(str(T0 + 10), str(now - 20))
                                                     .replace(str(T0 + 20), str(now - 10))))
    page.job_kind = "slurm"
    page.jobs = [job]
    view = page.jobs_view
    view.show_jobs()
    view.select("5")
    assert view.tree.topLevelItem(0).text(PROGRESS_COL).startswith("thermal 20 %")
    assert "peak 38.10 °C" in view.progress_panel.detail.text()
    assert not page.window.run_side.isTabVisible(page.window.run_side.indexOf(page.live_page))
    done = []
    page.job_completed.connect(done.append)
    finished = jobs.Job("5", "PS", state="COMPLETED", work_dir=str(tmp_path))
    page._jobs_polled([finished])
    assert done == ["5"]
    page._jobs_polled([finished])                  # once only
    assert done == ["5"]


def test_thermal_details_peak_so_far_and_current():
    assert P.describe_detail("t=2.00s Tmax=38.10C") == "peak 38.10 °C at t = 2.0 s"
    assert P.describe_detail("Tnow=38.10C") == "now 38.10 °C"


# ── BabelBrain sampler (functions/babelbrain/bb_progress.py) ─────────────

@pytest.fixture
def bbp(prestus_root, monkeypatch):
    monkeypatch.syspath_prepend(str(prestus_root / "functions" / "babelbrain"))
    import importlib
    mod = importlib.import_module("bb_progress")
    return importlib.reload(mod)


def _sample(bbp, monkeypatch, capsys, stage, frames, sampler=None):
    monkeypatch.setattr(bbp, "_frames", lambda: frames)
    sampler = sampler or bbp._Sampler(stage, "")
    sampler.sample()
    return sampler, capsys.readouterr().out


def test_bb_acoustic_sampler_prints_steps_per_fdtd_run(bbp, monkeypatch, capsys):
    args1, args2 = {"TimeSteps": 1000}, {"TimeSteps": 400}
    s, out = _sample(bbp, monkeypatch, capsys, "acoustic", {"_Execution": {"nStep": 249, "arguments": args1}})
    assert out.startswith("[progress] acoustic 250/1000 @") and out.rstrip().endswith("FDTD run 1")
    s, out = _sample(bbp, monkeypatch, capsys, "acoustic", {"_Execution": {"nStep": 249, "arguments": args1}}, s)
    assert out == ""                                     # unchanged: not printed again
    s, out = _sample(bbp, monkeypatch, capsys, "acoustic", {"_Execution": {"nStep": 9, "arguments": args2}}, s)
    assert "10/400" in out and "FDTD run 2" in out
    _, out = _sample(bbp, monkeypatch, capsys, "acoustic", {"other": {}})
    assert out == ""                                     # not in the solver loop: nothing


def test_bb_thermal_sampler_counts_on_and_off_phases(bbp, monkeypatch, capsys):
    p0, field = object(), object()
    cyc = {"nCurrent": 1, "TotalIterations": 3, "TotalDurationSteps": 100, "TotalDurationStepsOff": 50, "p0": p0}
    on = {"n": 9, "TotalDurationSteps": 100, "d_T0": None, "Pressure": field}
    _, out = _sample(bbp, monkeypatch, capsys, "thermal", {"BHTE": on, "RunBHTECycles": cyc})
    assert "[progress] thermal 160/450 @" in out         # 1 x 150 + 10
    off = {"n": 4, "TotalDurationSteps": 50, "d_T0": None, "Pressure": p0}
    _, out = _sample(bbp, monkeypatch, capsys, "thermal", {"BHTE": off, "RunBHTECycles": cyc})
    assert "[progress] thermal 255/450 @" in out         # 1 x 150 + 100 + 5
    multi = {"n": 0, "TotalDurationSteps": 100, "d_T0": None}
    _, out = _sample(bbp, monkeypatch, capsys, "thermal", {"BHTEMultiplePressureFields": multi})
    assert "[progress] thermal 1/100 @" in out           # without the cycle loop: this call
    _, out = _sample(bbp, monkeypatch, capsys, "thermal", {"BHTE": {"n": 3}})
    assert out == ""                                     # before the time loop (no d_T0)


def test_bb_sampler_lines_parse_in_the_nexus(bbp, monkeypatch, capsys):
    s, out1 = _sample(bbp, monkeypatch, capsys, "acoustic", {"_Execution": {"nStep": 99, "arguments": {"TimeSteps": 1000}}})
    p = P.RunProgress.from_text(_solver("babelbrain") + out1)
    assert p.estimate(time.time()).fraction == pytest.approx(0.1)


def test_bb_frames_need_cupy(bbp):
    np = pytest.importorskip("numpy")
    a = np.ones((4, 4, 4))
    assert bbp._cupy_max(a) is None and bbp._cupy_plane(a, lambda x: x) is None


def test_bb_thermal_children_get_a_sampler(bbp, monkeypatch):
    import sys
    import types
    started = []
    monkeypatch.setattr(bbp, "start", lambda stage: started.append(stage))
    mods = {}
    for name in ("BabelBrain", "BabelBrain.Babel_Thermal", "BabelBrain.Babel_Thermal.CalculateThermalProcess",
                 "ThermalModeling", "ThermalModeling.CalculateTemperatureEffects"):
        mods[name] = types.ModuleType(name)
        monkeypatch.setitem(sys.modules, name, mods[name])
    mods["BabelBrain.Babel_Thermal.CalculateThermalProcess"].SubProcess = lambda *a: "sub"
    mods["ThermalModeling.CalculateTemperatureEffects"].RunInProcess = lambda *a: "chunk"
    bbp.patch_thermal_children()
    bbp.patch_thermal_children()                         # idempotent
    assert mods["BabelBrain.Babel_Thermal.CalculateThermalProcess"].SubProcess(1) == "sub"
    assert mods["ThermalModeling.CalculateTemperatureEffects"].RunInProcess() == "chunk"
    assert started == ["thermal", "thermal"]


# ── jwave runner helpers (functions/jwave/run_jwave_sim.py) ─────────────
# The script imports jax and h5py at the top; its progress and frame helpers need only
# numpy, so they are taken from its source.

@pytest.fixture
def jw(prestus_root, tmp_path):
    import ast
    import types
    np = pytest.importorskip("numpy")
    src = (prestus_root / "functions" / "jwave" / "run_jwave_sim.py").read_text(encoding="utf-8")
    wanted = {"_progress", "_write_data", "_live_frame", "_progress_last", "_live_state",
              "LIVE_INTERVAL_S"}
    nodes = [n for n in ast.parse(src).body
             if (isinstance(n, ast.FunctionDef) and n.name in wanted)
             or (isinstance(n, ast.Assign) and any(getattr(t, "id", "") in wanted for t in n.targets))]
    mod = types.ModuleType("jw_helpers")
    mod.np, mod.os = np, os
    mod.LIVE_FILE = str(tmp_path / "live" / "acoustic_sub-001.json")
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "run_jwave_sim.py", "exec"), mod.__dict__)
    return mod


def test_jwave_progress_lines_rate_limited_and_parsed(jw, capsys):
    jw._progress(1, 1000, "jwave")
    for k in range(2, 1000):
        jw._progress(k, 1000, "jwave")             # fast loop: limited
    jw._progress(1000, 1000, "jwave")
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("[progress] acoustic 1/1000 @") and lines[-1].startswith("[progress] acoustic 1000/1000 @")
    assert len(lines) < 10
    p = P.RunProgress.from_text(_solver("jwave_pseudo_spectral") + "\n".join(lines) + "\n")
    assert p.estimate(time.time()).fraction == 1.0


def test_jwave_live_frame_is_a_data_frame_announced_once(jw, capsys):
    from prestus_nexus.ui.live_field import load_frame
    np = jw.np
    field = np.outer(np.hanning(20), np.hanning(30)) * 2e6          # |p| [Pa]
    jw._live_state["spacing_mm"] = [0.5, 0.5]
    jw._live_frame(field, force=True)
    jw._live_frame(field * 2, force=True)
    out = capsys.readouterr().out
    assert out.count("[live] acoustic ") == 1
    meta, plane = load_frame(jw.LIVE_FILE)
    assert meta["units"] == "MPa" and meta["spacing_mm"] == [0.5, 0.5] and plane.shape == (20, 30)
    assert abs(plane.max() - 2 * field.max() / 1e6) < 1e-3           # the second frame, in MPa
    assert not list(Path(jw.LIVE_FILE).parent.glob("*.part"))
    jw._live_frame(field)                              # within 5 s and not forced: skipped
    assert capsys.readouterr().out == ""


def test_jobs_view_live_frames_in_the_live_view_tab(page, tmp_path):
    from PyQt6.QtGui import QColor, QImage
    frame = tmp_path / "thermal_sub-001.png"
    img = QImage(30, 30, QImage.Format.Format_RGB32)
    img.fill(QColor("orange"))
    img.save(str(frame))
    now = time.time()
    log = P.relevant_lines(_solver("babelbrain", "thermal", t=now - 30) +
                           f"[progress] thermal 10/100 @{now - 5} Tnow=37.90C\n[live] thermal {frame}\n")
    page.job_kind = "slurm"
    page.jobs = [jobs.Job("8", "PS", state="RUNNING", work_dir=str(tmp_path), progress_log=log),
                 jobs.Job("9", "PS", state="RUNNING", work_dir=str(tmp_path))]
    view, side, live = page.jobs_view, page.window.run_side, page.live_page
    page.show_jobs_tab("8")
    view.show_jobs()
    view.select("8")
    assert side.isTabVisible(side.indexOf(live)) and side.currentWidget() is view   # shown, not switched to
    assert live.source == "8" and "Job 8" in live.view.heading.text()
    assert live.view.tabs.tabText(0) == "Temperature · live"
    view.select("9")                                    # no frames: the tab keeps job 8
    assert live.source == "8"
    view.select("8")
    view.progress_panel.b_live.click()                  # the card's button switches to it
    assert side.currentWidget() is live


def test_live_frame_margin_takes_the_gui_colour_but_not_the_plot(page):
    from PyQt6.QtCore import QRect
    from PyQt6.QtGui import QColor, QImage, QPainter

    from prestus_nexus.ui.progress_panel import blend_margin
    img = QImage(60, 40, QImage.Format.Format_RGB32)
    img.fill(QColor("white"))                          # figure background (margin)
    p = QPainter(img)
    p.setPen(QColor("black"))
    p.drawRect(QRect(10, 10, 30, 20))                  # axes box; inside stays white (field zero)
    p.end()
    gui = QColor("#e9edf1")
    out = blend_margin(img, gui)
    assert out.pixelColor(0, 0).name() == gui.name() and out.pixelColor(59, 39).name() == gui.name()
    assert out.pixelColor(25, 20).name() == "#ffffff"   # enclosed by the axes: unchanged
    odd = QImage(img)
    odd.setPixelColor(0, 0, QColor("red"))              # corners differ: no uniform margin
    assert blend_margin(odd, gui).pixelColor(59, 0).name() == "#ffffff"



# ── C++ phases, child jobs, series counters, temperature curve, whole run ─

def test_kwave_cpp_phases_before_and_after_the_time_loop():
    p = P.RunProgress.from_text(_solver("cpp_gpu") + "| FFT plans creation:                 Done |\n")
    assert "creating FFT plans" in p.estimate(T0 + 5).detail and p.estimate(T0 + 5).fraction is None
    p.feed("|    50%   |   1.0s  |   1.0s  | x |\n")
    assert "creating FFT plans" not in p.estimate(T0 + 5).detail     # the table replaces it
    p.feed("|   100%   |   2.0s  |   0.0s  | x |\n+ Creating output file:   Done\n")
    assert "writing the output" in p.estimate(T0 + 5).detail
    assert "Creating output file" in P.relevant_lines("+ Creating output file:   Done\nnoise")


def test_children_outer_counter_and_temperature_series():
    text = ('Job "PS_depth_2" (ID: 77) submitted successfully\n'
            f"[progress] calibration 2 of 5 @{T0} x\n".replace(" of ", "/") +
            _solver("kWaveDiffusion", "thermal") +
            f"[progress] thermal 5/50 @{T0 + 5} t=1.00s Tmax=37.40C\n"
            f"[progress] thermal 10/50 @{T0 + 10} t=2.00s Tmax=37.90C\n")
    p = P.RunProgress.from_text(text)
    assert p.children == [("PS_depth_2", "77")]
    assert p.estimate(T0 + 10).outer.startswith("Depth 2 of 5")
    assert p.temperature == [(1.0, 37.4, True, "Tmax"), (2.0, 37.9, True, "Tmax")]
    assert p.solver.stage == "thermal"                 # the counter is not a solver
    assert 'submitted successfully' in P.relevant_lines(text)


def test_whole_run_progress_from_predicted_stage_times():
    expected = {"preproc": 100.0, "acoustic": 300.0, "thermal": 600.0}
    text = ("▶ prestus_pipeline [#1] x\n▶ preproc [#1] x\n⏱ preproc 90.0s\n▶ acoustic [#1] x\n"
            f"[progress] acoustic 0/1 @{T0} solver cpp_gpu\n|    50%   |   10.0s  |   20.0s  | x |\n")
    p = P.RunProgress.from_text(text)
    frac, left = p.overall(expected, T0)
    assert frac == pytest.approx((100 + 150) / 1000) and left == pytest.approx(20 + 600)
    assert p.overall(None, T0) is None and p.overall({}, T0) is None


def test_timeout_risk_child_jobs_and_store_round_trip():
    running = jobs.Job("1", "PS", "RUNNING", time_left="0:10:00")
    assert jobs.timeout_risk(running, 3600) and not jobs.timeout_risk(running, 60)
    assert not jobs.timeout_risk(jobs.Job("2", state="PENDING", time_left="0:10:00"), 3600)
    parent = jobs.Job("5", "drive", "RUNNING", progress_log='Job "PS_d1" (ID: 6) submitted successfully')
    kids = jobs.child_jobs([parent])
    assert [(k.job_id, k.name, k.parent) for k in kids] == [("6", "PS_d1", "5")]
    assert jobs.child_jobs([parent] + kids) == []                    # once
    parent.expected = {"acoustic": 10.0}
    data = jobs.to_store([parent] + kids, "slurm", True)
    assert all("progress_log" not in row for row in data["jobs"])
    back, kind, remote = jobs.from_store(data)
    assert kind == "slurm" and remote and [j.job_id for j in back] == ["5", "6"]
    assert back[0].expected == {"acoustic": 10.0} and back[1].parent == "5"
    assert jobs.from_store("garbage") == ([], None, False)


def test_jobs_are_restored_in_a_new_window_and_cleared(schema, prestus_root, tmp_path):
    from PyQt6.QtWidgets import QApplication

    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import RUN, MainWindow
    app = QApplication.instance() or QApplication([])  # noqa: F841
    sess = Session(tmp_path / "s.yaml")
    sess.set(jobs.STORE_KEY, jobs.to_store([jobs.Job("31", "PS", "COMPLETED"),
                                            jobs.Job("32", "PS", "COMPLETED")], "slurm", False))
    sess.save()
    m = ConfigModel(schema)
    m.save(tmp_path / "config_r.yaml")
    page = MainWindow(m, Session(tmp_path / "s.yaml"), prestus_root)._page(RUN)
    assert [j.job_id for j in page.jobs] == ["31", "32"] and page.job_kind == "slurm"
    page.clear_finished_jobs()
    assert page.jobs == [] and jobs.from_store(page.window.session.get(jobs.STORE_KEY))[0] == []
    page.jobs_timer.stop()


def test_ended_jobs_notify_and_children_are_followed(page, monkeypatch):
    from prestus_nexus.ui import run_page as rp
    told = []
    monkeypatch.setattr(rp, "notify", lambda w, title, text: told.append((title, text)))
    monkeypatch.setattr(page, "refresh_jobs", lambda: None)
    page.job_kind = "slurm"
    page.jobs = [jobs.Job("40", "drive", "RUNNING")]
    log = 'Job "PS_d1" (ID: 41) submitted successfully'
    page._jobs_polled([jobs.Job("40", "drive", "RUNNING", progress_log=log)])
    assert [j.job_id for j in page.jobs] == ["40", "41"] and page.jobs[1].parent == "40"
    assert told == []
    page._jobs_polled([jobs.Job("40", "drive", "FAILED", progress_log=log)])
    assert told and "1 cluster job ended (1 not completed)" in told[0][0]


def test_card_whole_run_line_and_timeout_warning(page):
    from prestus_nexus.ui.jobs_view import timeout_text
    card = page.progress_panel
    p = P.RunProgress.from_text("▶ prestus_pipeline [#1] x\n▶ acoustic [#1] x\n"
                                f"[progress] acoustic 0/1 @{T0} solver cpp_gpu\n|    50%   |   10.0s  |   20.0s  | x |\n")
    card.show_progress(p, T0, expected={"acoustic": 40.0, "thermal": 60.0})
    assert not card.whole.isHidden() and "Whole run (predicted): ~20 %" in card.whole.text()
    card.show_progress(p, T0)
    assert card.whole.isHidden()
    job = jobs.Job("9", "PS", "RUNNING", time_left="0:00:30")
    assert "TIMEOUT" in timeout_text(job, 80.0) and timeout_text(job, 10.0) == ""


@pytest.mark.parametrize("order", ["F", "C"])
def test_data_frame_drawn_in_the_live_view(page, tmp_path, order):
    import json
    np = pytest.importorskip("numpy")
    plane = np.arange(12, dtype="<f4").reshape(3, 4) + 37
    raw = tmp_path / "thermal_sub-001.f32"
    (plane.ravel(order=order)).astype("<f4").tofile(raw)
    meta = {"quantity": "temperature", "units": "°C", "label": "Temperature", "spacing_mm": [0.5, 1.0],
            "markers": {"trans": [1, 0], "focus": [1, 3]}, "data": raw.name, "shape": [3, 4], "order": order,
            "min": 37.0, "max": 48.0}
    js = tmp_path / "thermal_sub-001.json"
    js.write_text(json.dumps(meta), encoding="utf-8")
    from prestus_nexus.ui.live_field import load_frame
    got = load_frame(str(js))
    assert got is not None and (got[1] == plane).all()
    now = time.time()
    page.task = "Run"
    page.progress = P.RunProgress.from_text(
        _solver("kWaveDiffusion", "thermal", t=now) + f"[progress] thermal 5/50 @{now} t=1.00s Tmax=37.40C\n"
        f"[progress] thermal 10/50 @{now} t=2.00s Tmax=37.90C\n[live] thermal {js}\n")
    page._show_progress()
    view = page.live_page.view
    assert view.field is not None and view.stack.currentWidget() is view.field
    field_ax, curve_ax = view.field.figure.axes          # the field and the peak temperature curve
    assert len(field_ax.child_axes) == 1                 # its colour bar, attached to it
    view.field.figure.canvas.draw()
    assert curve_ax.get_position().x0 > field_ax.get_position().x1    # beside it: full height kept
    assert curve_ax.get_ylabel() == "peak [°C]"
    assert "Temperature" in view.caption.text()
