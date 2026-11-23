"""Progress and live view of BabelBrain simulations for the PRESTUS Nexus.

BabelBrain prints no progress while its solvers run, and has no callback for it. Both
solvers (BabelViscoFDTD) step through time in a Python loop, though:

- acoustic: StaggeredFDTD_3D_With_Relaxation_<backend>._Execution, ``for nStep in
  range(TimeSteps)``; the running RMS (or peak) pressure accumulator is ArraysGPUOp['SqrAcc'];
- thermal: RayleighAndBHTE.BHTE / BHTEMultiplePressureFields, ``for n in
  range(TotalDurationSteps)``, inside CalculateTemperatureEffects.RunBHTECycles
  (``for nCurrent in range(..., TotalIterations)``: one on and one off call per
  repetition); the temperature is d_T0 / d_T1 (ping-pong).

A daemon thread started here reads those loop variables from the main thread's frames
every INTERVAL_S seconds (read-only, under the GIL) and prints

    [progress] <stage> <k>/<n> @<unix time> <detail>

the line PRESTUS' prestus_progress prints, which the Nexus parses. With
PRESTUS_LIVE_FILE set (simulation.live_view; babelbrain_wrapper and
babelbrain_thermal_wrapper set it to a .json path), it also saves the plane through the
field's maximum as a data frame there (float32 .f32 + JSON, as live_view.m 'data') and
prints ``[live] <stage> <path>`` once. Field frames need the
CUDA backend (CuPy arrays); with OpenCL, Metal, or MLX only the progress is reported.

Nothing here may break a simulation: every sample is wrapped in try/except.

BabelBrain runs the thermal solver in child processes (multiprocessing, fork on Linux):
patch_thermal_children wraps their targets so that each child starts its own sampler.
"""

from __future__ import annotations

import os
import sys
import threading
import time

INTERVAL_S = 5.0
LIVE_ENV = "PRESTUS_LIVE_FILE"

_started: set[tuple[int, str]] = set()        # (pid, stage): one sampler per process and stage


def start(stage: str) -> None:
    """Start the sampler for 'acoustic' or 'thermal' in this process (once)."""
    key = (os.getpid(), stage)
    if key in _started:
        return
    _started.add(key)
    sampler = _Sampler(stage, os.environ.get(LIVE_ENV, ""))
    threading.Thread(target=sampler.run, name=f"prestus-progress-{stage}", daemon=True).start()


def patch_thermal_children() -> None:
    """Start a thermal sampler in BabelBrain's thermal child processes: wrap the process
    targets CalculateThermalProcess.SubProcess and CalculateTemperatureEffects.RunInProcess
    (looked up when the child is created, so the wrappers are what the children run)."""
    for mod_name, fn_name in (("BabelBrain.Babel_Thermal.CalculateThermalProcess", "SubProcess"),
                              ("ThermalModeling.CalculateTemperatureEffects", "RunInProcess")):
        try:
            mod = __import__(mod_name, fromlist=[fn_name])
            orig = getattr(mod, fn_name)
            if getattr(orig, "_prestus_progress", False):
                continue
            setattr(mod, fn_name, _with_sampler(orig))
        except Exception as e:                  # older/newer BabelBrain: progress only where found
            print(f"[bb_progress] Thermal progress not hooked into {mod_name}.{fn_name}: {e}",
                  flush=True)


def _with_sampler(fn):
    def wrapped(*args, **kwargs):
        start("thermal")
        return fn(*args, **kwargs)
    wrapped._prestus_progress = True
    wrapped.__name__ = getattr(fn, "__name__", "wrapped")
    return wrapped


# ── sampling ────────────────────────────────────────────────────────────

def _frames():
    """Locals of the main thread's frames, innermost first: {function name: locals}."""
    frame = sys._current_frames().get(threading.main_thread().ident)
    found: dict[str, dict] = {}
    while frame is not None:
        name = frame.f_code.co_name
        if name not in found:
            found[name] = frame.f_locals
        frame = frame.f_back
    return found


class _Sampler:
    def __init__(self, stage: str, live_file: str):
        self.stage = stage
        self.live_file = live_file
        self.last = None             # last (k, n) printed
        self.run_ids: list[int] = []  # FDTD runs seen (forward, back-propagation, refocused)
        self.announced = False

    def run(self) -> None:
        while True:
            time.sleep(INTERVAL_S)
            try:
                self.sample()
            except Exception:
                pass

    def sample(self) -> None:
        loc = _frames()
        if self.stage == "acoustic":
            self._acoustic(loc.get("_Execution"), loc.get("StaggeredFDTD_3D_with_relaxation") or {})
        else:
            self._thermal(loc.get("BHTE") or loc.get("BHTEMultiplePressureFields"),
                          loc.get("RunBHTECycles"))

    def _print(self, k: int, n: int, detail: str) -> None:
        if (k, n) == self.last or n <= 0:
            return
        self.last = (k, n)
        print(f"[progress] {self.stage} {min(k, n)}/{n} @{time.time():.1f} {detail}".rstrip(), flush=True)

    # acoustic ---------------------------------------------------------
    def _acoustic(self, ex, model) -> None:
        if not ex or "nStep" not in ex or "arguments" not in ex:
            return
        args = ex["arguments"]
        if id(args) not in self.run_ids:
            self.run_ids.append(id(args))
        run = self.run_ids.index(id(args)) + 1
        self._print(int(ex["nStep"]) + 1, int(args["TimeSteps"]), f"FDTD run {run}")
        gpu = ex.get("ArraysGPUOp") or {}
        if self.live_file and "SqrAcc" in gpu:
            acc = gpu["SqrAcc"]
            plane = _cupy_plane(acc, lambda a: a[(...,) + (a.shape[3] - 1, 0)] if a.ndim == 5 else a)
            if plane is not None:                         # sum of p² so far (SelRMSorPeak = 1)
                import numpy as np
                plane = np.sqrt(np.clip(np.asarray(plane, dtype=float), 0, None))
                top = float(np.nanmax(plane)) if plane.size else 0.0
                self._save(plane / top if top > 0 else plane,
                           {"quantity": "pressure", "units": "relative",
                            "label": "RMS pressure so far (relative), plane through the maximum",
                            "spacing_mm": _spacing_mm(model.get("SpatialStep"))})

    # thermal ----------------------------------------------------------
    def _thermal(self, bh, cyc) -> None:
        if not bh or "TotalDurationSteps" not in bh or "d_T0" not in bh:
            return
        n_call, steps = int(bh.get("n", 0)), int(bh["TotalDurationSteps"])
        k, n = n_call + 1, steps
        if cyc and "TotalIterations" in cyc:
            on, off = int(cyc.get("TotalDurationSteps", steps)), int(cyc.get("TotalDurationStepsOff", 0))
            is_off = "p0" in cyc and bh.get("Pressure") is cyc["p0"]
            per = on + off
            k = int(cyc.get("nCurrent", 0)) * per + (on if is_off else 0) + n_call + 1
            n = int(cyc["TotalIterations"]) * per
        temp = bh.get("d_T1") if n_call % 2 == 0 else bh.get("d_T0")
        detail = ""
        tmax = _cupy_max(temp)
        if tmax is not None:
            detail = f"Tnow={tmax:.2f}C"
        self._print(k, n, detail)
        if self.live_file and temp is not None:
            plane = _cupy_plane(temp, lambda a: a)
            if plane is not None:
                self._save(plane, {"quantity": "temperature", "units": "°C",
                                   "label": "Temperature, plane through the maximum",
                                   "spacing_mm": _spacing_mm(bh.get("dx"))})

    # frames -----------------------------------------------------------
    def _save(self, plane, meta: dict) -> None:
        """A data frame as live_view.m writes it ('data'): float32 plane + JSON, JSON last."""
        try:
            import json

            import numpy as np
            a = np.ascontiguousarray(np.asarray(plane, dtype="<f4"))
            folder = os.path.dirname(self.live_file)
            if folder:
                os.makedirs(folder, exist_ok=True)
            raw = self.live_file[:-5] + ".f32"
            a.tofile(raw + ".part")
            os.replace(raw + ".part", raw)
            meta = dict(meta, data=os.path.basename(raw), shape=list(a.shape), order="C",
                        min=float(np.nanmin(a)), max=float(np.nanmax(a)), markers={})
            with open(self.live_file + ".part", "w", encoding="utf-8") as fh:
                json.dump(meta, fh)
            os.replace(self.live_file + ".part", self.live_file)
            if not self.announced:
                self.announced = True
                print(f"[live] {self.stage} {self.live_file}", flush=True)
        except Exception as e:
            if not self.announced:
                self.announced = True              # report once, then stay quiet
                print(f"[bb_progress] Live view frame not written: {e}", flush=True)


def _spacing_mm(step) -> list[float]:
    """BabelBrain's grid step [m] (a scalar) as [mm, mm]; [1, 1] when unknown."""
    try:
        mm = 1e3 * float(step)
        return [mm, mm] if mm > 0 else [1.0, 1.0]
    except (TypeError, ValueError):
        return [1.0, 1.0]


def _is_cupy(a) -> bool:
    return a is not None and type(a).__module__.startswith("cupy")


def _cupy_max(a) -> float | None:
    if not _is_cupy(a):
        return None
    import cupy as cp
    with cp.cuda.Device(a.device.id):
        return float(a.max())


def _cupy_plane(a, select):
    """The x-z plane (2nd axis fixed) through the maximum of select(a), copied to the host;
    None when a is not a CuPy array (OpenCL, Metal, MLX backends)."""
    if not _is_cupy(a):
        return None
    import cupy as cp
    import numpy as np
    with cp.cuda.Device(a.device.id):
        vol = select(a)
        if vol.ndim != 3:
            return None
        iy = int(np.unravel_index(int(cp.argmax(vol)), vol.shape)[1])
        return cp.asnumpy(vol[:, iy, :])
