"""Running part of the pipeline from a page (a placement method, the segmentation): saved
first, then submitted as a job where the Run tab submits to a cluster (scheduler or SSH),
else run in MATLAB here; the log shows on the Run tab."""

from __future__ import annotations


def run_partial(window, task: str, name: str, what: str, overrides: list, extra: str = "",
                local_only: bool = False) -> str | None:
    """'job' if submitted, 'local' if started in MATLAB here, None if nothing started.
    local_only: it needs this computer (e.g. a display), so it never goes to the cluster."""
    from .. import launcher
    from .main_window import RUN
    window._changed()
    if not window.save():
        return None
    run = window.show_special(RUN)
    if not local_only and run.submits_to_cluster():
        return "job" if run.submit_partial(name, what, overrides, extra) else None
    script = launcher.write_partial_script(window.model.path, window.prestus_root, name, what,
                                           overrides, extra, local=True)
    return "local" if run.run_script(task, script) else None


class RunFollower:
    """Calls ``done(ok)`` once a run started with run_partial has finished: the local MATLAB
    run's end, or (for a job) every update of the Run tab's job list, so a page can re-read
    the files the job writes."""

    def __init__(self, window, task: str, done):
        self.window, self.task, self.done = window, task, done
        self.waiting = False       # a local run is going
        self.watching = False      # a job was submitted

    def start(self, how: str | None) -> None:
        from .main_window import RUN
        run = self.window.pages[RUN]
        if how == "local" and not self.waiting:
            self.waiting = True
            run.finished.connect(self._finished)
        elif how == "job" and not self.watching:
            self.watching = True
            run.jobs_changed.connect(self._jobs_changed)

    def _finished(self, task: str, ok: bool) -> None:
        from .main_window import RUN
        if task != self.task:
            return
        self.window.pages[RUN].finished.disconnect(self._finished)
        self.waiting = False
        self.done(ok)

    def _jobs_changed(self) -> None:
        self.done(None)
