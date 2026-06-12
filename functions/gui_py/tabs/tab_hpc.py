"""Tab 7 – HPC."""

from __future__ import annotations
from widgets import TextEdit, FloatEdit, Dropdown, Checkbox
from tabs.base import BaseTab


class TabHPC(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        self._add_section("Scheduler")
        self._add_row("HPC profile",
                      Dropdown("hpc.name", ["default", "snellius"], "default", registry=r))
        self._add_row("Partition",
                      TextEdit("hpc.partition", "", registry=r))
        self._add_row("GPU request",
                      TextEdit("hpc.gpu", "", registry=r),
                      note="e.g.  nvidia_a100-sxm4-40gb:1")
        self._add_row("Reservation",
                      TextEdit("hpc.reservation", "", registry=r))
        self._add_row("Job prefix",
                      TextEdit("hpc.job_prefix", "PRESTUS", registry=r))

        self._add_section("Resources")
        self._add_row("Wall time limit",
                      TextEdit("hpc.timelimit", "04:00:00", registry=r))
        self._add_row("Memory limit",
                      FloatEdit("hpc.memorylimit", 20, registry=r),
                      note="GB")

        self._add_section("Behaviour")
        self._add_row("Wait for job",
                      Checkbox("hpc.wait_for_job", "Block until HPC job completes",
                               False, registry=r))
        self._add_row("Max wait checks",
                      FloatEdit("hpc.max_wait_checks", 540, registry=r),
                      note="At ~20 s per check; 540 ≈ 3 hours.")
        self._add_row("LD_LIBRARY_PATH",
                      TextEdit("hpc.ld_library_path", "", registry=r),
                      note='Set if SimNIBS shows "undefined symbol" errors on HPC.')

        self._finish()
