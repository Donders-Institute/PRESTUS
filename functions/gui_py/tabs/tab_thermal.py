"""Tab 6 – Thermal."""

from __future__ import annotations
from widgets import FloatEdit, Checkbox
from tabs.base import BaseTab
import param_io


class TabThermal(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        self._add_section("Sonication Protocol Timing")
        self._timing_widgets = []
        timing_fields = [
            ("timing.pd",              "Pulse duration (PD)",              None, "s"),
            ("timing.pri",             "Pulse repetition interval (PRI)",  None, "s"),
            ("timing.ptd",             "Pulse train duration (PTD)",       None, "s"),
            ("timing.ptri",            "Pulse train rep. interval (PTRI)", None, "s"),
            ("timing.ptrd",            "Pulse train rep. duration (PTRD)", None, "s"),
            ("timing.post_ptri_dur",   "Post-PTRI steady-state",           None, "s"),
            ("timing.pt_timestep",     "In-train timestep",               0.02, "s"),
            ("timing.post_pt_timestep","Post-train timestep",             1.0,  "s"),
        ]
        for tag, label, default, unit in timing_fields:
            w = FloatEdit(tag, default, registry=r)
            self._timing_widgets.append(w)
            self._add_row(label, w, note=unit)

        eq_step = Checkbox("timing.equal_step_duration",
                           "Force equal durations for on and off cycles",
                           False, registry=r)
        self._timing_widgets.append(eq_step)
        self._add_row("Equal on/off step durations", eq_step)

        self._add_section("CEM43 Thermal Dose")
        self._add_row("ISO CEM43",
                      Checkbox("thermal.cem43_iso",
                               "Use ISO definition (R=0.5 above 43°C; zero below 39°C; ∞ above 57°C)",
                               False, registry=r),
                      note="Default uses k-Wave kWaveDiffusion formula (Sapareto-Dewey).")

        self._add_section("Sensor & Recording")
        self._add_row("Sensor half-window size",
                      FloatEdit("thermal.sensor_xy_halfsize", 100, registry=r),
                      note="voxels – reduce if memory is limited.")
        self._add_row("Record T at every step",
                      Checkbox("thermal.record_t_at_every_step",
                               "Record full sensor window at every timestep (memory intensive)",
                               False, registry=r))

        self._add_section("Initial Temperature (°C)")
        for tissue in param_io.TISSUES:
            self._add_row(tissue,
                          FloatEdit(f"thermal.temp_0.{tissue}", 37, registry=r),
                          note="°C")

        self._finish()

    def set_timing_enabled(self, enabled: bool):
        for w in self._timing_widgets:
            w.setEnabled(enabled)
