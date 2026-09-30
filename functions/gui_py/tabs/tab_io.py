"""Tab 1 – I/O & Paths."""

from __future__ import annotations
from PyQt6.QtWidgets import QHBoxLayout, QWidget
from widgets import BrowseEdit, TextEdit, Dropdown, Checkbox
from tabs.base import BaseTab


class TabIO(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        # ── Subject ──────────────────────────────────────────────────────────
        self._add_section("Subject")
        self._add_row("* Subject ID", TextEdit("subject_id", "1", mandatory=True, registry=r))

        # ── Data Paths ───────────────────────────────────────────────────────
        self._add_section("Data Paths")
        self._add_row("* Structural MRI path",
                      BrowseEdit("path.anat", "dir", registry=r, mandatory=True))
        self._add_row("* Simulation output path",
                      BrowseEdit("path.sim", "dir", registry=r, mandatory=True))
        self._add_row("* Segmentation path",
                      BrowseEdit("path.seg", "dir", registry=r, mandatory=True))
        self._add_row("Localite raw path",
                      BrowseEdit("path.localite_raw", "dir", registry=r))
        self._add_row("T1 filename pattern",
                      TextEdit("path.t1_pattern", "sub-%1$03d_T1w.nii*", registry=r),
                      note="Relative to path.anat. Use %1$03d for subject ID substitution.")
        self._add_row("T2 filename pattern",
                      TextEdit("path.t2_pattern", "", registry=r),
                      note="Optional. Leave empty to run SimNIBS charm with T1 only.")
        self._add_row("pseudoCT base path",
                      BrowseEdit("path.pct", "dir", registry=r),
                      note="Leave empty to write pseudoCT into the SimNIBS m2m folder.")

        # ── Environment ───────────────────────────────────────────────────────
        self._add_section("Environment")
        self._add_row("* SimNIBS bin path",
                      BrowseEdit("startup.simnibs_bin_path", "dir", registry=r, mandatory=True))
        self._add_row("MATLAB bin path",
                      BrowseEdit("startup.matlab_bin_path", "file", registry=r),
                      note="Required for pseudoCT generation on HPC.")
        self._add_row("FSL bin path",
                      BrowseEdit("startup.fsl_bin_path", "dir", registry=r))
        self._add_row("ANTs bin path",
                      BrowseEdit("startup.ants_bin_path", "dir", registry=r))
        self._add_row("Extra paths (addpath)",
                      TextEdit("startup.paths_to_add", "", registry=r),
                      note="Semicolon-separated absolute paths.")
        self._add_row("Extra subpaths (genpath)",
                      TextEdit("startup.subpaths_to_add", "", registry=r),
                      note="Semicolon-separated; each added recursively.")

        # ── Output Settings ───────────────────────────────────────────────────
        self._add_section("Output Settings")
        self._add_row("Output affix",
                      TextEdit("io.output_affix", "", registry=r))
        self._add_row("Overwrite files",
                      Dropdown("io.overwrite_files", ["always", "force", "never", "ask"],
                               "always", registry=r))
        self._add_row("Overwrite SimNIBS",
                      Checkbox("io.overwrite_simnibs", "Re-run SimNIBS segmentation",
                               False, registry=r))

        # ── Matrix Saving ─────────────────────────────────────────────────────
        self._add_section("Matrix Saving")
        self._add_note(
            "Global flag is the fallback; individual flags take precedence when set."
        )
        self._add_row("Global save flag",
                      Checkbox("io.save_matrices", "Global fallback: save all intermediate matrices",
                               False, registry=r))
        self._add_row("Save source matrices",
                      Checkbox("io.save_source_matrices", "kwave_source.mat (shared across variants)",
                               True, registry=r))
        self._add_row("Save acoustic matrices",
                      Checkbox("io.save_acoustic_matrices",
                               "sensor_data + kgrid + medium (large)", False, registry=r))
        self._add_row("Save thermal matrices",
                      Checkbox("io.save_thermal_matrices",
                               "heating_res.mat — temperature timeseries", True, registry=r))
        self._add_row("Save heating video",
                      Checkbox("io.save_heatingvideo", "Save MP4 of incremental heating",
                               False, registry=r))

        # ── Advanced I/O ──────────────────────────────────────────────────────
        self._add_section("Advanced I/O")
        self._add_row("External acoustic NIfTI",
                      BrowseEdit("io.external_acoustic_nifti", "file", registry=r),
                      note="When set, skips acoustic simulation. NIfTI must contain p-max [Pa] in T1 space.")
        self._add_row("Save complex pressure",
                      Checkbox("io.save_p_complex",
                               "Save magnitude + phase of steady-state pressure field",
                               False, registry=r))
        self._add_row("Save property maps",
                      Checkbox("io.save_property_maps",
                               "Save per-tissue acoustic/thermal property NIfTIs",
                               False, registry=r))
        self._add_row("Save MNI outputs",
                      Checkbox("io.save_MNI",
                               "Save simulation outputs in MNI space in addition to native T1w",
                               True, registry=r))
        self._finish()
