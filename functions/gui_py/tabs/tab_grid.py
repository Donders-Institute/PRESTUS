"""Tab 4 – Grid."""

from __future__ import annotations
from widgets import FloatEdit, Dropdown, Checkbox, XYZWidget
from tabs.base import BaseTab


class TabGrid(BaseTab):
    def __init__(self, registry: dict, parent=None):
        super().__init__(registry, parent)
        r = registry

        self._add_section("Spatial Resolution")
        self._add_row("* Grid resolution",
                      FloatEdit("grid.resolution_mm", 0.5, registry=r),
                      note="mm")
        self._add_row("Minimum PPW",
                      FloatEdit("grid.min_ppw", 6, registry=r),
                      note="Minimum points per wavelength at transducer frequency.")

        self._add_section("Temporal Resolution")
        self._add_row("CFL number",
                      FloatEdit("grid.source_cfl", 0.15, registry=r),
                      note="Courant-Friedrichs-Lewy fraction. Default 0.15 for stability in heterogeneous skull.")
        self._add_row("PPW override",
                      FloatEdit("grid.source_ppw", None, registry=r),
                      note="Leave NaN to compute PPW from resolution and max sound speed.")
        self._add_row("Stability limit fraction",
                      FloatEdit("grid.source_limit_fraction", 0.9, registry=r),
                      note="If dt exceeds stability limit, rescale to this fraction. Set 0 to disable.")

        self._add_section("Grid Dimensions")
        self._add_row("Default dims (water/phantom)",
                      XYZWidget("grid.default_dims", (144, 144, 400), registry=r))
        self._add_row("PML size",
                      FloatEdit("grid.pml_size", 10, registry=r),
                      note="voxels")
        self._add_row("Max grid expand",
                      FloatEdit("grid.max_expand", 40, registry=r),
                      note="voxels")
        self._add_row("Axisymmetric mode",
                      Checkbox("grid.axisymmetric", "2D axisymmetric (kspaceFirstOrderAS)",
                               False, registry=r))
        self._add_row("Use kWaveArray",
                      Checkbox("grid.use_kWaveArray",
                               "Recommended for accurate transducer modelling",
                               True, registry=r))

        self._add_section("Grid Orientation & Thermal Grid")
        self._add_row("Grid orientation mode",
                      Dropdown("grid.mode", ["transducer_axis", "ras_plus"],
                                "transducer_axis", registry=r),
                      note="transducer_axis: rotate so focal axis aligns with z. ras_plus: keep RAS+ orientation.")
        self._add_row("Thermal resolution",
                      FloatEdit("grid.thermal_resolution_mm", None, registry=r),
                      note="mm – leave NaN to use acoustic grid resolution.")
        self._add_row("Thermal FOV",
                      FloatEdit("grid.thermal_fov_mm", None, registry=r),
                      note="mm – leave NaN to match acoustic simulation FOV.")

        self._finish()
