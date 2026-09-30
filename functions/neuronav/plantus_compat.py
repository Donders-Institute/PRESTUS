"""
plantus_compat.py  –  PRESTUS-side fixes applied to an imported PlanTUS module
=============================================================================
Shared by prestus_plantus_launcher.py and prestus_plantus_gui_launcher.py so
both modes behave identically and the upstream PlanTUS checkout can stay
unmodified.
"""

import numpy as np

# Compatibility shims for NumPy 2.0 removals used by upstream PlanTUS
# (applied on import, before either launcher imports PlanTUS).
if not hasattr(np, 'asfarray'):
    np.asfarray = lambda a, dtype=float: np.asarray(a, dtype=dtype)
if not hasattr(np, 'float'):
    np.float = float
if not hasattr(np, 'string_'):
    np.string_ = np.bytes_


def patch_compute_flhm(PlanTUS):
    """Make PlanTUS.compute_FLHM_for_focal_distance work with < 4 points.

    PlanTUS fits a cubic (4 coefficients) with curve_fit, which fails when
    fewer than 4 calibration points are given, e.g. the single point from
    the O'Neill geometric fallback in plantus_calib_from_calibration.
    Non-finite pairs (FLHM undeterminable for that focal depth) are dropped;
    below 4 points, FLHM is interpolated linearly (held constant outside the
    calibrated range).
    """
    orig = PlanTUS.compute_FLHM_for_focal_distance

    def compute_FLHM_for_focal_distance(focal_distance, focal_distance_list, flhm_list):
        fd = np.asarray(focal_distance_list, dtype=float).ravel()
        fl = np.asarray(flhm_list, dtype=float).ravel()
        valid = np.isfinite(fd) & np.isfinite(fl)
        fd, fl = fd[valid], fl[valid]
        if fd.size == 0:
            raise ValueError("No finite (focal_distance, FLHM) pairs in "
                             "focal_distance_list / flhm_list.")
        if fd.size >= 4:
            return orig(focal_distance, list(fd), list(fl))
        order = np.argsort(fd)
        return float(np.interp(focal_distance, fd[order], fl[order]))

    PlanTUS.compute_FLHM_for_focal_distance = compute_FLHM_for_focal_distance
