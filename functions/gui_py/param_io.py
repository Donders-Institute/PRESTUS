"""Parameter ↔ widget I/O helpers.

Mirrors the MATLAB GUI's apply_params_to_gui / collect_params logic.
All widget tags use dotted paths (e.g. ``transducer.annular.elem_n``).
XYZ triplet sub-fields are tagged with a ``_0`` / ``_1`` / ``_2`` suffix.
"""

from __future__ import annotations
import copy
import math
from typing import Any

from PyQt6.QtWidgets import QCheckBox, QComboBox, QLineEdit, QPlainTextEdit

from widgets import FloatEdit, TextEdit, Dropdown, Checkbox, CSVArrayEdit, BrowseEdit


# ── CSV array field tags (numeric arrays stored as comma-separated strings) ────
CSV_ARRAY_TAGS = {
    "transducer.annular.elem_id_mm",
    "transducer.annular.elem_od_mm",
    "transducer.annular.elem_phase_deg",
    "transducer.target_isppa_wcm2",
    "placement.plantus.focal_distance_list",
    "placement.plantus.flhm_list",
    "calibration.elem_phase_correction_deg",
    "calibration.focal_depths_wrt_exit_plane",
    "calibration.desired_intensities",
    "calibration.opt_limits",
    "layers.water",
    "layers.brain",
    "layers.skin",
    "layers.skull",
    "layers.skull_cortical",
    "layers.skull_trabecular",
}

# XYZ triplet base tags — sub-fields are base_tag_0 / _1 / _2
XYZ_TAGS = {
    "transducer.trans_pos",
    "transducer.focus_pos",
    "grid.default_dims",
    "placement.heuristic.mni_target_mm",
    "placement.plantus.mni_target_mm",
}

# Medium table tissues / properties
TISSUES = ["water", "brain", "skin", "skull", "skull_cortical", "skull_trabecular"]
MEDIUM_PROPS = [
    "sound_speed", "density", "alpha_coeff", "alpha_power",
    "thermal_conductivity", "specific_heat_capacity", "perfusion",
    "absorption_fraction",
]
MEDIUM_DEFAULTS = {
    "water":           [1500,  1000,  0.002, 2.0, 0.6,  4000, 0,   0.02],
    "brain":           [1546,  1046,  0.6,   1.3, 0.51, 3630, 559, 0.02],
    "skin":            [1537,  1116,  0.8,   1.3, 0.37, 3391, 106, 0.02],
    "skull":           [2800,  1912,  11.6,  1.3, 0.32, 1313, 0,   0.02],
    "skull_cortical":  [2800,  1912,  11.6,  1.3, 0.32, 1313, 0,   0.02],
    "skull_trabecular":[2300,  1178,  6.9,   1.3, 0.32, 1313, 0,   0.02],
}


# ── Flatten / unflatten nested dicts ─────────────────────────────────────────

def flatten(d: dict, prefix: str = "") -> dict:
    """Recursively flatten a nested dict to dotted-key flat dict."""
    out = {}
    for k, v in d.items():
        full_key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, full_key))
        else:
            out[full_key] = v
    return out


def unflatten(flat: dict) -> dict:
    """Reconstruct a nested dict from a flat dotted-key dict."""
    result = {}
    for key, val in flat.items():
        parts = key.split(".")
        d = result
        for part in parts[:-1]:
            d = d.setdefault(part, {})
        d[parts[-1]] = val
    return result


def set_nested(d: dict, path: str, val) -> dict:
    parts = path.split(".")
    cur = d
    for part in parts[:-1]:
        cur = cur.setdefault(part, {})
    cur[parts[-1]] = val
    return d


def get_nested(d: dict, path: str):
    try:
        cur = d
        for part in path.split("."):
            cur = cur[part]
        return cur
    except (KeyError, TypeError):
        return None


# ── Apply params dict → widgets ───────────────────────────────────────────────

def apply_params(params: dict, registry: dict, medium_table=None):
    """Populate all registered widgets from a (possibly nested) params dict."""
    flat = flatten(params)

    for tag, widget in registry.items():
        val = flat.get(tag)
        if val is None:
            continue
        _set_widget(widget, val)

    # XYZ triplets: source value is a list/tuple [x, y, z]
    for base_tag in XYZ_TAGS:
        val = get_nested(params, base_tag)
        if val is None:
            continue
        if not isinstance(val, (list, tuple)):
            continue
        for i, v in enumerate(val[:3]):
            sub_tag = f"{base_tag}_{i}"
            if sub_tag in registry:
                _set_widget(registry[sub_tag], v)

    # CSV array fields: convert list to comma-separated string
    for tag in CSV_ARRAY_TAGS:
        if tag not in registry:
            continue
        val = flat.get(tag)
        if val is None:
            continue
        widget = registry[tag]
        if isinstance(val, (list, tuple)):
            widget.set_value(val)
        elif isinstance(val, str):
            widget.set_value(val)

    # Medium properties table
    if medium_table is not None and "medium_properties" in params:
        _populate_medium_table(medium_table, params["medium_properties"])


def _set_widget(widget, val):
    if hasattr(widget, "set_value"):
        widget.set_value(val)
    elif isinstance(widget, QCheckBox):
        widget.setChecked(bool(val))
    elif isinstance(widget, QComboBox):
        idx = widget.findText(str(val))
        if idx >= 0:
            widget.setCurrentIndex(idx)
        else:
            widget.addItem(str(val))
            widget.setCurrentText(str(val))
    elif isinstance(widget, QPlainTextEdit):
        widget.setPlainText(str(val) if val is not None else "")
    elif isinstance(widget, QLineEdit):
        widget.setText(str(val))


# ── Collect params dict from widgets ─────────────────────────────────────────

def collect_params(registry: dict, last_params: dict | None = None,
                   medium_table=None) -> dict:
    """Read all registered widgets and return a nested params dict.

    Starts from last_params so un-visited lazy tabs still contribute their
    values.  Widget values take precedence.
    """
    import copy
    params = copy.deepcopy(last_params) if last_params else {}
    # Remove fields that load_parameters manages internally
    for f in ("startup", "state"):
        params.pop(f, None)

    for tag, widget in registry.items():
        # Skip XYZ sub-fields — handled separately below
        base = _xyz_base(tag)
        if base is not None:
            continue
        # Skip CSV array fields — handled separately
        if tag in CSV_ARRAY_TAGS:
            continue

        val = _get_widget_value(widget)
        if val is None:
            continue
        if isinstance(val, str) and val == "":
            continue
        if isinstance(val, float) and math.isnan(val):
            continue
        # Skip sentinel "not selected" serial
        if tag == "transducer.serial" and val == "(manual)":
            continue
        set_nested(params, tag, val)

    # XYZ triplets → list
    for base_tag in XYZ_TAGS:
        vals = []
        any_found = False
        for i in range(3):
            sub_tag = f"{base_tag}_{i}"
            if sub_tag in registry:
                any_found = True
                v = _get_widget_value(registry[sub_tag])
                vals.append(v if v is not None else float("nan"))
        if any_found:
            # Check axisymmetric
            ax_widget = registry.get("grid.axisymmetric")
            if ax_widget and _get_widget_value(ax_widget) and base_tag in (
                "transducer.trans_pos", "transducer.focus_pos", "grid.default_dims"
            ):
                vals = [vals[0], vals[2]]
            set_nested(params, base_tag, vals)

    # CSV array fields → numeric lists
    for tag in CSV_ARRAY_TAGS:
        if tag not in registry:
            continue
        widget = registry[tag]
        val = widget.value() if hasattr(widget, "value") else None
        if val is not None:
            set_nested(params, tag, val)

    # Medium table
    if medium_table is not None:
        mp = _read_medium_table(medium_table)
        if mp:
            params["medium_properties"] = mp

    # Derived: dir_output
    sim_path = get_nested(params, "path.sim")
    subject_id = get_nested(params, "subject_id")
    if sim_path:
        sub = int(subject_id) if subject_id else 1
        import os
        params.setdefault("io", {})["dir_output"] = os.path.join(
            str(sim_path), f"sub-{sub:03d}"
        )

    return params


def _xyz_base(tag: str):
    """Return base_tag if tag is a XYZ sub-field (ends in _0/_1/_2), else None."""
    for base in XYZ_TAGS:
        for i in range(3):
            if tag == f"{base}_{i}":
                return base
    return None


def _get_widget_value(widget):
    if hasattr(widget, "value"):
        return widget.value()
    if isinstance(widget, QCheckBox):
        return widget.isChecked()
    if isinstance(widget, QComboBox):
        return widget.currentText()
    if isinstance(widget, QPlainTextEdit):
        return widget.toPlainText()
    if isinstance(widget, QLineEdit):
        return widget.text()
    return None


# ── Medium table helpers ──────────────────────────────────────────────────────

def _populate_medium_table(table, mp: dict):
    from PyQt6.QtWidgets import QTableWidgetItem
    for ti, tissue in enumerate(TISSUES):
        tp = mp.get(tissue, {})
        for pi, prop in enumerate(MEDIUM_PROPS):
            val = tp.get(prop)
            if val is not None:
                item = table.item(ti, pi)
                if item:
                    item.setText(str(val))


def _read_medium_table(table) -> dict:
    mp = {}
    from PyQt6.QtWidgets import QTableWidgetItem
    for ti, tissue in enumerate(TISSUES):
        mp[tissue] = {}
        for pi, prop in enumerate(MEDIUM_PROPS):
            item = table.item(ti, pi)
            if item:
                try:
                    mp[tissue][prop] = float(item.text())
                except ValueError:
                    pass
    return mp
