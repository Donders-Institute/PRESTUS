"""The Results tab's reporting overview: one run's parameters grouped as the ITRUSST
consensus asks TUS studies to report them, from the files PRESTUS already writes. Qt-free.

Every row names where its value comes from, because the overview mixes sources:

- ``run CSV``          the run's metrics CSV (acoustic_analysis.m, thermal_analysis.m): as simulated
- ``free-water CSV``   the free-water (``water``) run with the same output affix, if any
- ``targeting CSV``    ``<tag>_<medium><affix>_targeting.csv`` (run_targeting_analysis)
- ``placement record`` ``<tag>_placement<affix>.json`` (write_placement_positions.m)
- ``config``           the open config, as it is now (it may differ from what the run used)
- ``library``          the equipment entry named by transducer.serial
- ``T1w intensity map (Nexus)`` focal widths measured by the Nexus on the run's T1w
  intensity map (focal_widths.py), for runs whose CSV has no focal_fwhm_* columns

Values are shown as written, at 4 significant digits; unit conversions (Pa to MPa, Hz to
kHz, s to ms) are the only arithmetic, apart from the pulse timing derived from the config
(PRF = 1/PRI, duty cycles) and the per-tissue ISPTA, which say so. What no output holds
reads 'not available'.

ISPTA follows acoustic_analysis.m: Ita = Ipa × PD/PRI, the temporal average over the pulse
repetition within a pulse train (the pulse-train duty cycle PTD/PTRI is not included). The
per-tissue ISPTA is the tissue's ISPPA times that duty cycle, which equals the tissue maximum
of Ita; the duty cycle is the run's own (Ispta/Isppa of its CSV), else PD/PRI of the config.

The non-significant-risk check uses the limits and the colour rule of PRESTUS's HTML report
(get_risk_limits.m, risk_color.m, build_safety_banner in generate_simulation_report.m):
within below 50 % of a limit, near from 50 % to 100 %, above beyond it.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from . import focal_widths as fwm
from . import results as rs
from .model import ConfigModel, is_empty
from .ui.errors import log_current

NA = "not available"
PENDING = "measuring…"
ISPTA_WINDOW = "pulse-train average, PD/PRI"

# get_risk_limits.m: finite ITRUSST NSR limits (informational Inf ones are left out)
LIMITS_LAYERED: dict[str, tuple[str, float, str]] = {
    "MI_tc": ("MI (transcranial)", 1.9, ""),
    "MI": ("MI (free water)", 1.9, ""),
    "MI_brain": ("MI (brain)", 1.9, ""),
    "MI_skull": ("MI (skull)", 1.9, ""),
    "MI_skin": ("MI (skin)", 1.9, ""),
    "riseT_brain": ("Temp rise (brain)", 2.0, "°C"),
    "riseT_skull": ("Temp rise (skull)", 2.0, "°C"),
    "riseT_skin": ("Temp rise (skin)", 2.0, "°C"),
    "riseT37_brain": ("ΔT from 37 °C (brain)", 2.0, "°C"),
    "riseT37_skull": ("ΔT from 37 °C (skull)", 2.0, "°C"),
    "riseT37_skin": ("ΔT from 37 °C (skin)", 2.0, "°C"),
    "CEM43_brain": ("CEM43 (brain)", 2.0, "min"),
    "CEM43_skull": ("CEM43 (skull)", 16.0, "min"),
    "CEM43_skin": ("CEM43 (skin)", 21.0, "min"),
    "maxT_brain": ("Max temp (brain)", 39.0, "°C"),
    "maxT_skull": ("Max temp (skull)", 39.0, "°C"),
    "maxT_skin": ("Max temp (skin)", 39.0, "°C"),
}
LIMITS_WATER: dict[str, tuple[str, float, str]] = {
    "MI": ("MI (free water)", 1.9, ""),
    "MI_tc": ("MI (transcranial)", 1.9, ""),
}
ISO_COLUMNS = {"CEM43_brain": "CEM43iso_brain", "CEM43_skull": "CEM43iso_skull",
               "CEM43_skin": "CEM43iso_skin"}
TISSUES = ("brain", "skull", "skin")


@dataclass
class Row:
    label: str
    value: str                       # with its unit, or NA
    source: str = ""
    status: str = ""                 # 'within' | 'near' | 'above' | '' (no limit / no value)
    limit: str = ""                  # e.g. '≤ 1.9'


@dataclass
class Section:
    title: str
    rows: list[Row] = field(default_factory=list)
    note: str = ""


@dataclass
class Report:
    run_label: str
    medium: str
    sections: list[Section]
    verdict: str                     # 'within' | 'review' | 'none'
    above: list[str] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)


# ── value helpers ──────────────────────────────────────────────────────────

def _num(v) -> float | None:
    if v is None or isinstance(v, bool):
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def fmt(x: float | None, unit: str = "", scale: float = 1.0) -> str:
    """x / scale at 4 significant digits with its unit; NA when there is no number."""
    if x is None:
        return NA
    v = x / scale
    s = f"{v:.4g}"
    return f"{s} {unit}".rstrip() if unit else s


def _vec(v) -> str:
    if isinstance(v, (list, tuple)) and v and all(_num(c) is not None for c in v):
        return "[" + ", ".join(f"{float(c):.4g}" for c in v) + "]"
    if isinstance(v, str) and v.strip():
        return v.strip()
    return NA


def risk_status(value: float | None, limit: float) -> str:
    """risk_color.m: within (< 50 %), near (50–100 %), above (> limit); '' without a value."""
    if value is None:
        return ""
    if value > limit:
        return "above"
    if value >= 0.5 * limit:
        return "near"
    return "within"


def _get(model: ConfigModel, key: str, ti: int = 0):
    try:
        return model.get(key, ti)
    except (KeyError, IndexError):  # a key or transducer this config does not have
        return None


def _tx(model: ConfigModel, key: str, ti: int):
    """transducer.<key> of transducer ti; else the library entry's geometry."""
    v = _get(model, f"transducer.{key}", ti)
    if not is_empty(v):
        return v
    dev = None
    try:
        dev = model.device(ti)
    except (KeyError, IndexError):  # no such transducer: no library geometry
        dev = None
    node = dev.geometry if dev is not None else None
    for part in key.split("."):
        node = node.get(part) if isinstance(node, dict) else None
    return node


def read_table(path: Path) -> dict[str, str]:
    """Last row of a CSV as {column: value} ({} if missing)."""
    if not path.is_file():
        return {}
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    return dict(zip(rows[0], rows[-1])) if len(rows) >= 2 else {}


def placement_record(run: rs.Run) -> dict:
    for affix in dict.fromkeys([run.affix, rs._TRANS.sub("", run.affix)]):
        p = run.base / f"{run.tag}_placement{affix}.json"
        try:
            data = json.loads(p.read_text(encoding="utf-8"))
        except FileNotFoundError:  # best-effort: no record under this affix: try the next
            continue
        except (OSError, ValueError):
            log_current(f"ITRUSST report: cannot read {p}")
            continue
        if isinstance(data, dict):
            return data
    return {}


def free_water_run(runs: list[rs.Run], run: rs.Run) -> rs.Run | None:
    """The free-water run with the same output affix (the free-water baseline or a water run)."""
    if run.medium == "water":
        return run
    return next((r for r in runs if r.medium == "water" and r.affix == run.affix and r.csv), None)


# ── sections ───────────────────────────────────────────────────────────────

def _transducer_section(model: ConfigModel, metrics: dict, ti: int, n: int) -> Section:
    title = "Transducer" + (f" {ti + 1}" if n > 1 else "")
    s = Section(title)
    dev = None
    try:
        dev = model.device(ti)
    except (KeyError, IndexError):  # no such transducer: no library entry
        dev = None
    serial = _tx(model, "serial", ti)
    name = _tx(model, "name", ti)
    s.rows.append(Row("Manufacturer", dev.manufact if dev is not None and dev.manufact else NA, "library"))
    s.rows.append(Row("Model / serial", " · ".join(str(x) for x in (
        dev.name if dev is not None and dev.name else None, serial if not is_empty(serial) else None,
        name if not is_empty(name) and name != serial else None) if x) or NA, "config / library"))
    kind = _tx(model, "type", ti)
    s.rows.append(Row("Type", str(kind) if not is_empty(kind) else NA, "config"))
    elem_n = _num(_tx(model, "annular.elem_n", ti)) if kind != "matrix" else None
    if elem_n is None and dev is not None and dev.n_elem:
        elem_n = float(dev.n_elem)
    s.rows.append(Row("Number of elements", fmt(elem_n), "config / library"))
    od = _tx(model, "annular.elem_od_mm", ti) if kind != "matrix" else _tx(model, "matrix.outer_diameter_mm", ti)
    ods = [_num(x) for x in (od if isinstance(od, list) else [od])]
    ods = [x for x in ods if x is not None]
    s.rows.append(Row("Aperture (outer diameter)", fmt(max(ods) if ods else None, "mm"), "config / library"))
    sub = "matrix" if kind == "matrix" else "annular"
    s.rows.append(Row("Radius of curvature", fmt(_num(_tx(model, f"{sub}.curv_radius_mm", ti)), "mm"),
                      "config / library"))
    s.rows.append(Row("Exit plane to geometric focus", fmt(_num(_tx(model, f"{sub}.dist_geom_ep_mm", ti)), "mm"),
                      "config / library"))
    f_run = _num(metrics.get("freq_Hz"))
    f_cfg = _num(_tx(model, "freq_hz", ti))
    s.rows.append(Row("Fundamental frequency", fmt(f_run if f_run is not None else f_cfg, "kHz", 1e3),
                      "run CSV" if f_run is not None else "config"))
    s.rows.append(Row("Focal depth setting (from exit plane)",
                      fmt(_num(_tx(model, "focal_distance_ep", ti)), "mm"), "config"))
    s.rows.append(Row("Free-water ISPPA the drive is scaled to",
                      fmt(_num(_tx(model, "target_isppa_wcm2", ti)), "W/cm²"), "config"))
    return s


def _placement_section(model: ConfigModel, run: rs.Run, metrics: dict, targeting: dict) -> Section:
    s = Section("Placement and targeting")
    rec = placement_record(run)
    mode = rec.get("mode") or targeting.get("placement_mode")
    s.rows.append(Row("Placement method", str(mode) if mode else NA,
                      "placement record" if rec.get("mode") else ("targeting CSV" if mode else "")))
    txs = rec.get("transducers")
    txs = [txs] if isinstance(txs, dict) else (txs if isinstance(txs, list) else [])
    for i, t in enumerate(txs):
        t = t if isinstance(t, dict) else {}
        k = f" (transducer {i + 1})" if len(txs) > 1 else ""
        s.rows.append(Row("Transducer position, bowl apex" + k, _vec(t.get("trans_pos")) + " T1 voxel",
                          "placement record") if t.get("trans_pos") else
                      Row("Transducer position, bowl apex" + k, NA, "placement record"))
        s.rows.append(Row("Planned focus (target)" + k, _vec(t.get("focus_pos")) + " T1 voxel",
                          "placement record") if t.get("focus_pos") else
                      Row("Planned focus (target)" + k, NA, "placement record"))
    if not txs:
        s.rows.append(Row("Transducer and focus positions", NA, "placement record"))
    mni = [targeting.get(f"target_mni_{a}") for a in "xyz"]
    if all(_num(v) is not None for v in mni):
        s.rows.append(Row("Intended target (MNI)", _vec([float(v) for v in mni]) + " mm",
                          "targeting CSV" + (f" ({targeting['target_mni_source']})"
                                             if targeting.get("target_mni_source") else "")))
    s.rows.append(Row("Simulated focus: distance from exit plane to peak",
                      fmt(_num(metrics.get("real_focal_distance_ep_mm")), "mm"), "run CSV"))
    s.rows.append(Row("Simulated focus: FLHM centre from exit plane",
                      fmt(_num(metrics.get("focal_flhm_center_ep_mm")), "mm"), "run CSV"))
    s.rows.append(Row("Focus error (focal centroid to planned focus)",
                      fmt(_num(targeting.get("focus_error_mm")), "mm"), "targeting CSV"))
    if targeting:
        s.rows.append(Row("Focus error, axial / lateral",
                          f"{fmt(_num(targeting.get('focus_error_axial_mm')), 'mm')} / "
                          f"{fmt(_num(targeting.get('focus_error_lateral_mm')), 'mm')}", "targeting CSV"))
    if not targeting:
        s.note = ("Focus errors need the targeting analysis (modules.run_targeting_analysis), "
                  "which writes the targeting CSV.")
    return s


def _timing_section(model: ConfigModel) -> Section:
    s = Section("Pulse timing", note="From the open config; PRF and duty cycles are derived from it.")
    pd, pri = _num(_get(model, "timing.pd")), _num(_get(model, "timing.pri"))
    ptd, ptri, ptrd = (_num(_get(model, f"timing.{k}")) for k in ("ptd", "ptri", "ptrd"))
    pos = lambda x: x if x is not None and x > 0 else None      # noqa: E731 — 0 = not set
    pd, pri, ptd, ptri, ptrd = map(pos, (pd, pri, ptd, ptri, ptrd))
    s.rows.append(Row("Pulse duration (PD)", fmt(pd, "ms", 1e-3), "config"))
    s.rows.append(Row("Pulse repetition interval (PRI)", fmt(pri, "ms", 1e-3), "config"))
    s.rows.append(Row("Pulse repetition frequency (PRF = 1/PRI)", fmt(1 / pri if pri else None, "Hz"), "derived"))
    dc = pd / pri if pd and pri else None
    s.rows.append(Row("Pulse duty cycle (PD/PRI)", fmt(dc * 100 if dc is not None else None, "%"), "derived"))
    s.rows.append(Row("Pulse train duration (PTD)", fmt(ptd, "s"), "config"))
    s.rows.append(Row("Pulse train repetition interval (PTRI)", fmt(ptri, "s"), "config"))
    s.rows.append(Row("Pulse train repetition duration (PTRD)", fmt(ptrd, "s"), "config"))
    tdc = dc * ptd / ptri if dc is not None and ptd and ptri else None
    s.rows.append(Row("Overall duty cycle (PD/PRI × PTD/PTRI)", fmt(tdc * 100 if tdc is not None else None, "%"),
                      "derived"))
    total = ptrd if ptrd else (ptd if ptd else None)
    s.rows.append(Row("Total sonication duration", fmt(total, "s"),
                      "config (PTRD)" if ptrd else ("config (PTD)" if ptd else "config")))
    return s


def _limit_row(label: str, key: str, metrics: dict, source: str, limits: dict, unit: str,
               scale: float = 1.0, column: str | None = None) -> Row:
    v = _num(metrics.get(column or key))
    lim = limits.get(key)
    status = risk_status(v, lim[1]) if lim else ""
    limit = f"≤ {fmt(lim[1], unit, scale)}" if lim else ""
    return Row(label, fmt(v, unit, scale), source, status, limit)


def _width_rows(metrics: dict, measured) -> list[Row]:
    """−3 dB focal widths: the CSV's focal_fwhm_* columns, else what the Nexus measured on
    the T1w intensity map (``measured``: FocalWidths, PENDING or None)."""
    fw = fwm.from_csv(metrics)
    if fw is None and measured == PENDING:
        return [Row("Focal width, axial (−3 dB FWHM)", PENDING, fwm.SOURCE_MAP),
                Row("Focal width, lateral (−3 dB FWHM), narrowest / widest", PENDING, fwm.SOURCE_MAP)]
    if fw is None and isinstance(measured, fwm.FocalWidths):
        fw = measured
    if fw is None:
        return [Row("Focal width, axial (−3 dB FWHM)", NA, ""),
                Row("Focal width, lateral (−3 dB FWHM), narrowest / widest", NA, "")]
    src = fw.source + (f" ({fw.note})" if fw.note else "")
    lat = NA if fw.lateral_min_mm is None or fw.lateral_max_mm is None else \
        f"{fmt(fw.lateral_min_mm)} / {fmt(fw.lateral_max_mm, 'mm')}"
    return [Row("Focal width, axial (−3 dB FWHM)", fmt(fw.axial_mm, "mm"), src),
            Row("Focal width, lateral (−3 dB FWHM), narrowest / widest", lat, src)]


def run_duty_cycle(metrics: dict, model: ConfigModel | None) -> tuple[float | None, str]:
    """(PD/PRI, where it comes from): the run's own (Ispta / Isppa in its CSV, as
    acoustic_analysis.m wrote Ita = Ipa × PD/PRI), else the open config's PD/PRI."""
    ispta, isppa = _num(metrics.get("Ispta")), _num(metrics.get("Isppa"))
    if ispta is not None and isppa:
        return ispta / isppa, "run CSV (Ispta/Isppa)"
    if model is not None:
        pd, pri = _num(_get(model, "timing.pd")), _num(_get(model, "timing.pri"))
        if pd and pd > 0 and pri and pri > 0:
            return pd / pri, "config (PD/PRI)"
    return None, ""


def tissue_ispta(metrics: dict, tissue: str, model: ConfigModel | None) -> tuple[float | None, str]:
    """ISPTA in a tissue [W/cm²]: its ISPPA × PD/PRI (= the tissue maximum of Ita), and
    the source label."""
    isppa = _num(metrics.get(f"Isppa_{tissue}"))
    dc, dc_src = run_duty_cycle(metrics, model)
    if isppa is None or dc is None:
        return None, "derived"
    return isppa * dc, f"derived: ISPPA × duty cycle from {dc_src}"


def _free_water_section(water: dict, metrics: dict, measured=None) -> Section:
    s = Section("Free-water acoustic parameters (simulated)")
    src = "free-water CSV"
    s.rows.append(_limit_row("Peak pressure (spatial peak, temporal peak)", "Psptp", water, src, LIMITS_WATER,
                             "MPa", 1e6))
    s.rows.append(Row("ISPPA", fmt(_num(water.get("Isppa")), "W/cm²"), src))
    s.rows.append(Row(f"ISPTA ({ISPTA_WINDOW})", fmt(_num(water.get("Ispta")), "W/cm²"), src))
    s.rows.append(_limit_row("Mechanical index", "MI", water, src, LIMITS_WATER, ""))
    s.rows.append(Row("Focal volume, −3 dB / −6 dB",
                      f"{fmt(_num(water.get('focal_vol_3dB_mm3')), 'mm³')} / "
                      f"{fmt(_num(water.get('focal_vol_6dB_mm3')), 'mm³')}", src))
    s.rows.append(Row("FLHM centre from exit plane", fmt(_num(water.get("focal_flhm_center_ep_mm")), "mm"), src))
    s.rows.append(Row("Free-water ISPPA reference of this run", fmt(_num(metrics.get("freefield_isppa_wcm2")), "W/cm²"),
                      "run CSV"))
    s.rows += _width_rows(water, measured) if water else _width_rows({}, None)
    if not water:
        s.note = ("No free-water run with this output affix (modules.run_water_baseline writes one).")
    return s


def _in_situ_section(metrics: dict, layered: bool, model: ConfigModel | None = None,
                     measured=None) -> Section:
    s = Section("Estimated in situ acoustic parameters (simulated)")
    src = "run CSV"
    lim = LIMITS_LAYERED if layered else LIMITS_WATER
    if layered:
        for t in TISSUES:
            s.rows.append(_limit_row(f"Peak pressure ({t})", "Psptp", metrics, src, {}, "MPa", 1e6,
                                     column=f"Psptp_{t}"))
        s.rows.append(_limit_row("Peak pressure (whole grid)", "Psptp", metrics, src, lim, "MPa", 1e6))
        for t in TISSUES:
            s.rows.append(Row(f"ISPPA ({t})", fmt(_num(metrics.get(f"Isppa_{t}")), "W/cm²"), src))
        for t in TISSUES:
            v, how = tissue_ispta(metrics, t, model)
            s.rows.append(Row(f"ISPTA ({t}, {ISPTA_WINDOW})", fmt(v, "W/cm²"), how))
        s.rows.append(Row(f"ISPTA (whole grid, {ISPTA_WINDOW})", fmt(_num(metrics.get("Ispta")), "W/cm²"), src))
        s.rows.append(_limit_row("MI, transcranial (MItc)", "MI_tc", metrics, src, lim, ""))
        for t in TISSUES:
            s.rows.append(_limit_row(f"MI ({t})", f"MI_{t}", metrics, src, lim, ""))
        s.rows.append(Row("Focal volume in brain, −3 dB / −6 dB",
                          f"{fmt(_num(metrics.get('focal_vol_3dB_mm3')), 'mm³')} / "
                          f"{fmt(_num(metrics.get('focal_vol_6dB_mm3')), 'mm³')}", src))
        s.rows += _width_rows(metrics, measured)
        s.rows.append(Row("IPA at the target", fmt(_num(metrics.get("Ipa_target")), "W/cm²"), src))
        s.rows.append(Row("Peak pressure at the target", fmt(_num(metrics.get("Ptp_target")), "MPa", 1e6), src))
    else:
        s.note = "A free-water run has no in situ estimate; see the free-water parameters."
    return s


def _thermal_section(metrics: dict, layered: bool, iso: bool) -> Section:
    s = Section("Thermal (simulated)")
    src = "run CSV"
    has = any(k.startswith(("maxT", "riseT", "CEM43")) for k in metrics)
    if not has:
        s.note = "This run has no thermal results (modules.run_heating_sims)."
    if layered:
        for t in TISSUES:
            s.rows.append(_limit_row(f"Max temperature rise ({t})", f"riseT_{t}", metrics, src, LIMITS_LAYERED, "°C"))
        for t in TISSUES:
            s.rows.append(_limit_row(f"Max temperature ({t})", f"maxT_{t}", metrics, src, LIMITS_LAYERED, "°C"))
        for t in TISSUES:
            key = f"CEM43_{t}"
            s.rows.append(_limit_row(f"CEM43 ({t}, {'ISO' if iso else 'k-Wave'} form)", key, metrics, src,
                                     LIMITS_LAYERED, "min", column=ISO_COLUMNS[key] if iso else key))
    else:
        s.rows.append(Row("Max temperature", fmt(_num(metrics.get("maxT")), "°C"), src))
        s.rows.append(Row("Max CEM43", fmt(_num(metrics.get("maxCEM43iso" if iso else "maxCEM43")), "min"), src))
    s.rows.append(Row("Cranial thermal index (TIC, informational)", fmt(_num(metrics.get("TIC"))), src))
    return s


def _safety_section(metrics: dict, water: dict, layered: bool, iso: bool) -> tuple[Section, list[str], dict]:
    """build_safety_banner: every finite limit with a value in the run's CSV."""
    s = Section("ITRUSST non-significant-risk check",
                note="Limits as in PRESTUS's HTML report (get_risk_limits.m). Above a limit means: "
                     "inspect adherence to the NSR criteria that govern the protocol; some criteria "
                     "are alternatives (temperature rise, absolute temperature, or CEM43).")
    limits = LIMITS_LAYERED if layered else LIMITS_WATER
    above, counts = [], {"within": 0, "near": 0, "above": 0}
    for key, (label, limit, unit) in limits.items():
        column = ISO_COLUMNS.get(key, key) if iso else key
        v = _num(metrics.get(column))
        if v is None:
            continue
        scale = 1e6 if unit == "Pa" else 1.0
        shown_unit = "MPa" if unit == "Pa" else unit
        status = risk_status(v, limit)
        counts[status] += 1
        if status == "above":
            above.append(label)
        s.rows.append(Row(label, fmt(v, shown_unit, scale), "run CSV", status, f"≤ {fmt(limit, shown_unit, scale)}"))
    if layered and _num(water.get("MI")) is not None:
        v = _num(water.get("MI"))
        s.rows.append(Row("MI (free water, from the free-water run)", fmt(v), "free-water CSV",
                          risk_status(v, 1.9), "≤ 1.9"))
    if not s.rows:
        s.note = "No metric with an ITRUSST limit in this run's CSV. " + s.note
    return s, above, counts


def width_runs(run: rs.Run, runs: list[rs.Run]) -> list[rs.Run]:
    """The runs whose focal widths the overview of ``run`` needs measured on a T1w map: the
    run itself (if layered) and its free-water run, when their CSV has no focal_fwhm_*
    columns and they have an intensity map and a recorded transducer axis."""
    water = free_water_run(runs, run)
    out = []
    for r in ([run] if run.medium != "water" else []) + [water]:
        if r is None or any(r is o for o in out):
            continue
        if fwm.from_csv(rs.read_metrics(r.csv)) is None and fwm.inputs_for(r) is not None:
            out.append(r)
    return out


def build(model: ConfigModel, run: rs.Run, runs: list[rs.Run],
          widths: dict | None = None) -> Report:
    """The overview of ``run``. ``widths`` maps a run's stem to the focal widths the Nexus
    measured on its T1w map (FocalWidths), or PENDING while that runs; runs whose CSV has
    the focal_fwhm_* columns use those."""
    widths = widths or {}
    metrics = rs.read_metrics(run.csv)
    water_run = free_water_run(runs, run)
    water = rs.read_metrics(water_run.csv) if water_run is not None and water_run is not run else (
        metrics if run.medium == "water" else {})
    targeting = read_table(run.base / f"{run.stem}_targeting.csv")
    layered = run.medium != "water"
    iso = bool(_get(model, "thermal.cem43_iso"))
    safety, above, counts = _safety_section(metrics, water, layered, iso)
    n = max(1, model.n_transducers)
    sections = [safety]
    sections += [_transducer_section(model, metrics, ti, n) for ti in range(n)]
    sections += [_placement_section(model, run, metrics, targeting), _timing_section(model),
                 _free_water_section(water, metrics,
                                     widths.get(water_run.stem) if water_run is not None else None),
                 _in_situ_section(metrics, layered, model, widths.get(run.stem) if layered else None),
                 _thermal_section(metrics, layered, iso)]
    verdict = "none" if not any(counts.values()) else ("review" if above else "within")
    return Report(run.label, run.medium, sections, verdict, above, counts)


def to_markdown(report: Report) -> str:
    """The overview as Markdown tables (for the Copy button)."""
    out = [f"# TUS reporting overview: {report.run_label}", ""]
    out.append({"within": "Within ITRUSST non-significant-risk limits.",
                "review": f"Above an ITRUSST NSR reference level: {', '.join(report.above)}. Inspect NSR adherence.",
                "none": "No ITRUSST-limited metric in this run's outputs."}[report.verdict])
    for sec in report.sections:
        out += ["", f"## {sec.title}", ""]
        if sec.note:
            out += [sec.note, ""]
        if not sec.rows:
            continue
        out += ["| Item | Value | Limit | Status | Source |", "|---|---|---|---|---|"]
        for r in sec.rows:
            out.append(f"| {r.label} | {r.value} | {r.limit} | {r.status} | {r.source} |")
    return "\n".join(out) + "\n"
