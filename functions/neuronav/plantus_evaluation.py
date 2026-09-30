#!/usr/bin/env python3
"""
plantus_evaluation.py  –  criteria evaluation of a PlanTUS placement
=====================================================================
Reads plantus_metrics.mat (written by both PRESTUS PlanTUS launchers) and
evaluates the placed vertex:

  * value and percentile of each criterion at the placed vertex,
  * a weight sweep (configured, equal, one criterion heavy, one criterion
    only) re-scored with plantus_common.compute_composite_score, so each
    row is what PRESTUS would select for those weights,
  * diagnostics that flag criteria or settings that cannot discriminate
    between candidate positions,
  * two figures (criterion maps; weight sensitivity), rendered as 2D
    projections so they take seconds.

Outputs, written next to the metrics file unless --out_dir is given:
  plantus_evaluation.json            summary read by generate_placement_report.m
  plantus_weight_sensitivity.csv     one row per weight set
  plantus_criteria_maps.png
  plantus_weight_sensitivity.png

The launchers call evaluate() after saving the metrics (placement.plantus.
evaluate, default true); failures there are logged and never abort a
placement.  Stand-alone use:

  python plantus_evaluation.py <plantus_metrics.mat> [--out_dir DIR]
                               [--label TEXT] [--no_figures]
"""

import argparse
import csv
import json
import os
import sys
import time

import numpy as np
import scipy.io

from plantus_common import (compute_composite_score, EVALUATION_FILENAME,
                            SWEEP_FILENAME, MAPS_FILENAME, SENSITIVITY_FILENAME)

# key, label, unit, higher_is_better (order used by compute_composite_score)
CRITERIA = [
    ("dist",  "Distance deviation |d − d_opt|", "mm",  False),
    ("ang",   "Skin–target angle",              "deg", False),
    ("inter", "Beam–target intersection",       "mm",  True),
    ("skull", "Skin–skull angle",               "deg", False),
    ("thick", "Skull thickness",                "mm",  False),
]
CRIT_KEYS = [c[0] for c in CRITERIA]
N_LEGACY = 4  # metrics files written before skull thickness hold only four

# Diagnostic thresholds
DIST_FILTER_MIN_EXCLUDED = 0.05  # distance range should exclude ≥5% of non-avoided vertices
INTERSECTION_MIN_FRACTION = 0.02 # ≥2% of eligible vertices should hit the target
NAN_FRACTION_WARN = 0.25         # criterion undefined at >25% of eligible vertices
SENSITIVITY_MM = 20.0            # "heavy" weight set moving the spot further than this
TRADEOFF_PCT = 50.0              # placed spot worse than this % of eligible spots


# ---------------------------------------------------------------------------
# Loading and scoring
# ---------------------------------------------------------------------------

def load_metrics(path):
    """Return a dict with the arrays and settings needed for evaluation."""
    m = scipy.io.loadmat(path, squeeze_me=True)
    d = np.asarray(m["skin_target_distances"], float)
    dopt = float(m["optimal_distance"])
    raw = {
        "dist":  np.abs(d - dopt),
        "ang":   np.asarray(m["skin_target_angles"], float),
        "inter": np.asarray(m["skin_target_intersections"], float),
        "skull": np.asarray(m["skin_skull_angles"], float),
    }
    criteria = CRITERIA[:N_LEGACY]
    if "skull_thickness" in m:
        raw["thick"] = np.asarray(m["skull_thickness"], float)
        criteria = CRITERIA
    best = int(m["best_vertex"]) if "best_vertex" in m else -1
    return {
        "coords":    np.asarray(m["skin_coordinates"], float),
        "faces":     np.asarray(m["skin_faces"], int),
        "target":    np.atleast_1d(np.asarray(m.get("target_center", [np.nan] * 3), float)),
        "d":         d,
        "raw":       raw,
        "criteria":  criteria,
        "avoidance": np.asarray(m["avoidance"], float),
        "dmin":      float(m["min_distance"]),
        "dmax":      float(m["max_distance"]),
        "dopt":      dopt,
        "weights":   np.atleast_1d(np.asarray(m["weights"], float)),
        "placed":    int(m["selected_vertex"]) if "selected_vertex" in m else best,
        "selection_mode": str(m["selection_mode"]) if "selection_mode" in m else "auto",
    }


def make_scorer(M):
    raw = M["raw"]

    def score_for(w):
        return compute_composite_score(
            M["d"], raw["ang"], raw["inter"], raw["skull"], M["avoidance"],
            *w, M["dmax"], min_distance=M["dmin"], optimal_distance=M["dopt"],
            skull_thickness=raw.get("thick"))
    return score_for


def build_weight_sets(configured, n_crit):
    """Ordered (name, [w_dist, w_ang, w_inter, w_skull, w_thick]) pairs.

    Only the first n_crit criteria are swept (n_crit = 4 for legacy files).
    "heavy" gives one criterion 0.55 and splits 0.45 over the others.
    """
    keys = CRIT_KEYS[:n_crit]
    pad = [0.0] * (5 - n_crit)
    sets = [("configured", [float(x) for x in configured]),
            ("equal", [1.0 / n_crit] * n_crit + pad)]
    for i, key in enumerate(keys):
        w = [0.45 / (n_crit - 1)] * n_crit + pad
        w[i] = 0.55
        sets.append((f"{key}_heavy", w))
    for i, key in enumerate(keys):
        w = [0.0] * 5
        w[i] = 1.0
        sets.append((f"{key}_only", w))
    return sets


def percentile_rank(values, idx, valid, higher_is_better):
    """% of valid vertices strictly better than vertex idx (0 = best).

    NaN counts as worst, as in the composite score."""
    v = values[valid]
    x = values[idx]
    if not np.isfinite(x):
        better = np.isfinite(v)
    else:
        better = (v > x) if higher_is_better else (v < x)
    return 100.0 * better.sum() / max(len(v), 1)


def _num(x):
    """JSON-safe float (NaN/inf → None)."""
    x = float(x)
    return x if np.isfinite(x) else None


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def run_sweep(M, score_for):
    """Weight sweep. Returns (rows, chosen, ref_score, valid)."""
    coords, raw = M["coords"], M["raw"]
    ref_score = score_for(M["weights"])
    valid = np.isfinite(ref_score)
    placed = M["placed"]
    rows, chosen = [], []
    for name, w in build_weight_sets(M["weights"], len(M["criteria"])):
        s = score_for(w)
        b = int(np.argmin(s)) if np.isfinite(s).any() else -1
        chosen.append((name, w, s, b))
        if b < 0:
            continue
        row = {
            "weight_set": name,
            **{f"w_{k}": w[i] for i, k in enumerate(CRIT_KEYS)},
            "vertex": b,
            "x_ras": coords[b, 0], "y_ras": coords[b, 1], "z_ras": coords[b, 2],
            "skin_target_distance_mm": M["d"][b],
            "score": s[b],
            "rank_under_configured": int((ref_score[valid] < ref_score[b]).sum()) + 1,
            "dist_to_placed_mm": float(np.linalg.norm(coords[b] - coords[placed])),
        }
        for key, _, _, hib in M["criteria"]:
            row[f"{key}_value"] = raw[key][b]
            row[f"{key}_pct_better"] = percentile_rank(raw[key], b, valid, hib)
        rows.append(row)
    return rows, chosen, ref_score, valid


def diagnostics(M, ref_score, valid, rows):
    """List of {level, code, message} flags for the report."""
    out = []

    def flag(level, code, msg):
        out.append({"level": level, "code": code, "message": msg})

    d, raw, placed = M["d"], M["raw"], M["placed"]
    free = M["avoidance"] != 0
    in_range = (d >= M["dmin"]) & (d <= M["dmax"])
    n_free = int(free.sum())
    if n_free:
        excl = 1.0 - (free & in_range).sum() / n_free
        if excl < DIST_FILTER_MIN_EXCLUDED:
            flag("warning", "distance_range_unconstrained",
                 f"The distance range [{M['dmin']:.1f}, {M['dmax']:.1f}] mm excludes only "
                 f"{100 * excl:.1f}% of the skin outside the avoidance mask, so it does not "
                 "limit where the transducer can go. It is set from focal_distance_list ± "
                 "max(flhm_list); check the FLHM values.")

    n_valid = int(valid.sum())
    if n_valid == 0:
        flag("error", "no_eligible", "No eligible vertex (all excluded).")
        return out

    frac_hit = float((raw["inter"][valid] > 0).mean())
    if frac_hit < INTERSECTION_MIN_FRACTION:
        flag("warning", "intersection_sparse",
             f"The beam–target intersection is non-zero at only {100 * frac_hit:.1f}% of "
             "eligible vertices, so this criterion acts almost as an on/off flag. "
             "A larger target mask (placement.plantus.mask_radius_mm) makes it graded.")

    for key, lab, _, _ in M["criteria"]:
        if key == "dist":
            continue
        nan_frac = float((~np.isfinite(raw[key][valid])).mean())
        if nan_frac == 1.0:
            flag("warning", f"{key}_unavailable",
                 f"{lab} is undefined at every eligible vertex; its weight is ignored.")
        elif nan_frac > NAN_FRACTION_WARN:
            flag("info", f"{key}_nan",
                 f"{lab} is undefined at {100 * nan_frac:.0f}% of eligible vertices "
                 "(scored as worst).")

    if not np.isfinite(ref_score[placed]):
        reasons = []
        if not free[placed]:
            reasons.append("inside the avoidance mask")
        if d[placed] < M["dmin"]:
            reasons.append(f"closer than min_distance ({d[placed]:.1f} mm)")
        if d[placed] > M["dmax"]:
            reasons.append(f"further than max_distance ({d[placed]:.1f} mm)")
        flag("warning", "placed_ineligible",
             "The placed vertex is not eligible under the configured scoring: "
             + (", ".join(reasons) or "excluded") + ".")
    else:
        rank = int((ref_score[valid] < ref_score[placed]).sum()) + 1
        best = int(np.argmin(ref_score))
        if placed != best:
            flag("info", "placed_not_optimum",
                 f"The placed vertex ({M['selection_mode']} selection) ranks {rank} of "
                 f"{n_valid} eligible vertices, "
                 f"{np.linalg.norm(M['coords'][placed] - M['coords'][best]):.1f} mm from "
                 "the composite optimum.")

    w = M["weights"]
    for i, (key, lab, _, hib) in enumerate(M["criteria"]):
        if w[i] <= 0:
            continue
        pct = percentile_rank(raw[key], placed, valid, hib)
        if pct > TRADEOFF_PCT:
            flag("info", f"tradeoff_{key}",
                 f"At the placed vertex, {lab.lower()} is worse than at {pct:.0f}% of "
                 "eligible vertices (traded off against the other criteria).")

    heavy = [r for r in rows if r["weight_set"].endswith("_heavy")]
    moved = [r for r in heavy if r["dist_to_placed_mm"] > SENSITIVITY_MM]
    if moved:
        names = ", ".join(f"{r['weight_set']} ({r['dist_to_placed_mm']:.0f} mm)" for r in moved)
        flag("warning", "weight_sensitive",
             f"Emphasising a single criterion moves the chosen spot by more than "
             f"{SENSITIVITY_MM:.0f} mm: {names}. The placement depends on the weighting.")
    return out


def evaluate(metrics_path, out_dir=None, label="", figures=True, log=print):
    """Evaluate a PlanTUS placement; write JSON, CSV and (optionally) figures.

    Returns the summary dict that is written to plantus_evaluation.json.
    """
    t0 = time.time()
    out_dir = out_dir or os.path.dirname(os.path.abspath(metrics_path))
    os.makedirs(out_dir, exist_ok=True)
    if not label:
        parts = os.path.abspath(out_dir).split(os.sep)
        label = " ".join(parts[-2:])

    M = load_metrics(metrics_path)
    score_for = make_scorer(M)
    rows, chosen, ref_score, valid = run_sweep(M, score_for)
    flags = diagnostics(M, ref_score, valid, rows)
    placed, coords, raw = M["placed"], M["coords"], M["raw"]
    n_valid = int(valid.sum())

    placed_info = {
        "vertex": placed,
        "selection_mode": M["selection_mode"],
        "ras_mm": [_num(x) for x in coords[placed]],
        "skin_target_distance_mm": _num(M["d"][placed]),
        "score": _num(ref_score[placed]),
        "rank": (int((ref_score[valid] < ref_score[placed]).sum()) + 1
                 if np.isfinite(ref_score[placed]) else None),
        "criteria": {},
    }
    crit_summary = []
    for i, (key, lab, unit, hib) in enumerate(M["criteria"]):
        v = raw[key][valid]
        fin = v[np.isfinite(v)]
        crit_summary.append({
            "key": key, "label": lab, "unit": unit,
            "higher_is_better": hib, "weight": _num(M["weights"][i]),
            "eligible_min": _num(fin.min()) if len(fin) else None,
            "eligible_median": _num(np.median(fin)) if len(fin) else None,
            "eligible_max": _num(fin.max()) if len(fin) else None,
            "eligible_defined_fraction": _num(len(fin) / max(len(v), 1)),
        })
        placed_info["criteria"][key] = {
            "value": _num(raw[key][placed]),
            "pct_better": _num(percentile_rank(raw[key], placed, valid, hib)),
        }

    heavy_only = [r for r in rows if r["weight_set"] != "configured"]
    summary = {
        "label": label,
        "metrics_file": os.path.abspath(metrics_path),
        "n_vertices": int(len(coords)),
        "n_eligible": n_valid,
        "distance_range_mm": [M["dmin"], M["dmax"]],
        "optimal_distance_mm": M["dopt"],
        "target_ras_mm": [_num(x) for x in M["target"]],
        "criteria": crit_summary,
        "placed": placed_info,
        "composite_optimum_vertex": int(np.argmin(ref_score)) if n_valid else None,
        "sweep": [{k: (_num(v) if isinstance(v, (float, np.floating)) else
                       int(v) if isinstance(v, (np.integer,)) else v)
                   for k, v in r.items()} for r in rows],
        "sweep_distinct_spots": len({r["vertex"] for r in rows}),
        "sweep_max_dist_to_placed_mm": _num(max((r["dist_to_placed_mm"] for r in heavy_only),
                                                default=0.0)),
        "flags": flags,
        "figures": {},
    }

    with open(os.path.join(out_dir, SWEEP_FILENAME), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["weight_set"])
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.4g}" if isinstance(v, (float, np.floating)) else v)
                         for k, v in r.items()})

    if figures:
        try:
            summary["figures"] = render_figures(M, ref_score, valid, rows, chosen,
                                                out_dir, label)
        except Exception as exc:  # noqa: BLE001 — figures are optional
            log(f"Warning: evaluation figures not written: {exc}")

    summary["elapsed_s"] = round(time.time() - t0, 2)
    with open(os.path.join(out_dir, EVALUATION_FILENAME), "w") as fh:
        json.dump(summary, fh, indent=1)

    n_warn = sum(f["level"] == "warning" for f in flags)
    log(f"Evaluation: placed vertex {placed} rank {placed_info['rank']} of {n_valid}; "
        f"{summary['sweep_distinct_spots']} distinct spots across {len(rows)} weight sets "
        f"(max {summary['sweep_max_dist_to_placed_mm'] or 0:.0f} mm from placed); "
        f"{n_warn} warning(s). [{summary['elapsed_s']:.1f} s]")
    for f in flags:
        if f["level"] != "info":
            log(f"  {f['level']}: {f['message']}")
    return summary


# ---------------------------------------------------------------------------
# Figures (orthographic 2D projection: fast, no mplot3d)
# ---------------------------------------------------------------------------

def _face_normals(coords, faces):
    v0, v1, v2 = coords[faces[:, 0]], coords[faces[:, 1]], coords[faces[:, 2]]
    n = np.cross(v1 - v0, v2 - v0)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    centres = (v0 + v1 + v2) / 3.0
    n[np.einsum("ij,ij->i", n, centres - coords.mean(axis=0)) < 0] *= -1
    return n, centres


def _vertex_normals(coords, faces, fn):
    vn = np.zeros_like(coords)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    return vn / (np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12)


def _camera(view):
    """(view_dir, up) unit vectors; view_dir points from head to camera."""
    v = np.asarray(view, float)
    v /= np.linalg.norm(v)
    up = np.array([0.0, 0.0, 1.0])
    if abs(v @ up) > 0.95:          # looking from above/below: anterior up
        up = np.array([0.0, 1.0, 0.0])
    right = np.cross(-v, up)
    right /= np.linalg.norm(right)
    up = np.cross(right, -v)
    return v, right, up


def _draw(ax, geo, values, cmap, norm, view, markers=()):
    """Draw the skin coloured by per-vertex values (NaN = grey) seen from view."""
    from matplotlib.collections import PolyCollection
    coords, faces, fn, vn, centre, radius = geo
    v, right, up = _camera(view)
    vis = fn @ v > 0
    f = faces[vis]
    depth = (coords[f].mean(axis=1) @ v)
    order = np.argsort(depth)       # far first
    f = f[order]
    xy = np.stack([(coords - centre) @ right, (coords - centre) @ up], axis=1)
    fv = values[f]
    invalid = ~np.isfinite(fv).all(axis=1)
    rgba = cmap(norm(np.where(invalid, 0.0, np.nan_to_num(fv).mean(axis=1))))
    rgba[invalid] = (0.80, 0.80, 0.80, 1.0)
    shade = 0.55 + 0.45 * np.clip(fn[vis][order] @ v, 0, 1)
    rgba[:, :3] *= shade[:, None]
    ax.add_collection(PolyCollection(xy[f], facecolors=rgba, edgecolors=rgba,
                                     linewidths=0.2))
    for vidx, kw in markers:
        if vidx is None or vidx < 0 or vn[vidx] @ v <= 0:
            continue                # far side of the head
        ax.scatter(*xy[vidx], zorder=5, **kw)
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_aspect("equal")
    ax.set_axis_off()


def _norm(vals):
    from matplotlib import colors
    fin = vals[np.isfinite(vals)]
    if len(fin) == 0:
        return colors.Normalize(0, 1)
    lo, hi = np.percentile(fin, [1, 99])
    if hi <= lo:
        lo, hi = fin.min(), fin.max()
    if hi <= lo:
        hi = lo + 1
    return colors.Normalize(lo, hi)


def render_figures(M, ref_score, valid, rows, chosen, out_dir, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.lines import Line2D

    coords, faces = M["coords"], M["faces"]
    fn, _ = _face_normals(coords, faces)
    vn = _vertex_normals(coords, faces, fn)
    centre = coords.mean(axis=0)
    radius = 0.5 * (coords.max(axis=0) - coords.min(axis=0)).max() * 1.02
    geo = (coords, faces, fn, vn, centre, radius)
    placed = M["placed"]
    best = int(np.argmin(ref_score)) if valid.any() else -1

    views = [("placement", coords[placed] - centre), ("right", [1, 0, 0]),
             ("left", [-1, 0, 0]), ("top", [0, 0, 1])]
    marks = [(placed, dict(s=110, c="red", marker="*", edgecolors="k", linewidths=0.6))]
    if best != placed:
        marks.append((best, dict(s=55, c="white", marker="o", edgecolors="k", linewidths=0.8)))

    # --- criterion maps ---------------------------------------------------
    maps = []
    for key, lab, unit, hib in M["criteria"]:
        v = M["raw"][key].copy()
        v[~valid] = np.nan
        maps.append((f"{lab} [{unit}]\n{'higher' if hib else 'lower'} is better",
                     v, "viridis" if hib else "viridis_r"))
    comp = ref_score.copy()
    comp[~valid] = np.nan
    maps.append(("Composite score\n(configured weights)", comp, "viridis_r"))

    nr, nc = len(maps), len(views)
    fig, axes = plt.subplots(nr, nc, figsize=(2.6 * nc + 2.4, 2.5 * nr), squeeze=False)
    for r, (title, vals, cmap_name) in enumerate(maps):
        cmap, norm = plt.get_cmap(cmap_name), _norm(vals)
        for c, (vname, view) in enumerate(views):
            _draw(axes[r, c], geo, vals, cmap, norm, view, marks)
            if r == 0:
                axes[r, c].set_title(vname, fontsize=9)
        axes[r, 0].text(-0.08, 0.5, title, transform=axes[r, 0].transAxes, rotation=90,
                        ha="right", va="center", fontsize=8)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[r, :].tolist(),
                     fraction=0.02, pad=0.01)
    legend = "grey = excluded   ★ placed" + ("   ○ composite optimum" if best != placed else "")
    fig.suptitle(f"{label} — PlanTUS criteria ({legend})", fontsize=10, y=0.995)
    fig.subplots_adjust(top=0.96)
    p_maps = os.path.join(out_dir, MAPS_FILENAME)
    fig.savefig(p_maps, dpi=100, bbox_inches="tight")
    plt.close(fig)

    # --- weight sensitivity ----------------------------------------------
    crit = M["criteria"]
    keys = [c[0] for c in crit]
    # first row: the placed vertex, for comparison with every weight set
    placed_row = {"weight_set": f"placed ({M['selection_mode']})", "dist_to_placed_mm": 0.0}
    for key, _, _, hib in crit:
        placed_row[f"{key}_value"] = M["raw"][key][placed]
        placed_row[f"{key}_pct_better"] = percentile_rank(M["raw"][key], placed, valid, hib)
    table = [placed_row] + rows
    pct = np.array([[r[f"{k}_pct_better"] for k in keys] for r in table])
    fig = plt.figure(figsize=(12.5, 0.42 * len(table) + 2.2))
    ax = fig.add_axes([0.13, 0.14, 0.43, 0.74])
    im = ax.imshow(pct, cmap="magma_r", vmin=0, vmax=100, aspect="auto")
    for i, r in enumerate(table):
        for j, (key, _, unit, _) in enumerate(crit):
            val = r[f"{key}_value"]
            txt = "n/a" if not np.isfinite(val) else f"{val:.1f} {unit}"
            ax.text(j, i, f"{txt}\n{pct[i, j]:.0f}%", ha="center", va="center",
                    fontsize=6.5, color="white" if pct[i, j] > 55 else "black")
    ax.set_xticks(range(len(crit)))
    ax.set_xticklabels([c[1].replace(" ", "\n", 1) for c in crit], fontsize=7)
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels([f"{r['weight_set']}  (Δ {r['dist_to_placed_mm']:.0f} mm)"
                        for r in table], fontsize=7)
    ax.axhline(0.5, color="k", lw=1.2)
    ax.set_title("Placed spot and spot chosen per weight set (Δ = distance to placed):\n"
                 "value and % of eligible spots that are better", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, label="% eligible spots better (0 = best)")

    tab = plt.get_cmap("tab10")
    spot_marks = [(r["vertex"], dict(s=48, color=tab(i % 10), marker="o",
                                     edgecolors="k", linewidths=0.5))
                  for i, r in enumerate(rows)]
    spot_marks.append((placed, dict(s=130, c="red", marker="*", edgecolors="k",
                                    linewidths=0.6)))
    base = np.zeros(len(coords))
    base[~valid] = np.nan
    base_cmap = colors.ListedColormap(["#e9e4d8"])
    for k, (vname, view) in enumerate([views[0], views[3]]):   # placement, top
        axk = fig.add_axes([0.62 + 0.19 * k, 0.2, 0.18, 0.66])
        _draw(axk, geo, base, base_cmap, colors.Normalize(0, 1), view, spot_marks)
        axk.set_title(vname, fontsize=8)
    handles = [Line2D([], [], ls="", marker="o", color=tab(i % 10), markeredgecolor="k",
                      label=r["weight_set"]) for i, r in enumerate(rows)]
    handles.append(Line2D([], [], ls="", marker="*", color="red", markeredgecolor="k",
                          markersize=10, label="placed"))
    fig.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.99, 0.0),
               fontsize=6.5, frameon=False, ncol=4)
    fig.suptitle(f"{label} — weight sensitivity "
                 f"({len({r['vertex'] for r in rows})} distinct spots; far-side spots hidden)",
                 fontsize=10)
    p_sens = os.path.join(out_dir, SENSITIVITY_FILENAME)
    fig.savefig(p_sens, dpi=100)
    plt.close(fig)
    return {"criteria_maps": p_maps, "weight_sensitivity": p_sens}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description="Evaluate a PlanTUS placement.")
    p.add_argument("metrics_mat")
    p.add_argument("--out_dir", default=None)
    p.add_argument("--label", default="")
    p.add_argument("--no_figures", action="store_true")
    a = p.parse_args()
    if not os.path.isfile(a.metrics_mat):
        sys.exit(f"Not found: {a.metrics_mat}")
    evaluate(a.metrics_mat, a.out_dir, a.label, figures=not a.no_figures)


if __name__ == "__main__":
    main()
