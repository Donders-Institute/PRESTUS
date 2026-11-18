#!/usr/bin/env python3
"""Write the data of the load calculator in the docs (doc_hpc.md, Resource Estimates).

The calculator (doc/javascripts/load_predictor.js) shows the time of each stage of a run
for a GPU, solver, resolution, and medium, as the Nexus predicts it (workload.py) from
nexus/prestus_nexus/throughput.yaml. This script resolves the rates per GPU and solver
(the table's own entries, else the reference solver times the code-type and precision
factors, else for a GPU without entries the same-generation entries scaled by memory
bandwidth) and writes them, with the stage times and the grid constants, to
doc/javascripts/load_model_data.js. Needs PyYAML only (no Nexus imports): the constants
it shares with workload.py are read from that file's source.

Usage:
    python tests/docs/gen_load_data.py           # (re)write doc/javascripts/load_model_data.js
    python tests/docs/gen_load_data.py --check   # exit 1 if it is missing or out of date
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import statistics
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "nexus" / "prestus_nexus" / "throughput.yaml"
WORKLOAD = ROOT / "nexus" / "prestus_nexus" / "workload.py"
DEFAULTS = ROOT / "config" / "config_default.yaml"
OUT = ROOT / "doc" / "javascripts" / "load_model_data.js"

SOLVERS = {"gpu": ("cpp_gpu", "matlab_gpu"), "cpu": ("cpp_cpu", "matlab_cpu")}
PRECISIONS = ("single", "double")
STAGES = ("placement", "preproc", "preproc_shared", "medium", "nifti_medium", "source",
          "acoustic_analysis", "nifti_acoustic", "thermal_analysis", "nifti_thermal",
          "segmentation", "pct", "report")


def constants() -> dict:
    """REF_HEAD_MM, REF_HEAD_L, DEFAULT_T0_S from workload.py."""
    out = {}
    for node in ast.parse(WORKLOAD.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name in ("REF_HEAD_MM", "REF_HEAD_L", "DEFAULT_T0_S"):
                out[name] = ast.literal_eval(node.value)
    return out


def norm(text: str) -> str:
    """workload._norm"""
    t = re.sub(r"[\s_]+", "-", str(text).lower())
    return t.replace("nvidia-", "").replace("tesla-", "")


def best_match(items: list[dict], hint: str) -> dict | None:
    """workload._best_match: the item whose longest match key occurs in hint."""
    best, length = None, 0
    for e in items:
        for m in e.get("match") or []:
            k = norm(m)
            if k and k in hint and len(k) > length:
                best, length = e, len(k)
    return best


def entry_rate(e: dict, factor: float = 1.0, note: str = "") -> dict:
    """workload._entry_rate"""
    lo, hi = (e.get("range") or [e["ns"], e["ns"]])[:2]
    return {"ns": e["ns"] * factor, "lo": lo * factor, "hi": hi * factor, "t0": float(e.get("t0_s") or 0.0),
            "basis": f"{e['device']} (n={e.get('n', '?')}){note}"}


def fallbacks(pool: list[dict], want: dict, kind: str, factors: dict) -> list[tuple[list[dict], float, str]]:
    """workload._fallbacks: this solver, then the reference solver x factors."""
    fits = lambda e, w: all(k not in e or e[k] == v for k, v in w.items())   # noqa: E731
    out = [([e for e in pool if fits(e, want)], 1.0, "")]
    over, notes, factor = {}, [], 1.0
    ref = "cpp_gpu" if kind == "gpu" else "cpp_cpu"
    fx = (factors.get("code_type") or {}).get(want["code_type"])
    if want["code_type"] != ref and isinstance(fx, dict) and fx.get("x"):
        over["code_type"], factor = ref, factor * float(fx["x"])
        notes.append(f"x{float(fx['x']):.1f} for {want['code_type']}")
    fp = (factors.get("precision") or {}).get(want["precision"])
    if want["precision"] != "single" and isinstance(fp, dict) and fp.get("x"):
        over["precision"], factor = "single", factor * float(fp["x"])
        notes.append(f"x{float(fp['x']):.1f} for {want['precision']} precision")
    if over:
        out.append(([e for e in pool if fits(e, {**want, **over})], factor, "; " + ", ".join(notes)))
    return [c for c in out if c[0]]


def bandwidth_rate(cands: list[dict], hint: str, bandwidth: list[dict]) -> dict | None:
    """workload._bandwidth_rate"""
    target = best_match(bandwidth, hint)
    if target is None or not target.get("gbs"):
        return None
    k, devs, t0 = [], [], []
    for e in cands:
        spec = best_match(bandwidth, " ".join(norm(m) for m in e.get("match") or [e["device"]]))
        if spec is not None and spec.get("arch") == target.get("arch") and spec.get("gbs"):
            k.append(float(e["ns"]) * float(spec["gbs"]))
            devs.append(e["device"])
            t0.append(float(e.get("t0_s") or 0.0))
    if not k:
        return None
    g = float(target["gbs"])
    name = target.get("name") or (target.get("match") or ["this GPU"])[0]
    return {"ns": statistics.median(k) / g, "lo": min(k) / g * 0.8, "hi": max(k) / g * 1.25,
            "t0": statistics.median(t0), "basis": f"{name}: scaled by memory bandwidth from {', '.join(devs)}"}


def resolve(table: dict, kind: str, solver: str, precision: str, hint: str) -> dict | None:
    """workload.table_rate_for a kwave, non-axisymmetric config on the device named by hint."""
    want = {"code_type": solver, "backend": "kwave", "axisymmetric": 0, "precision": precision}
    pool = [e for e in table.get("acoustic") or [] if e.get("kind") == kind and e.get("ns")]
    options = fallbacks(pool, want, kind, table.get("factors") or {})
    if not options:
        return None
    for cands, factor, note in options:
        hit = best_match(cands, hint)
        if hit is not None:
            return entry_rate(hit, factor, note)
    if kind == "gpu":
        for cands, factor, note in options:
            bw = bandwidth_rate(cands, hint, table.get("bandwidth") or [])
            if bw is not None:
                return {**bw, "ns": bw["ns"] * factor, "lo": bw["lo"] * factor, "hi": bw["hi"] * factor,
                        "basis": bw["basis"] + note}
    cands, factor, note = options[0]
    if len(cands) == 1:
        return entry_rate(cands[0], factor, note)
    ns = [float(e["ns"]) for e in cands]
    runs = sum(int(e.get("n") or 0) for e in cands)
    what = (f"cluster CPU runs (n={runs}; the core count is not modelled)" if kind == "cpu"
            else f"any GPU ({', '.join(e['device'] for e in cands)})")
    return {"ns": statistics.median(ns) * factor, "lo": min(ns) * factor, "hi": max(ns) * factor,
            "t0": statistics.median(float(e.get("t0_s") or 0.0) for e in cands), "basis": what + note}


def thermal_rate(table: dict, kind: str, hint: str) -> dict | None:
    """workload.table_rate for the thermal stage (no solver keys)."""
    cands = [e for e in table.get("thermal") or [] if e.get("kind") == kind and e.get("ns")]
    hit = best_match(cands, hint)
    if hit is None:
        if len(cands) != 1:
            return None
        hit = cands[0]
    r = entry_rate(hit)
    r.pop("t0")
    return r


def devices(table: dict) -> list[dict]:
    """Every GPU with an acoustic entry or a bandwidth spec, and the CPU, with resolved rates."""
    bandwidth = table.get("bandwidth") or []
    out, seen = [], set()
    for e in table.get("acoustic") or []:
        if e.get("kind") == "gpu" and e.get("backend", "kwave") == "kwave" and not e.get("axisymmetric") \
                and e["device"] not in seen:
            seen.add(e["device"])
            hint = " ".join(norm(m) for m in e.get("match") or [e["device"]])
            out.append({"id": e["device"], "label": re.sub(r"^(NVIDIA|Tesla) ", "", e["device"]),
                        "kind": "gpu", "measured": True, "hint": hint})
    covered = {id(best_match(bandwidth, d["hint"])) for d in out}
    for b in bandwidth:
        if id(b) not in covered:
            name = b.get("name") or b["match"][0]
            out.append({"id": name, "label": f"{name} (estimated from bandwidth)", "kind": "gpu",
                        "measured": False, "hint": norm(b["match"][0])})
    out.sort(key=lambda d: (not d["measured"], d["label"].lower()))
    out.append({"id": "CPU", "label": "CPU", "kind": "cpu", "measured": True, "hint": ""})
    for d in out:
        d["rates"] = {}
        for solver in SOLVERS[d["kind"]]:
            for precision in PRECISIONS:
                r = resolve(table, d["kind"], solver, precision, d["hint"])
                if r is not None:
                    d["rates"][f"{solver}|{precision}"] = {k: (round(v, 4) if isinstance(v, float) else v)
                                                           for k, v in r.items()}
        d["thermal"] = thermal_rate(table, d["kind"], d["hint"])
    # a GPU of a generation without measurements gets only the span of all GPUs: leave it out
    return [d for d in out if d["rates"] and (d["measured"] or not any(
        r["basis"].startswith("any ") for r in d["rates"].values()))]


def speeds(defaults: dict) -> dict:
    """(c_min, c_max) of the layered head (layers in config_default.yaml, plus water) and of water."""
    props = defaults.get("medium_properties") or {}
    layers = set(defaults.get("layers") or {}) | {"water"}
    c = [float(v["sound_speed"]) for k, v in props.items()
         if k in layers and isinstance(v, dict) and isinstance(v.get("sound_speed"), (int, float))]
    water = float(props["water"]["sound_speed"])
    return {"head": [min(c), max(c)], "water": [water, water]}


def build() -> str:
    table = yaml.safe_load(TABLE.read_text(encoding="utf-8")) or {}
    defaults = yaml.safe_load(DEFAULTS.read_text(encoding="utf-8")) or {}
    grid = defaults.get("grid") or {}
    const = constants()
    stages = {}
    for e in table.get("stages") or []:
        if e.get("stage") in STAGES and e.get("medium") in ("head", "water"):
            stages.setdefault(e["medium"], {})[e["stage"]] = {
                k: e[k] for k in ("s_1mm", "exponent", "ratio", "n")}
    data = {
        "source": "nexus/prestus_nexus/throughput.yaml",
        "ref_head_mm": list(const["REF_HEAD_MM"]), "ref_head_l": list(const["REF_HEAD_L"]),
        "t0_default": const["DEFAULT_T0_S"],
        "cfl": float(grid.get("source_cfl") or 0.15), "transits": 1.0,
        "water_margin_mm": float(grid.get("water_margin_mm") or 5.0),
        "z_beyond_focus_mm": float(grid.get("z_beyond_focus_mm") or 40.0),
        "speeds": speeds(defaults),
        "devices": devices(table),
        "stages": stages,
    }
    body = json.dumps(data, indent=1, sort_keys=True)
    return ("// Generated by tests/docs/gen_load_data.py from nexus/prestus_nexus/throughput.yaml; do not edit.\n"
            f"window.PRESTUS_LOAD_MODEL = {body};\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="exit 1 if the file is missing or out of date")
    args = ap.parse_args()
    text = build()
    if args.check:
        if not OUT.is_file() or OUT.read_text(encoding="utf-8") != text:
            print(f"{OUT.relative_to(ROOT)} is out of date: run python tests/docs/gen_load_data.py", file=sys.stderr)
            return 1
        print(f"{OUT.relative_to(ROOT)} is up to date")
        return 0
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(json.loads(text.split('=', 1)[1].rstrip().rstrip(';'))['devices'])} devices")
    return 0


if __name__ == "__main__":
    sys.exit(main())
