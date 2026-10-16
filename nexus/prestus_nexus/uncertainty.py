"""Medium-property sets of the uncertainty variants (functions/core/uncertainty_medium_config.m).

With simulation.uncertainty on, PRESTUS runs the default medium and two variants: conservative
(lower bound of intensity and heating) and liberal (upper bound). Each variant's file is merged
onto the run's config (uncertainty_pipeline.m, make_sim_params), so a tissue or property the
file does not list keeps its medium_properties value.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .merge import is_empty, mergestruct
from .yamlio import load_file

VARIANTS = ("conservative", "liberal")
VARIANT_LABEL = {"default": "Default", "conservative": "Conservative (lower bound)",
                 "liberal": "Liberal (upper bound)"}
BRAIN = ("brain",)
SKULL = ("skull", "skull_trabecular", "skull_cortical")
TISSUE_LABEL = {"brain": "Brain", "skin": "Skin", "water": "Water", "skull": "Skull (single layer)",
                "skull_trabecular": "Skull trabecular", "skull_cortical": "Skull cortical"}
# (key, label, unit); 'attenuation' is derived: alpha_coeff * f^alpha_power at the transducer frequency
PROPERTIES = [
    ("sound_speed", "Sound speed", "m/s"),
    ("density", "Density", "kg/m³"),
    ("attenuation", "Attenuation", "dB/cm"),
    ("alpha_coeff", "Attenuation coefficient", "dB/(cm·MHz^y)"),
    ("alpha_power", "Attenuation power law y", ""),
    ("thermal_conductivity", "Thermal conductivity", "W/m/°C"),
    ("specific_heat_capacity", "Specific heat", "J/kg/°C"),
    ("perfusion", "Perfusion", "mL/min/kg"),
    ("absorption_fraction", "Absorption fraction", ""),
]
PLOTTED = ("sound_speed", "density", "attenuation", "thermal_conductivity",
           "specific_heat_capacity", "perfusion", "absorption_fraction")
FALLBACK_MHZ = 0.5          # the ITRUSST property database compares attenuation at 0.5 MHz


@dataclass
class PropertySets:
    sets: dict[str, dict] = field(default_factory=dict)      # variant -> {tissue: {prop: value}}
    files: dict[str, Path] = field(default_factory=dict)     # variant -> resolved file
    listed: dict[str, dict] = field(default_factory=dict)    # variant -> the file's own medium_properties
    errors: dict[str, str] = field(default_factory=dict)     # variant -> why its file cannot be used
    freq_mhz: float = FALLBACK_MHZ
    freq_from_transducer: bool = False

    def value(self, variant: str, tissue: str, prop: str) -> Any:
        props = self.sets.get(variant, {}).get(tissue, {})
        if prop == "attenuation":
            a, y = props.get("alpha_coeff"), props.get("alpha_power")
            return attenuation_db_cm(a, y, self.freq_mhz)
        return props.get(prop)

    def inherited(self, variant: str, tissue: str, prop: str) -> bool:
        """The variant's value comes from medium_properties, not from its own file."""
        own = self.listed.get(variant, {}).get(tissue, {})
        if prop == "attenuation":
            return "alpha_coeff" not in own and "alpha_power" not in own
        return prop not in own


def attenuation_db_cm(alpha_coeff: Any, alpha_power: Any, f_mhz: float) -> float | None:
    try:
        return float(alpha_coeff) * float(f_mhz) ** float(alpha_power)
    except (TypeError, ValueError):
        return None


def variant_file(model, prestus_root: str | Path, variant: str) -> Path:
    """The file uncertainty_medium_config.m resolves for a variant."""
    folder = Path(prestus_root) / "config" / "uncertainty"
    value = model.get(f"simulation.uncertainty_{variant}")
    if is_empty(value):
        return folder / f"config_medium_{variant}.yaml"
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() or _is_windows_absolute(str(value)) else folder / p


def _is_windows_absolute(s: str) -> bool:
    return len(s) > 2 and s[1] == ":" and s[2] in "\\/" or s.startswith("\\")


def default_properties(model) -> dict:
    """medium_properties as this config sets it (defaults, project defaults, study)."""
    node = model.schema.get("medium_properties")
    out: dict = {}
    for tissue in (node.children if node else []):
        if tissue.is_group:
            out[tissue.name] = {p.name: model.get(p.path) for p in tissue.children if not p.is_group}
    return out


def transducer_mhz(model) -> float | None:
    f = model.get("transducer.freq_hz", 0)
    if isinstance(f, (list, tuple)):
        f = f[0] if f else None
    try:
        return float(f) / 1e6 if f and float(f) > 0 else None
    except (TypeError, ValueError):
        return None


def property_sets(model, prestus_root: str | Path) -> PropertySets:
    base = default_properties(model)
    f = transducer_mhz(model)
    out = PropertySets(sets={"default": base}, freq_mhz=f or FALLBACK_MHZ, freq_from_transducer=f is not None)
    for variant in VARIANTS:
        path = variant_file(model, prestus_root, variant)
        out.files[variant] = path
        if not path.is_file():
            out.errors[variant] = f"file not found: {path}"
            continue
        try:
            data = load_file(path) or {}
        except Exception as exc:
            out.errors[variant] = f"cannot read {path.name}: {exc}"
            continue
        listed = data.get("medium_properties") if isinstance(data, dict) else None
        if not isinstance(listed, dict) or not listed:
            out.errors[variant] = f"{path.name} has no medium_properties block"
            continue
        out.listed[variant] = listed
        out.sets[variant] = mergestruct(base, listed)
    return out
