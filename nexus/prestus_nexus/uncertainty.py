"""Medium-property sets of the uncertainty variants (functions/uncertainty/uncertainty_medium_layer.m).

With simulation.uncertainty on, PRESTUS runs the default medium and four variants in two pairs
(functions/uncertainty/uncertainty_variants.m). Temp-H / Temp-L take, per property, the ITRUSST
bound that maximises / minimises brain heating. Trans-H / Trans-L take, per acoustic property,
the bound that maximises / minimises brain intensity; they run acoustic only, so their files
list no thermal properties. Each variant takes the shipped ITRUSST bounds
(config/uncertainty/config_medium_<variant>.yaml) with the study's edits in
simulation.uncertainty_overrides merged on top, and that layer is merged onto the run's config
(uncertainty_pipeline.m, make_sim_params), so a tissue or property neither lists keeps its
medium_properties value.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .merge import is_empty, mergestruct
from .yamlio import load_file

# ITRUSST tissue-property overview: the reference values and bounds the variants take
ITRUSST_PROPERTIES_URL = "https://itrusst.github.io/itrusst_properties/"

VARIANTS = ("temp_high", "temp_low", "trans_high", "trans_low")
VARIANT_LABEL = {"default": "Default", "temp_high": "Temp-H (max. heating)",
                 "temp_low": "Temp-L (min. heating)", "trans_high": "Trans-H (max. intensity)",
                 "trans_low": "Trans-L (min. intensity)"}
VARIANT_SHORT = {"default": "Default", "temp_high": "Temp-H", "temp_low": "Temp-L",
                 "trans_high": "Trans-H", "trans_low": "Trans-L"}
# pair -> (high variant, low variant); the trans pair runs without heating
PAIRS = {"temp": ("temp_high", "temp_low"), "trans": ("trans_high", "trans_low")}
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
    ("absorption", "Absorption", "dB/cm"),
    ("absorption_fraction", "Absorption fraction", ""),
    ("absorption_coeff", "Absorption coefficient", "dB/(cm·MHz^y)"),
    ("thermal_conductivity", "Thermal conductivity", "W/m/°C"),
    ("specific_heat_capacity", "Specific heat", "J/kg/°C"),
    ("perfusion", "Perfusion", "mL/min/kg"),
]
THERMAL = ("thermal_conductivity", "specific_heat_capacity", "perfusion", "absorption",
           "absorption_fraction", "absorption_coeff")
ABSORPTION = ("absorption_fraction", "absorption_coeff")
PLOTTED = ("sound_speed", "density", "attenuation", "thermal_conductivity",
           "specific_heat_capacity", "perfusion", "absorption")
FALLBACK_MHZ = 0.5          # the ITRUSST property database compares attenuation at 0.5 MHz
OVERRIDES_KEY = "simulation.uncertainty_overrides"


@dataclass
class PropertySets:
    sets: dict[str, dict] = field(default_factory=dict)      # variant -> {tissue: {prop: value}}
    files: dict[str, Path] = field(default_factory=dict)     # variant -> resolved file
    listed: dict[str, dict] = field(default_factory=dict)    # variant -> its bounds layer (file + edits)
    edits: dict[str, dict] = field(default_factory=dict)     # variant -> {tissue: {prop: value}} edited
    shipped: dict[str, dict] = field(default_factory=dict)   # variant -> the shipped bounds, unedited
    errors: dict[str, str] = field(default_factory=dict)     # variant -> why it cannot be used
    freq_mhz: float = FALLBACK_MHZ
    freq_from_transducer: bool = False

    def value(self, variant: str, tissue: str, prop: str) -> Any:
        props = self.sets.get(variant, {}).get(tissue, {})
        if prop == "attenuation":
            a, y = props.get("alpha_coeff"), props.get("alpha_power")
            return attenuation_db_cm(a, y, self.freq_mhz)
        if prop == "absorption":
            # heat deposition at the transducer frequency: the coefficient, else the fraction of the attenuation
            y = props.get("alpha_power")
            if not is_empty(props.get("absorption_coeff")):
                return attenuation_db_cm(props["absorption_coeff"], y, self.freq_mhz)
            att = attenuation_db_cm(props.get("alpha_coeff"), y, self.freq_mhz)
            frac = props.get("absorption_fraction")
            return att * float(frac) if att is not None and isinstance(frac, (int, float)) else None
        value = props.get(prop)
        return None if is_empty(value) else value

    def inherited(self, variant: str, tissue: str, prop: str) -> bool:
        """The variant's value comes from medium_properties, not from its bounds or edits."""
        own = self.listed.get(variant, {}).get(tissue, {})
        if prop == "attenuation":
            return "alpha_coeff" not in own and "alpha_power" not in own
        if prop in ("absorption", *ABSORPTION):
            return not any(not is_empty(own.get(k)) for k in ABSORPTION)
        return prop not in own

    def edited(self, variant: str, tissue: str, prop: str) -> bool:
        return prop in self.edits.get(variant, {}).get(tissue, {})

    def bound(self, variant: str, tissue: str, prop: str) -> Any:
        """The shipped ITRUSST bound (None when the bounds do not list it)."""
        value = self.shipped.get(variant, {}).get(tissue, {}).get(prop)
        return None if is_empty(value) else value

    def n_edits(self) -> int:
        return sum(len(p) for t in self.edits.values() for p in t.values())


def absorption_layer(merged: dict, layer: dict) -> dict:
    """medium_absorption_layer.m: a layer that sets one absorption field of a tissue clears the other."""
    for tissue, props in layer.items():
        if not isinstance(props, dict) or not isinstance(merged.get(tissue), dict):
            continue
        for key, other in (ABSORPTION, ABSORPTION[::-1]):
            if not is_empty(props.get(key)) and other in merged[tissue]:
                merged[tissue][other] = []
    return merged


def attenuation_db_cm(alpha_coeff: Any, alpha_power: Any, f_mhz: float) -> float | None:
    try:
        return float(alpha_coeff) * float(f_mhz) ** float(alpha_power)
    except (TypeError, ValueError):
        return None


def variant_file(prestus_root: str | Path, variant: str) -> Path:
    """The shipped ITRUSST bounds of a variant."""
    return Path(prestus_root) / "config" / "uncertainty" / f"config_medium_{variant}.yaml"


def heating(variant: str) -> bool:
    """uncertainty_variants.m: the Trans variants run acoustic only."""
    return variant not in PAIRS["trans"]


def editable(variant: str, prop: str) -> bool:
    """A key the edits may set: a config property (not a derived row) the variant uses."""
    return prop not in ("attenuation", "absorption") and (heating(variant) or prop not in THERMAL)


def overrides(model) -> dict:
    """simulation.uncertainty_overrides as {variant: {tissue: {prop: value}}} ({} when unset)."""
    value = model.get(OVERRIDES_KEY)
    return value if isinstance(value, dict) else {}


def override_errors(model) -> list[str]:
    """The checks uncertainty_medium_layer.m makes on simulation.uncertainty_overrides."""
    errs = []
    ov = model.get(OVERRIDES_KEY)
    if not is_empty(ov) and not isinstance(ov, dict):
        return ["must map variants to tissues and properties"]
    base = default_properties(model)
    for variant, tissues in overrides(model).items():
        if variant not in VARIANTS:
            errs.append(f"{variant}: not a variant (use {', '.join(VARIANTS)})")
            continue
        if not isinstance(tissues, dict):
            errs.append(f"{variant} must list tissues")
            continue
        for tissue, props in tissues.items():
            if tissue not in base:
                errs.append(f"{variant}.{tissue}: not a tissue of medium_properties")
                continue
            if not isinstance(props, dict):
                errs.append(f"{variant}.{tissue} must list properties")
                continue
            for prop, value in props.items():
                key = f"{variant}.{tissue}.{prop}"
                if prop not in base[tissue]:
                    errs.append(f"{key}: not a medium property")
                elif isinstance(value, bool) or not isinstance(value, (int, float)):
                    errs.append(f"{key} must be a number")
                elif not editable(variant, prop):
                    errs.append(f"{key}: this variant runs acoustic only, so thermal properties are not used")
            if all(not is_empty(props.get(k)) for k in ABSORPTION):
                errs.append(f"{variant}.{tissue} sets both absorption_fraction and absorption_coeff")
    return errs


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


def signature(model, prestus_root: str | Path) -> tuple:
    """Changes when the property sets would: bounds files, edits, medium_properties, frequency."""
    files = []
    for variant in VARIANTS:
        f = variant_file(prestus_root, variant)
        try:
            files.append((str(f), f.stat().st_mtime))
        except OSError:
            files.append((str(f), None))
    base = repr(sorted((t, sorted(p.items())) for t, p in default_properties(model).items()))
    return tuple(files), base, repr(model.get(OVERRIDES_KEY)), transducer_mhz(model)


def property_sets(model, prestus_root: str | Path) -> PropertySets:
    base = default_properties(model)
    f = transducer_mhz(model)
    out = PropertySets(sets={"default": base}, freq_mhz=f or FALLBACK_MHZ, freq_from_transducer=f is not None)
    edits = overrides(model)
    for variant in VARIANTS:
        path = variant_file(prestus_root, variant)
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
        out.shipped[variant] = listed
        own = edits.get(variant)
        if isinstance(own, dict) and own:
            own = {t: p for t, p in own.items() if isinstance(p, dict) and p}
            out.edits[variant] = own
            listed = absorption_layer(mergestruct(listed, own), own)
        out.listed[variant] = listed
        out.sets[variant] = absorption_layer(mergestruct(base, listed), listed)
    return out
