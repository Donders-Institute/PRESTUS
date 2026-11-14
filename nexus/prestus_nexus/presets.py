"""Named hpc.* override bundles (run_presets.yaml) for the Run tab's "Run profile"
picker. Applying a preset only sets keys on the open study via ConfigModel.set -- it
never writes a new file and never submits anything; see launcher.py for the (manual)
local-run path and doc_hpc.md for the Donders/Snellius scheduler workflows a preset
is meant to prepare a study for."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import yamlio
from .model import ConfigModel, values_equal

_PRESETS_FILE = Path(__file__).parent / "run_presets.yaml"


@dataclass(frozen=True)
class Preset:
    name: str
    description: str
    overrides: dict[str, object] = field(default_factory=dict)   # dotted key -> value


def load_presets(path: str | Path | None = None) -> list[Preset]:
    data = yamlio.load_file(path or _PRESETS_FILE) or {}
    out = []
    for p in data.get("presets") or []:
        out.append(Preset(name=str(p.get("name", "")), description=str(p.get("description", "")),
                          overrides=dict(p.get("overrides") or {})))
    return out


def apply(model: ConfigModel, preset: Preset) -> None:
    """Set every hpc.* override, the same way any other field write does (so the ●
    'changed' marker and Issues panel pick it up for free). Keys outside hpc.* are
    rejected -- a preset only ever touches hpc.*, and this guards against a typo in
    run_presets.yaml silently writing somewhere else."""
    for key in preset.overrides:
        if not _allowed(key):
            raise ValueError(f"run preset {preset.name!r}: override key {key!r} is not under hpc.* "
                             "or platform")
    for key, value in preset.overrides.items():
        model.set(key, value)


def _allowed(key: str) -> bool:
    return key in ("hpc", "platform") or key.startswith("hpc.")


def matches(model: ConfigModel, preset: Preset) -> bool:
    """The study already has every value the preset sets (and the preset sets some)."""
    return bool(preset.overrides) and all(
        values_equal(model.get(k), v) for k, v in preset.overrides.items())


def current(model: ConfigModel, all_presets: list[Preset]) -> Preset | None:
    """The first preset whose values the study has, if any."""
    return next((p for p in all_presets if matches(model, p)), None)
