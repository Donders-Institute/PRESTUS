"""Per-user Nexus state (~/.config/prestus_nexus/last_session.yaml)."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import yaml

MAX_RECENT = 10


def session_file() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config"
    return Path(base) / "prestus_nexus" / "last_session.yaml"


class Session:
    """Small key/value store; never raises on unreadable or unwritable files."""

    def __init__(self, path: Path | None = None):
        self.path = path or session_file()
        self.data: dict[str, Any] = {}
        try:
            loaded = yaml.safe_load(self.path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                self.data = loaded
        except (OSError, yaml.YAMLError):
            pass

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    def set(self, key: str, value: Any) -> None:
        self.data[key] = value

    def add_recent(self, study_file: str | Path) -> None:
        path = str(Path(study_file).resolve())
        recent = [p for p in self.data.get("recent", []) if p != path]
        self.data["recent"] = [path] + recent[:MAX_RECENT - 1]

    def recent(self) -> list[str]:
        return [p for p in self.data.get("recent", []) if Path(p).is_file()]

    def save(self) -> None:
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(yaml.safe_dump(self.data, sort_keys=True), encoding="utf-8")
        except OSError:
            pass
