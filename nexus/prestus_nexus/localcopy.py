"""Run the Nexus from a local copy of its own code.

PRESTUS often lives on a network drive (e.g. /project mounted over SMB). Python then
reads every module of this package over the network, and the Nexus imports most of
them while it builds its window: on such a mount that alone took ~30 s of a ~45 s
start-up. ``local_main`` mirrors the package's files (.py, .yaml) into the user's
cache folder, copying only files whose size or modification time changed, and runs
the Nexus from there. The checkout stays the only source: any edit there is copied
on the next start.

The original folder is passed on in $PRESTUS_NEXUS_SOURCE (see ``paths.package_dir``),
so PRESTUS is still found next to it. Set PRESTUS_NEXUS_LOCAL_COPY=0 to run from the
checkout directly (e.g. while developing the Nexus itself).

Uses only the standard library: it runs before anything else of the package is
imported.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Callable

SOURCE_ENV = "PRESTUS_NEXUS_SOURCE"
OPT_OUT_ENV = "PRESTUS_NEXUS_LOCAL_COPY"
PACKAGE = "prestus_nexus"
SUFFIXES = (".py", ".yaml")
MANIFEST = "manifest.json"


def cache_root() -> Path:
    """The Nexus's cache folder ($PRESTUS_NEXUS_CACHE, else the platform's user cache)."""
    if os.environ.get("PRESTUS_NEXUS_CACHE"):
        return Path(os.environ["PRESTUS_NEXUS_CACHE"])
    if sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    elif os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache")
    return base / "prestus_nexus"


def _signature(src: Path) -> dict[str, list[int]]:
    """relative path -> [size, mtime_ns] of the package files (one directory listing each)."""
    sig: dict[str, list[int]] = {}
    for folder in (src, src / "ui"):
        with os.scandir(folder) as it:
            for entry in it:
                if entry.is_file() and entry.name.endswith(SUFFIXES):
                    st = entry.stat()
                    rel = entry.name if folder == src else f"ui/{entry.name}"
                    sig[rel] = [st.st_size, st.st_mtime_ns]
    return sig


def mirror(src: Path, dest_root: Path) -> Path:
    """Bring ``dest_root/prestus_nexus`` in line with ``src``; returns ``dest_root``."""
    sig = _signature(src)
    dest = dest_root / PACKAGE
    manifest = dest_root / MANIFEST
    try:
        old = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        old = {}
    if old == sig and all((dest / rel).is_file() for rel in sig):
        return dest_root
    (dest / "ui").mkdir(parents=True, exist_ok=True)
    for rel, stamp in sig.items():
        if old.get(rel) != stamp or not (dest / rel).is_file():
            shutil.copy2(src / rel, dest / rel)
    for rel in set(old) - set(sig):
        (dest / rel).unlink(missing_ok=True)
    manifest.write_text(json.dumps(sig), encoding="utf-8")    # last: a broken copy is redone
    return dest_root


def local_main() -> Callable[[list[str]], int] | None:
    """``app.main`` of the local copy, or None to go on from the checkout: when already
    running from a copy, in a frozen app, when opted out, or if the copy fails."""
    here = Path(__file__).resolve().parent
    if (getattr(sys, "frozen", False) or os.environ.get(SOURCE_ENV)
            or os.environ.get(OPT_OUT_ENV, "1").strip().lower() in ("0", "false", "no", "off")):
        return None
    try:
        dest_root = cache_root() / ("code-" + hashlib.sha1(str(here).encode()).hexdigest()[:10])
        if here.is_relative_to(dest_root):
            return None
        mirror(here, dest_root)
    except OSError as exc:
        print(f"PRESTUS Nexus: running from {here} (no local copy: {exc})", file=sys.stderr)
        return None
    os.environ[SOURCE_ENV] = str(here)
    for name in [n for n in sys.modules if n == PACKAGE or n.startswith(PACKAGE + ".")]:
        del sys.modules[name]
    sys.path.insert(0, str(dest_root))
    importlib.invalidate_caches()
    return importlib.import_module(f"{PACKAGE}.app").main
