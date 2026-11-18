"""PRESTUS usage-statistics consent, as the MATLAB side records it (Tools ▸ Usage statistics…).

The Nexus runs MATLAB with ``matlab -batch``, where PRESTUS never asks for consent
(telemetry_setup.m: non-interactive sessions only print a notice). This module reads and
writes the same record that telemetry_set_consent.m and telemetry_setup_reset.m use:

    ~/.prestus/telemetry.json   {"opt_in":true,"decided_on":"YYYY-MM-DD","prestus_ver":"<hash>"}

(one line of MATLAB jsonencode output; prestus_ver is the short git hash that
prestus_version() returns first). No record means undecided, which sends nothing. The
Nexus itself never sends anything; PRESTUS sends events from MATLAB (track_usage.m) only
when the record says opt_in: true.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .fsio import write_text_atomic

DOC_URL = "https://github.com/Donders-Institute/PRESTUS/blob/main/doc/doc_telemetry.md"
ENV = "PRESTUS_TELEMETRY"


def prefdir() -> Path:
    """~/.prestus (prefdir_prestus.m: Java's user.home, i.e. the home folder)."""
    return Path.home() / ".prestus"


def consent_file() -> Path:
    return prefdir() / "telemetry.json"


@dataclass
class Consent:
    state: str                  # 'undecided' | 'opted_in' | 'opted_out' | 'unreadable'
    decided_on: str = ""
    prestus_ver: str = ""

    @property
    def opted_in(self) -> bool:
        return self.state == "opted_in"


def read_consent(path: Path | None = None) -> Consent:
    """The recorded decision. As track_usage.m: only opt_in == true counts as consent; an
    unreadable file counts as no consent (and, unlike a missing one, as a decision)."""
    path = path or consent_file()
    if not path.is_file():
        return Consent("undecided")
    try:
        cfg = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return Consent("unreadable")
    if not isinstance(cfg, dict):
        return Consent("unreadable")
    v = cfg.get("opt_in")              # isequal(cfg.opt_in, true): true, or the number 1
    yes = v is True or (isinstance(v, (int, float)) and not isinstance(v, bool) and v == 1)
    state = "opted_in" if yes else "opted_out"
    return Consent(state, str(cfg.get("decided_on") or ""), str(cfg.get("prestus_ver") or ""))


def prestus_hash(prestus_root: Path | None) -> str:
    """Short git hash of the PRESTUS checkout ('unknown' as prestus_version.m)."""
    if prestus_root is None:
        return "unknown"
    kwargs = {"creationflags": getattr(subprocess, "CREATE_NO_WINDOW", 0)} if os.name == "nt" else {}
    try:
        out = subprocess.run(["git", "-C", str(prestus_root), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, encoding="utf-8", errors="replace",
                             timeout=10, check=False, **kwargs)
    except Exception:
        return "unknown"
    return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else "unknown"


NOTICE_VERSION = 2   # telemetry_notice_version.m: 2 adds the skull property medians (2026-11-17)


def encode_record(opt_in: bool, decided_on: str, prestus_ver: str) -> str:
    """The file text as MATLAB's jsonencode(struct) + newline writes it (with the notice
    version the decision was made under, as telemetry_set_consent.m)."""
    cfg = {"opt_in": bool(opt_in), "decided_on": decided_on, "prestus_ver": prestus_ver,
           "notice": NOTICE_VERSION}
    return json.dumps(cfg, separators=(",", ":"), ensure_ascii=False) + "\n"


def set_consent(opt_in: bool, prestus_root: Path | None = None, path: Path | None = None,
                today: str | None = None) -> Path:
    """telemetry_set_consent(opt_in): record the decision; returns the file."""
    path = path or consent_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    today = today or _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")
    write_text_atomic(path, encode_record(opt_in, today, prestus_hash(prestus_root)))
    return path


def reset_consent(path: Path | None = None) -> bool:
    """telemetry_setup_reset(): forget the decision; True if a record was removed."""
    path = path or consent_file()
    try:
        path.unlink()
        return True
    except FileNotFoundError:
        return False


def env_override() -> str:
    """'on' / 'off' when PRESTUS_TELEMETRY is set for runs started from the Nexus, else ''."""
    v = os.environ.get(ENV, "").strip().lower()
    return {"on": "on", "1": "on", "true": "on", "off": "off", "0": "off", "false": "off"}.get(v, "")


NOTICE = (
    "PRESTUS can optionally collect anonymous usage statistics to help the developers "
    "understand which features are used and on which platforms. One event is sent per "
    "pipeline run, from MATLAB; the Nexus itself sends nothing.\n\n"
    "Collected (examples): PRESTUS, MATLAB, and k-Wave versions, OS platform; CPU model, core "
    "count, RAM, GPU model, installed toolboxes, HPC name; which pipeline modules are enabled; "
    "transducer type, frequency, and element count; simulation medium, layer names, whether "
    "pCT is used; grid resolution, grid size, time steps, simulation stage durations and peak "
    "RAM; the median sound speed, density, and attenuation of the modelled skull (all skull "
    "voxels, and on the beam path); requested job resources, queue wait, and solver settings; "
    "run start, duration, success or failure, the error type and failing stage (no message "
    "text); a random local ID.\n\n"
    "Never collected: subject IDs, file paths, coordinates, or any free text; IP addresses or "
    "host names; pressure values, simulation results, or tissue property maps.")


def _state_text(c: Consent, path: Path) -> str:
    if c.state == "opted_in":
        text = f"Current decision: <b>allowed</b> (since {c.decided_on or '?'})."
    elif c.state == "opted_out":
        text = f"Current decision: <b>declined</b> (since {c.decided_on or '?'})."
    elif c.state == "unreadable":
        text = "Current decision: the record cannot be read, so nothing is sent."
    else:
        text = "No decision recorded: nothing is sent."
    env = env_override()
    if env == "off":
        text += f"<br>{ENV}=off is set: runs started from this Nexus send nothing."
    elif env == "on" and c.state == "undecided":
        text += f"<br>{ENV}=on is set: runs started from this Nexus send events."
    return text + f"<br>Stored in <tt>{path}</tt>."


def show_consent_dialog(parent=None, prestus_root: Path | None = None) -> Consent:
    """Tools ▸ Usage statistics…: show the notice and the recorded decision; record what
    the user picks (Allow / Decline / Forget my decision). Closing changes nothing."""
    from PyQt6.QtWidgets import QMessageBox
    path = consent_file()
    current = read_consent(path)
    box = QMessageBox(parent)
    box.setWindowTitle("PRESTUS usage statistics")
    box.setIcon(QMessageBox.Icon.Question)
    box.setText(_state_text(current, path))
    box.setInformativeText(NOTICE + f"\n\nFull details: {DOC_URL}")
    allow = box.addButton("Allow", QMessageBox.ButtonRole.AcceptRole)
    decline = box.addButton("Decline", QMessageBox.ButtonRole.RejectRole)
    forget = box.addButton("Forget my decision", QMessageBox.ButtonRole.ResetRole)
    close = box.addButton(QMessageBox.StandardButton.Close)
    forget.setEnabled(current.state != "undecided")
    box.setDefaultButton(close)
    box.setEscapeButton(close)
    box.exec()
    clicked = box.clickedButton()
    try:
        if clicked is allow:
            set_consent(True, prestus_root, path)
        elif clicked is decline:
            set_consent(False, prestus_root, path)
        elif clicked is forget:
            reset_consent(path)
    except OSError as exc:
        from .ui.errors import show_error
        show_error(parent, "Usage statistics", f"Cannot write {path}: {exc.strerror or exc}")
    return read_consent(path)
