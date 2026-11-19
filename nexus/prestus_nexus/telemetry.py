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
    notice: int = 1             # notice version the decision was made under (1 when absent)

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
    n = cfg.get("notice")              # track_usage.m: a numeric scalar, else version 1
    notice = int(n) if isinstance(n, (int, float)) and not isinstance(n, bool) else 1
    return Consent(state, str(cfg.get("decided_on") or ""), str(cfg.get("prestus_ver") or ""),
                   notice)


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
# what each notice version added (telemetry_notice_version.m), shown when consent predates it
NOTICE_ADDS = {2: "the medians of the modelled skull sound speed, density, and attenuation"}


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


def outdated(c: Consent) -> bool:
    """Allowed under an earlier notice: fields the later notices add are not sent."""
    return c.opted_in and c.notice < NOTICE_VERSION


def status(c: Consent, env: str | None = None) -> tuple[str, str, str]:
    """(short text, tooltip, colour role) for the status-bar indicator: 'green' when runs
    started from this Nexus send events, 'amber' when they do not."""
    env = env_override() if env is None else env
    since = f" since {c.decided_on}" if c.decided_on else ""
    click = " Click to change."
    if env == "off":
        return ("Not sharing usage statistics",
                f"{ENV}=off is set, so runs started here send nothing." + click, "amber")
    if c.state == "opted_in":
        if outdated(c):
            return ("Sharing usage statistics · notice updated",
                    f"You allowed{since}. The notice has changed since; "
                    "click to read what is new and confirm." , "green")
        return "Sharing usage statistics", f"You allowed{since}." + click, "green"
    if c.state == "undecided" and env == "on":
        return ("Sharing usage statistics (this session)",
                f"{ENV}=on is set; no decision is recorded." + click, "green")
    if c.state == "undecided":
        return ("Not sharing usage statistics · choose", "You have not decided yet; "
                "nothing is sent until you allow it." + click, "amber")
    if c.state == "unreadable":
        return ("Not sharing usage statistics",
                "The saved decision cannot be read, so nothing is sent." + click, "amber")
    return "Not sharing usage statistics", f"You declined{since}." + click, "amber"


def should_prompt(c: Consent, env: str | None = None, asked_notice: int = 0) -> bool:
    """Ask at launch? As telemetry_setup.m with a GUI: no decision recorded and no
    PRESTUS_TELEMETRY set (asked at every launch until answered: 'Ask me later'). In
    addition, an Allow given under an earlier notice is asked again, once per notice
    version (asked_notice: the version last asked about, remembered in the session)."""
    env = env_override() if env is None else env
    if env == "off":
        return False
    if c.state == "undecided":
        return env != "on"
    return outdated(c) and asked_notice < NOTICE_VERSION


NOTICE = (
    "With your permission, each simulation run sends one anonymous record from MATLAB, so "
    "the developers can see which features and platforms PRESTUS is used on. The Nexus "
    "itself sends nothing.<br><br>"
    "<b>Sent:</b> software versions and hardware; which modules, transducer, and medium "
    "were used; grid size and run times; skull property medians; whether the run "
    "succeeded; a random ID.<br>"
    "<b>Never sent:</b> subject IDs, file paths, coordinates, free text, IP addresses, "
    "host names, or simulation results.")

def _home_short(path: Path) -> str:
    try:
        return "~/" + path.relative_to(Path.home()).as_posix()
    except ValueError:
        return str(path)


def _headline(c: Consent) -> str:
    """The one bold line: the decision as it stands."""
    since = f" (since {c.decided_on})" if c.decided_on else ""
    return {"opted_in": f"You share usage statistics{since}.",
            "opted_out": f"You do not share usage statistics{since}.",
            "unreadable": "Your saved decision cannot be read, so nothing is shared.",
            "undecided": "Share anonymous usage statistics?"}[c.state]


def _notes(c: Consent, path: Path) -> str:
    """Lines above the notice: what changed, overrides, and where the decision is kept."""
    out = []
    if outdated(c):
        adds = "; ".join(NOTICE_ADDS[v] for v in sorted(NOTICE_ADDS) if v > c.notice)
        out.append(f"<b>New since you allowed:</b> {adds}. These are sent only once you "
                   "allow again.")
    env = env_override()
    if env == "off":
        out.append(f"{ENV}=off is set, so runs started here send nothing.")
    elif env == "on" and c.state == "undecided":
        out.append(f"{ENV}=on is set, so runs started here send records.")
    return "".join(f"{line}<br><br>" for line in out)


def _footer(c: Consent, path: Path) -> str:
    kept = "" if c.state == "undecided" else f"Saved in {_home_short(path)}. "
    return (f"<br><br><span style='color:gray'>{kept}"
            f"<a href='{DOC_URL}'>Full details</a></span>")


def show_consent_dialog(parent=None, prestus_root: Path | None = None,
                        at_launch: bool = False) -> Consent:
    """Tools ▸ Usage statistics… (and the launch prompt, at_launch): show the notice and
    the recorded decision; record what the user picks (Allow / Decline / Forget my
    decision). Closing ('Ask me later' at launch) changes nothing."""
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QLabel, QMessageBox
    path = consent_file()
    current = read_consent(path)
    box = QMessageBox(parent)
    box.setWindowTitle("PRESTUS usage statistics")
    box.setIcon(QMessageBox.Icon.Question)
    box.setTextFormat(Qt.TextFormat.RichText)
    box.setText(f"<b>{_headline(current)}</b>")
    box.setInformativeText(_notes(current, path) + NOTICE + _footer(current, path))
    for label in box.findChildren(QLabel):
        label.setOpenExternalLinks(True)
    allow = box.addButton("Allow", QMessageBox.ButtonRole.AcceptRole)
    decline = box.addButton("Decline", QMessageBox.ButtonRole.RejectRole)
    forget = None
    if not at_launch:
        forget = box.addButton("Forget my decision", QMessageBox.ButtonRole.ResetRole)
        forget.setEnabled(current.state != "undecided")
    if at_launch:
        close = box.addButton("Ask me later", QMessageBox.ButtonRole.RejectRole)
    else:
        close = box.addButton(QMessageBox.StandardButton.Close)
    box.setDefaultButton(close)
    box.setEscapeButton(close)
    box.exec()
    clicked = box.clickedButton()
    try:
        if clicked is allow:
            set_consent(True, prestus_root, path)
        elif clicked is decline:
            set_consent(False, prestus_root, path)
        elif forget is not None and clicked is forget:
            reset_consent(path)
    except OSError as exc:
        from .ui.errors import show_error
        show_error(parent, "Usage statistics", f"Cannot write {path}: {exc.strerror or exc}")
    return read_consent(path)
