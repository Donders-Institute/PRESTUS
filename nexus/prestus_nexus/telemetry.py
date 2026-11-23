"""PRESTUS usage-statistics consent, as the MATLAB side records it (Tools ▸ Usage statistics…).

The Nexus runs MATLAB with ``matlab -batch``, where PRESTUS never asks for consent
(telemetry_setup.m: non-interactive sessions only print a notice). This module reads and
writes the same record that telemetry_write_consent.m and telemetry_setup_reset.m use:

    ~/.prestus/telemetry.json
    {"opt_in":true,"levels":["basic","performance","outcomes"],"decided_on":"YYYY-MM-DD",
     "prestus_ver":"<hash>","notice":3}

(one line of MATLAB jsonencode output; prestus_ver is the short git hash that
prestus_version() returns first). Consent is per level (LEVELS, as
telemetry_level_info.m); a record without levels predates them and shares what its
notice covered: usage and performance (telemetry_levels.m). No record means undecided,
which sends nothing. The Nexus itself never sends anything; PRESTUS sends events from
MATLAB (track_usage.m) for the shared levels only.
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


# (name, label, description), in prompt order, as telemetry_level_info.m
LEVELS = (
    ("basic", "Usage and errors",
     "Which features, platforms, and versions are used, and whether runs succeed"),
    ("performance", "Performance",
     "Hardware, grid size, run times, and memory: improves run-time estimates"),
    ("outcomes", "Summary outcomes (anonymous)",
     "Rounded summary statistics: brain intensity, temperature rise, targeting error, and skull "
     "properties (no maps, coordinates, or targets)"),
)
LEVEL_NAMES = tuple(name for name, _label, _desc in LEVELS)
LEGACY_LEVELS = ("basic", "performance")    # an opt-in recorded before levels existed


@dataclass
class Consent:
    state: str                  # 'undecided' | 'opted_in' | 'opted_out' | 'unreadable'
    decided_on: str = ""
    prestus_ver: str = ""
    notice: int = 1             # notice version the decision was made under (1 when absent)
    levels: tuple = ()          # shared levels (LEVEL_NAMES order); empty unless opted in

    @property
    def opted_in(self) -> bool:
        return self.state == "opted_in"


def _known_levels(v) -> tuple:
    """The known level names in v, in prompt order (telemetry_levels.m)."""
    if isinstance(v, str):
        v = [v]
    if not isinstance(v, (list, tuple)):
        return ()
    asked = {str(x).strip().lower() for x in v}
    return tuple(n for n in LEVEL_NAMES if n in asked)


def read_consent(path: Path | None = None) -> Consent:
    """The recorded decision. As telemetry_levels.m: only opt_in == true with at least one
    level counts as consent; a record without levels shares LEGACY_LEVELS; an unreadable
    file counts as no consent (and, unlike a missing one, as a decision)."""
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
    n = cfg.get("notice")              # a numeric scalar, else version 1
    notice = int(n) if isinstance(n, (int, float)) and not isinstance(n, bool) else 1
    levels = (_known_levels(cfg["levels"]) if "levels" in cfg else LEGACY_LEVELS) if yes else ()
    state = "opted_in" if levels else "opted_out"
    return Consent(state, str(cfg.get("decided_on") or ""), str(cfg.get("prestus_ver") or ""),
                   notice, levels)


def prestus_hash(prestus_root: Path | None) -> str:
    """Short git hash of the PRESTUS checkout ('unknown' as prestus_version.m)."""
    if prestus_root is None:
        return "unknown"
    kwargs = {"creationflags": getattr(subprocess, "CREATE_NO_WINDOW", 0)} if os.name == "nt" else {}
    try:
        out = subprocess.run(["git", "-C", str(prestus_root), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, encoding="utf-8", errors="replace",
                             timeout=10, check=False, **kwargs)
    except (OSError, subprocess.SubprocessError, ValueError):
        return "unknown"
    return out.stdout.strip() if out.returncode == 0 and out.stdout.strip() else "unknown"


NOTICE_VERSION = 3   # telemetry_notice_version.m: 3 adds consent per level and outcomes
# what each notice version added (telemetry_notice_version.m), shown when consent predates it
NOTICE_ADDS = {2: "the medians of the modelled skull sound speed, density, and attenuation",
               3: "a choice per level, and a new level: summary outcomes (rounded, anonymous)"}


def encode_record(levels, decided_on: str, prestus_ver: str) -> str:
    """The file text as telemetry_write_consent.m writes it (MATLAB jsonencode + newline):
    the shared levels, and the notice version the decision was made under."""
    levels = _known_levels(list(levels))
    cfg = {"opt_in": bool(levels), "levels": list(levels), "decided_on": decided_on,
           "prestus_ver": prestus_ver, "notice": NOTICE_VERSION}
    return json.dumps(cfg, separators=(",", ":"), ensure_ascii=False) + "\n"


def set_consent(choice, prestus_root: Path | None = None, path: Path | None = None,
                today: str | None = None) -> Path:
    """telemetry_set_consent(choice): True shares every level, False none, a list of level
    names those levels. Records the decision; returns the file."""
    if isinstance(choice, bool):
        levels = LEVEL_NAMES if choice else ()
    else:
        levels = _known_levels(list(choice))
    path = path or consent_file()
    path.parent.mkdir(parents=True, exist_ok=True)
    today = today or _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")
    write_text_atomic(path, encode_record(levels, today, prestus_hash(prestus_root)))
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
    """'on' / 'off' when PRESTUS_TELEMETRY is set for runs started from the Nexus, else ''
    (telemetry_env_override.m: on/1/true, or a comma-separated list of levels, is 'on')."""
    v = os.environ.get(ENV, "").strip().lower()
    mode = {"on": "on", "1": "on", "true": "on", "off": "off", "0": "off", "false": "off"}.get(v, "")
    if not mode and v and all(p.strip() in LEVEL_NAMES for p in v.split(",")):
        mode = "on"
    return mode


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
        shared = ", ".join(label.lower() for name, label, _d in LEVELS if name in c.levels)
        if outdated(c):
            return ("Sharing usage statistics · notice updated",
                    f"You share {shared}{since}. The notice has changed since; "
                    "click to read what is new and choose." , "green")
        part = "" if len(c.levels) == len(LEVEL_NAMES) else f" ({len(c.levels)} of {len(LEVEL_NAMES)})"
        return f"Sharing usage statistics{part}", f"You share {shared}{since}." + click, "green"
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
    "Each simulation run can send one anonymous record from MATLAB to the PRESTUS "
    "developers. It guides development, improves run-time estimates, and lets us publish "
    "statistics across sites, such as skull transmission. "
    "Every record holds a random ID, the PRESTUS version, and whether the run succeeded; "
    "each level below adds its fields.")
NEVER_SENT = ("<b>Never stored in the usage statistics:</b> subject IDs, file paths, "
              "coordinates, target names, free text, IP addresses, host names, or any map or "
              "image of the simulation. Like any web service, the hosting provider (Supabase) "
              "keeps a request log with the sender's IP address and approximate location for "
              "1 day; it is not linked to the statistics.")

def _home_short(path: Path) -> str:
    try:
        return "~/" + path.relative_to(Path.home()).as_posix()
    except ValueError:
        return str(path)


def _headline(c: Consent) -> str:
    """The one bold line: the decision as it stands."""
    since = f" (since {c.decided_on})" if c.decided_on else ""
    return {"opted_in": f"You share usage statistics{since}.",
            "opted_out": f"You share no usage statistics{since}.",
            "unreadable": "Your saved choice cannot be read, so nothing is shared.",
            "undecided": "Share anonymous usage statistics?"}[c.state]


def initial_levels(c: Consent) -> tuple:
    """The boxes ticked when the dialog opens: every level for a new choice (undecided, or
    an Allow under an earlier notice), else the recorded levels."""
    if c.state in ("undecided", "unreadable") or outdated(c):
        return LEVEL_NAMES
    return c.levels


def _notes(c: Consent, path: Path) -> str:
    """Lines above the notice: what changed, overrides, and where the decision is kept."""
    out = []
    if outdated(c):
        adds = "; ".join(NOTICE_ADDS[v] for v in sorted(NOTICE_ADDS) if v > c.notice)
        out.append(f"<b>New since you allowed:</b> {adds}. Until you save a choice, runs "
                   "keep sharing what you allowed then.")
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


class ConsentDialog:
    """The consent dialog's widgets and logic, built on demand (PyQt6 imported here):
    "Share all (recommended)" in bold, ticked, with one indented box per level below it;
    unticking every level shares nothing. There is no Decline button."""

    def __init__(self, parent, current: Consent, path: Path, at_launch: bool):
        from PyQt6.QtCore import Qt
        from PyQt6.QtWidgets import (QCheckBox, QDialog, QDialogButtonBox, QLabel,
                                     QVBoxLayout, QWidget)
        self.dialog = d = QDialog(parent)
        d.setWindowTitle("PRESTUS usage statistics")
        d.setMinimumWidth(520)
        v = QVBoxLayout(d)
        head = QLabel(f"<b>{_headline(current)}</b>")
        v.addWidget(head)
        text = QLabel(_notes(current, path) + NOTICE)
        text.setWordWrap(True)
        text.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(text)
        self.all_box = QCheckBox("Share all (recommended)")
        f = self.all_box.font()
        f.setBold(True)
        self.all_box.setFont(f)
        v.addWidget(self.all_box)
        self.boxes = {}
        ticked = initial_levels(current)
        for name, label, desc in LEVELS:
            row = QWidget()
            h = QVBoxLayout(row)
            h.setContentsMargins(24, 0, 0, 2)
            h.setSpacing(0)
            box = QCheckBox(label)
            box.setChecked(name in ticked)
            note = QLabel(desc)
            note.setWordWrap(True)
            note.setStyleSheet("color: gray;")
            note.setContentsMargins(22, 0, 0, 0)
            h.addWidget(box)
            h.addWidget(note)
            v.addWidget(row)
            self.boxes[name] = box
            box.toggled.connect(self._sync_all)
        self.all_box.clicked.connect(self._set_all)
        self._sync_all()
        never = QLabel(NEVER_SENT)
        never.setWordWrap(True)
        never.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(never)
        foot = QLabel(_footer(current, path))
        foot.setOpenExternalLinks(True)
        foot.setTextFormat(Qt.TextFormat.RichText)
        v.addWidget(foot)
        buttons = QDialogButtonBox()
        self.save = buttons.addButton("Save", QDialogButtonBox.ButtonRole.AcceptRole)
        self.forget = None
        if at_launch:
            self.close = buttons.addButton("Ask me later", QDialogButtonBox.ButtonRole.RejectRole)
        else:
            self.forget = buttons.addButton("Forget my decision",
                                            QDialogButtonBox.ButtonRole.ResetRole)
            self.forget.setEnabled(current.state != "undecided")
            self.close = buttons.addButton(QDialogButtonBox.StandardButton.Close)
        self.save.setDefault(True)
        self.clicked = None
        buttons.clicked.connect(self._clicked)
        v.addWidget(buttons)

    def _set_all(self, checked: bool) -> None:
        for box in self.boxes.values():
            box.setChecked(checked)

    def _sync_all(self, *_args) -> None:
        from PyQt6.QtCore import Qt
        n = sum(b.isChecked() for b in self.boxes.values())
        state = (Qt.CheckState.Checked if n == len(self.boxes) else
                 Qt.CheckState.Unchecked if n == 0 else Qt.CheckState.PartiallyChecked)
        self.all_box.blockSignals(True)
        self.all_box.setTristate(state == Qt.CheckState.PartiallyChecked)
        self.all_box.setCheckState(state)
        self.all_box.blockSignals(False)

    def levels(self) -> tuple:
        return tuple(n for n in LEVEL_NAMES if self.boxes[n].isChecked())

    def _clicked(self, button) -> None:
        self.clicked = button
        self.dialog.accept() if button is self.save else self.dialog.reject()


def show_consent_dialog(parent=None, prestus_root: Path | None = None,
                        at_launch: bool = False) -> Consent:
    """Tools ▸ Usage statistics… (and the launch prompt, at_launch): show the notice with
    the levels; Save records the ticked levels (none ticked: nothing is shared), Forget my
    decision removes the record. Closing ('Ask me later' at launch) changes nothing."""
    path = consent_file()
    ui = ConsentDialog(parent, read_consent(path), path, at_launch)
    ui.dialog.exec()
    try:
        if ui.clicked is ui.save:
            set_consent(list(ui.levels()), prestus_root, path)
        elif ui.forget is not None and ui.clicked is ui.forget:
            reset_consent(path)
    except OSError as exc:
        from .ui.errors import show_error
        show_error(parent, "Usage statistics", f"Cannot write {path}: {exc.strerror or exc}")
    return read_consent(path)
