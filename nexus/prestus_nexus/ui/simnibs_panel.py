"""SimNIBS environment, version and existing segmentation, on the Segmentation card.

The environment is startup.simnibs_bin_path (also on Setup › Environment); its
version follows from it and is shown read-only. Segmenting again is the card's
'Segment again…' button (segmentation_panel); io.overwrite_simnibs stays on Setup ›
Input / output and is flagged here when on. The custom charm settings
file (segmentation.settings_file) starts as a copy of the install's charm.ini.

When runs are submitted over SSH (Run tab), *Scan cluster* lists the installs on the
cluster (remote.simnibs_scan_command); they are kept per host in the session.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from PyQt6.QtCore import QProcess, Qt, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (QComboBox, QFileDialog, QGridLayout, QHBoxLayout, QLabel,
                             QMessageBox, QPushButton, QSizePolicy, QToolButton, QVBoxLayout, QWidget)

from .. import remote as rmt
from .. import simnibs as sn
from ..model import ConfigModel, is_empty
from ..paths import m2m_dir
from . import icons, style

BIN = "startup.simnibs_bin_path"
RERUN = "io.overwrite_simnibs"
SETTINGS = "segmentation.settings_file"
SESSION_KEY = "simnibs_remote_installs"     # {host: [Install.to_dict()]}


def load_remote_installs(session) -> None:
    """Installs found on cluster hosts in earlier sessions."""
    for host, rows in (session.get(SESSION_KEY, {}) or {}).items():
        try:
            sn.set_remote_installs(host, [sn.Install.from_dict(r) for r in rows])
        except (KeyError, TypeError):  # best-effort: session rows in an older format are dropped
            continue


class SimnibsPanel(QWidget):
    changed = pyqtSignal()

    def __init__(self, model: ConfigModel, label_width: int = 290, session=None,
                 parent: QWidget | None = None):
        super().__init__(parent)
        self.model = model
        self.session = session              # SSH settings (Run tab) and the cluster installs found
        self._proc: QProcess | None = None
        if session is not None:
            load_remote_installs(session)
        self.page = None                    # the card: its 'advanced' switch shows the settings row
        self._installs: list[sn.Install] | None = None    # scanned on first show
        self._updating = False

        g = QGridLayout(self)
        g.setContentsMargins(0, 0, 0, 8)
        g.setHorizontalSpacing(14)
        g.setVerticalSpacing(8)
        g.setColumnMinimumWidth(0, label_width)
        g.setProperty("label_inset", 290 - label_width)   # follows the card (section.mark_label_grid)
        g.setColumnStretch(1, 1)

        def label(text: str, tip: str) -> QLabel:
            lab = QLabel(text)
            lab.setObjectName("FieldLabel")
            lab.setToolTip(tip)
            return lab

        def row(widgets) -> QWidget:
            w = QWidget()
            h = QHBoxLayout(w)
            h.setContentsMargins(0, 0, 0, 0)
            for wd, stretch in widgets:
                h.addWidget(wd, stretch)
            return w

        env_tip = ("The SimNIBS 4 install whose charm segments the T1 (startup.simnibs_bin_path, also on "
                   "Setup › Environment). Found on PATH, in conda environments, and in the default install "
                   "folders; '…' picks another bin folder.")
        self.env = QComboBox()
        self.env.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.env.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.env.setToolTip(env_tip)
        self.env.activated.connect(self._env_chosen)
        scan = QToolButton(text="Scan")
        scan.setToolTip("Look for SimNIBS installs again")
        scan.clicked.connect(self.rescan)
        browse = icons.PathButton("dir", "Choose a SimNIBS bin folder (the one that holds charm)")
        browse.clicked.connect(self._browse)
        self.scan_remote = QToolButton(text="Scan cluster")
        self.scan_remote.setToolTip("List the SimNIBS installs on the cluster runs are submitted to (Run tab, "
                                    "over SSH): conda environments, ~/SimNIBS*, /opt/simnibs*, and the "
                                    "'simnibs' environment modules")
        self.scan_remote.clicked.connect(self.scan_cluster)
        self.scan_note = QLabel()
        self.scan_note.setObjectName("CardSub")
        self.scan_note.setWordWrap(True)
        g.addWidget(label("SimNIBS environment", env_tip), 0, 0)
        env_box = QWidget()
        ev = QVBoxLayout(env_box)
        ev.setContentsMargins(0, 0, 0, 0)
        ev.setSpacing(2)
        ev.addWidget(row([(self.env, 1), (scan, 0), (self.scan_remote, 0), (browse, 0)]))
        ev.addWidget(self.scan_note)
        g.addWidget(env_box, 0, 1)

        ver_tip = "Set by the SimNIBS environment. PRESTUS segments with charm from SimNIBS 4 or later."
        self.version = QComboBox()
        self.version.setEnabled(False)          # follows the environment
        self.version.setToolTip(ver_tip)
        self.version_note = QLabel()
        self.version_note.setObjectName("CardSub")
        self.version_note.setWordWrap(True)
        self.linked = QLabel()
        self.linked.setObjectName("CardSub")
        self.linked.setWordWrap(True)
        self.linked.setTextFormat(Qt.TextFormat.RichText)
        self.linked.setToolTip("Settings that follow the SimNIBS environment while empty: "
                               "hpc.ld_library_path (the install's libraries for charm and the MNI "
                               "conversion) and placement.plantus.env_path (the Python PlanTUS runs in).")
        ver_box = QWidget()
        vv = QVBoxLayout(ver_box)
        vv.setContentsMargins(0, 0, 0, 0)
        vv.setSpacing(2)
        vv.addWidget(row([(self.version, 0), (self.version_note, 1)]))
        vv.addWidget(self.linked)
        g.addWidget(label("SimNIBS version", ver_tip), 1, 0)
        g.addWidget(ver_box, 1, 1)

        seg_tip = ("The segmentation in <path.seg>/m2m_sub-XXX for this config's subject: the SimNIBS "
                   "version and charm settings that made it, compared with the ones selected here.")
        self.state = QLabel()
        self.state.setWordWrap(True)
        self.state.setTextFormat(Qt.TextFormat.RichText)
        self.report = QPushButton("Open charm report")
        self.report.setToolTip("charm_report.html in the m2m folder: the segmentation and registration checks")
        self.report.clicked.connect(self._open_report)
        g.addWidget(label("Existing segmentation", seg_tip), 2, 0)
        g.addWidget(row([(self.state, 1), (self.report, 0)]), 2, 1)

        set_tip = ("Advanced: change charm's settings (e.g. denoising, atlas initialisation, threads). "
                   "Starts from a copy of the install's charm.ini, which you then edit; charm reads it "
                   "with --usesettings. Report the settings you changed.")
        self.custom = QPushButton("Customise charm settings…")
        self.custom.setToolTip(set_tip)
        self.custom.clicked.connect(self._customise)
        self.edit_settings = QPushButton("Edit settings file")
        self.edit_settings.clicked.connect(self._edit_settings)
        self.settings_label = label("charm settings", set_tip)
        self.settings_row = row([(self.custom, 0), (self.edit_settings, 0), (QWidget(), 1)])
        g.addWidget(self.settings_label, 3, 0)
        g.addWidget(self.settings_row, 3, 1)

    # ── state ───────────────────────────────────────────────────────────
    def _host(self) -> rmt.RemoteHost | None:
        return rmt.host_from_session(self.session) if self.session is not None else None

    def installs(self) -> list[sn.Install]:
        """Installs on this computer, then those found on the cluster runs are submitted to."""
        if self._installs is None:
            self._installs = sn.find_installs()
        host = self._host()
        return self._installs + (sn.remote_installs(host.host) if host else [])

    def rescan(self) -> None:
        self._installs = None
        self.refresh()

    def scan_cluster(self) -> None:
        host = self._host()
        if host is None or (self._proc is not None and self._proc.state() != QProcess.ProcessState.NotRunning):
            return
        argv = rmt.simnibs_scan_command(host)
        self._proc = QProcess(self)
        self._proc.setProcessChannelMode(QProcess.ProcessChannelMode.SeparateChannels)
        self._proc.finished.connect(lambda code, _s, h=host.host: self._scan_done(h, code))
        self._proc.errorOccurred.connect(lambda _e, h=host.host: self._scan_done(h, -1))
        self.scan_note.setText(f"Scanning {host.host}…")
        self.scan_remote.setEnabled(False)
        self._proc.start(argv[0], argv[1:])

    def _scan_done(self, host: str, code: int) -> None:
        proc, self._proc = self._proc, None
        self.scan_remote.setEnabled(True)
        if proc is None:
            return
        out = bytes(proc.readAllStandardOutput()).decode(errors="replace")
        found = rmt.parse_simnibs_scan(out, host)
        if found is None:
            err = bytes(proc.readAllStandardError()).decode(errors="replace").strip().splitlines()
            self.scan_note.setText(f"Scan of {host} failed{': ' + err[-1] if err else ''}. "
                                   "Log in with Connect… on the Run tab, then scan again.")
            return
        sn.set_remote_installs(host, found)
        if self.session is not None:
            saved = dict(self.session.get(SESSION_KEY, {}) or {})
            saved[host] = [i.to_dict() for i in found]
            self.session.set(SESSION_KEY, saved)
        self.scan_note.setText(f"{len(found)} SimNIBS install{'s' if len(found) != 1 else ''} found on {host}.")
        self.refresh()
        self.changed.emit()                 # the Issues panel checks the selection again

    def _selected(self) -> sn.Install | None:
        b = self.model.get(BIN)
        return None if is_empty(b) else sn.lookup(b)

    def _m2m(self) -> Path | None:
        sid = self.model.get("subject_id")
        seg = self.model.get("path.seg")
        root = seg if not is_empty(seg) else self.model.get("path.anat")
        if is_empty(root) or isinstance(sid, bool) or not isinstance(sid, (int, float)) \
                or float(sid) != int(sid) or sid < 0:
            return None
        return m2m_dir(root, int(sid))

    # ── actions ─────────────────────────────────────────────────────────
    def _set(self, path: str, value) -> None:
        self.model.set(path, value)
        self.refresh()
        self.changed.emit()

    def _env_chosen(self, index: int) -> None:
        if self._updating:
            return
        b = self.env.itemData(index)
        if b is not None and str(b) != str(self.model.get(BIN) or ""):
            self._choose(str(b))

    def _browse(self) -> None:
        start = str(self.model.get(BIN) or Path.home())
        d = QFileDialog.getExistingDirectory(self, "SimNIBS bin folder (holds charm)", start)
        if not d:
            return
        if not sn.is_simnibs_bin(d):
            QMessageBox.warning(self, "SimNIBS", f"No charm in {d}. Choose the bin folder of a SimNIBS 4 install.")
            return
        self._choose(d)

    def _choose(self, bin_path: str) -> None:
        """Select an environment; linked settings that pointed into the previous install are
        emptied, so they follow the new one (sn.LINKED_KEYS)."""
        old = self._selected()
        self.model.set(BIN, bin_path)
        cleared = []
        for key in sn.LINKED_KEYS:
            v = self.model.get(key)
            if old is not None and isinstance(v, str) and v and old.owns(v):
                self.model.set(key, "")
                cleared.append(key)
        if cleared:
            self.scan_note.setText("Now following the new environment: " + ", ".join(cleared) + ".")
        self.refresh()
        self.changed.emit()

    def _open_report(self) -> None:
        m2m = self._m2m()
        st = sn.segmentation_state(m2m) if m2m else None
        if st and st.report:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(st.report)))

    def _default_ini(self, inst: sn.Install | None) -> str | None:
        """The text of the install's charm.ini (fetched over SSH for a cluster install)."""
        if inst is None:
            return None
        if inst.default_settings is not None:
            return inst.default_settings.read_text(errors="replace")
        host = self._host()
        if inst.host and host is not None and host.host == inst.host and inst.package is not None:
            try:
                r = subprocess.run(rmt.cat_command(host, f"{inst.package}/charm.ini"),
                                   capture_output=True, text=True, timeout=60)
            except (OSError, subprocess.TimeoutExpired):
                return None
            return r.stdout if r.returncode == 0 and r.stdout.strip() else None
        return None

    def _customise(self) -> None:
        inst = self._selected()
        text = self._default_ini(inst)
        if text is None:
            QMessageBox.warning(self, "charm settings",
                                "Select a SimNIBS environment first: its charm.ini is the starting point."
                                + (" (A cluster install's file is read over SSH: Connect… on the Run tab.)"
                                   if inst is not None and inst.host else ""))
            return
        folder = self.model.path.parent if self.model.path else Path.home()
        target, _ = QFileDialog.getSaveFileName(self, "Save charm settings as", str(folder / "charm_settings.ini"),
                                                "charm settings (*.ini)")
        if not target:
            return
        m2m = self._m2m()
        if m2m is not None and Path(target).resolve().parent == m2m.resolve():
            QMessageBox.warning(self, "charm settings",
                                "Save it outside the m2m folder: charm deletes that folder when it re-runs.")
            return
        Path(target).write_text(text)
        self._set(SETTINGS, target)
        self._edit_settings()

    def _edit_settings(self) -> None:
        s = self.model.get(SETTINGS)
        if not is_empty(s) and Path(s).is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(s)))

    # ── display ─────────────────────────────────────────────────────────
    def _linked_text(self, inst: sn.Install | None, c: dict) -> str:
        """The linked settings: following the environment (empty) or set, and to what."""
        if inst is None:
            return ""
        def short(path: str) -> str:
            parts = str(path).replace("\\", "/").rstrip("/").split("/")
            return path if len(parts) <= 3 else "…/" + "/".join(parts[-3:])

        wrong = dict(sn.linked_mismatches(self.model.get, inst))
        parts, full = [], []
        for key, name, derived in (("hpc.ld_library_path", "Libraries", inst.ld_library_path),
                                   ("placement.plantus.env_path", "PlanTUS Python", inst.python)):
            v = self.model.get(key)
            if isinstance(v, str) and v:
                col = c["amber"] if key in wrong else c["ink3"]
                note = " (another SimNIBS install)" if key in wrong else ""
                parts.append(f"<span style='color:{col}'>{name}: set to {short(v)}{note}</span>")
                full.append(f"{key}: {v}")
            else:
                parts.append(f"{name}: {'follow the environment' if derived else 'none found'}")
                full.append(f"{key} (empty) → {derived or 'none found'}")
        self.linked.setToolTip("Settings that follow the SimNIBS environment while empty:\n" + "\n".join(full))
        return "<br>".join(parts)                    # one linked setting per line

    def refresh(self) -> None:
        c = style.active()
        self._updating = True
        try:
            current = self.model.get(BIN)
            inst = self._selected()
            host = self._host()
            self.scan_remote.setVisible(host is not None)
            if host is None:
                self.scan_note.setText("")
            self.env.clear()
            listed = False
            for i in self.installs():
                self.env.addItem(i.label() + ("" if i.supported is not False else "  (not supported)"), str(i.bin))
                if not is_empty(current) and str(current).rstrip("/") == str(i.bin):
                    self.env.setCurrentIndex(self.env.count() - 1)
                    listed = True
            if not is_empty(current) and not listed:
                self.env.insertItem(0, inst.label() if inst else
                                    f"{current} ({'not found here; Scan cluster' if host else 'no charm here'})",
                                    str(current))
                self.env.setCurrentIndex(0)
            elif is_empty(current):
                self.env.insertItem(0, "Not set: choose a SimNIBS install" if self.env.count()
                                    else "No SimNIBS found: choose its bin folder with …", None)
                self.env.setCurrentIndex(0)

            self.version.clear()
            self.version.addItem(f"SimNIBS {inst.version}" if inst and inst.version else "–")
            if inst is None:
                note, col = "", c["ink3"]
            elif inst.supported is False:
                note, col = "not supported: PRESTUS needs SimNIBS 4 or later", c["red"]
            elif inst.supported is None:
                note, col = "version not found in this environment", c["ink3"]
            else:
                note, col = "follows the environment" + (f" · on {inst.host}" if inst.host else ""), c["ink3"]
            self.version_note.setText(note)
            self.version_note.setStyleSheet(f"color: {col};")
            self.linked.setText(self._linked_text(inst, c))

            m2m = self._m2m()
            st = sn.segmentation_state(m2m) if m2m else None
            settings = self.model.get(SETTINGS)
            rerun = self.model.get(RERUN) in (1, True)
            if st is None:
                self.state.setText(f"<span style='color:{c['ink3']}'>Set subject_id and path.seg to find it.</span>")
            elif not st.exists:
                self.state.setText(f"{st.folder.name}: none yet; charm makes it on the first run.")
            else:
                made = f"SimNIBS {st.version}" if st.version else "SimNIBS (version not logged)"
                used = ""
                defaults = inst.default_lines() if inst else None
                if st.settings and defaults is not None:
                    same = sn.normalised_settings(st.settings) == defaults
                    used = ", default charm settings" if same else ", custom charm settings"
                diffs = sn.differences(st, inst, None if is_empty(settings) else str(settings))
                text = f"{st.folder.name}: {made}{used}."
                if diffs:
                    hint = "it will be replaced" if rerun else "reused as is; Segment again… replaces it"
                    text += (f"<br><span style='color:{c['amber']}'>Differs from the selection: "
                             f"{'; '.join(diffs)} ({hint}).</span>")
                if rerun:
                    text += (f"<br><span style='color:{c['amber']}'>Re-run SimNIBS segmentation is on "
                             "(Setup › Input / output): every run of this config deletes "
                             f"{st.folder.name} and segments again.</span>")
                self.state.setText(text)
            self.report.setEnabled(bool(st and st.report))
            self.report.setVisible(bool(st and st.exists))

            show_settings = bool(self.page is not None and self.page.show_advanced) or not is_empty(settings)
            self.settings_label.setVisible(show_settings)
            self.settings_row.setVisible(show_settings)
            self.edit_settings.setVisible(not is_empty(settings) and Path(str(settings)).is_file())
        finally:
            self._updating = False
