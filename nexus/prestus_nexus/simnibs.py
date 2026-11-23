"""SimNIBS installs and existing segmentations, read from the file system.

Nothing here runs SimNIBS. An install is a folder that holds ``charm``
(startup.simnibs_bin_path); its version comes from the simnibs package metadata
(``simnibs-X.Y.Z.dist-info``) in the environment that charm calls. An m2m folder
records what made it: ``charm_log.html`` (version on its first line),
``settings.ini`` (the charm settings used) and ``charm_report.html``.
MATLAB mirrors the version and settings checks (simnibs_install_info.m,
segmentation_check_existing.m).
"""

from __future__ import annotations

import glob
import os
import re
import shutil
import sys
from dataclasses import dataclass, replace
from pathlib import Path, PurePosixPath

MIN_MAJOR = 4                       # charm segmentation (labels in charm_seg_labels.m)
_LOG_VERSION = re.compile(r"INFO:\s*simnibs\s*version\s*([\d.]+)", re.IGNORECASE)
_WRAPPED = re.compile(r"""["']?([^"'\s]*site-packages[/\\]simnibs)[/\\]cli[/\\]charm\.py""")


def charm_names() -> tuple[str, ...]:
    return ("charm.cmd", "charm.exe", "charm.bat", "charm") if os.name == "nt" else ("charm",)


def charm_file(folder: str | Path) -> Path | None:
    for n in charm_names():
        p = Path(folder) / n
        if p.is_file():
            return p
    return None


def is_simnibs_bin(folder: str | Path) -> bool:
    return charm_file(folder) is not None


# ── one install ──────────────────────────────────────────────────────────

def package_dir(bin_path: str | Path) -> Path | None:
    """The simnibs package folder (…/site-packages/simnibs) that this bin folder's charm runs."""
    charm = charm_file(bin_path)
    if charm is None:
        return None
    # installer and postinstall links: a small wrapper naming …/site-packages/simnibs/cli/charm.py
    try:
        if charm.stat().st_size < 64_000:
            m = _WRAPPED.search(charm.read_text(errors="ignore"))
            if m and Path(m.group(1)).is_dir():
                return Path(m.group(1))
    except OSError:  # best-effort: not a readable wrapper: look for a conda environment next
        pass
    # conda environment: charm is a console script of the environment itself
    bins = [Path(bin_path), charm.resolve().parent]
    for b in bins:
        for root in (b.parent, b.parent / "simnibs_env"):
            for pat in ("lib/python3*/site-packages/simnibs", "Lib/site-packages/simnibs"):
                for hit in sorted(glob.glob(str(root / pat))):
                    if Path(hit).is_dir():
                        return Path(hit)
    return None


def package_version(pkg: Path | None) -> str | None:
    if pkg is None:
        return None
    for hit in glob.glob(str(pkg.parent / "simnibs-*.dist-info")):
        m = re.match(r"simnibs-([\d.]+[\w.+-]*)\.dist-info$", Path(hit).name)
        if m:
            return m.group(1)
    ver = pkg / "_version.py"
    if ver.is_file():
        m = re.search(r"""__version__\s*=\s*['"]([^'"]+)['"]""", ver.read_text(errors="ignore"))
        if m:
            return m.group(1)
    return None


def version_major(version: str | None) -> int | None:
    m = re.match(r"(\d+)", version or "")
    return int(m.group(1)) if m else None


def is_supported(version: str | None) -> bool | None:
    """True for SimNIBS >= 4, False below, None when the version is unknown."""
    major = version_major(version)
    return None if major is None else major >= MIN_MAJOR


@dataclass(frozen=True)
class Install:
    """A SimNIBS install: on this computer (host ''), or on the cluster as found by
    remote.simnibs_scan_command (paths there, with its charm.ini already read)."""
    bin: Path | PurePosixPath
    version: str | None
    package: Path | PurePosixPath | None
    host: str = ""                                  # '' = this computer
    module: str = ""                                # environment module that provides it
    pydeface: bool | None = None                    # pydeface in the bin folder (segmentation.deface)
    settings_lines: tuple[str, ...] | None = None   # remote: its charm.ini, normalised
    python: str = ""                                # interpreter charm runs (PlanTUS without env_path)
    ld_library_path: str = ""                       # <package>/external/lib/linux, if present

    @property
    def env_root(self) -> PurePosixPath | Path | None:
        """Environment folder that holds the simnibs package (…/lib/python3.x/site-packages/simnibs)."""
        if self.package is None:
            return None
        root = self.package.parent.parent               # lib/python3.x  (Windows: Lib)
        if root.name.startswith("python"):
            root = root.parent
        return root.parent

    def owns(self, path: str) -> bool | None:
        """Whether a folder belongs to this install: inside its environment, or the install folder
        that holds it (installer: <root>/simnibs_env). None when the install's folder is unknown."""
        root = self.env_root
        if root is None or not path:
            return None
        p, r = _norm(path, bool(self.host)), _norm(str(root), bool(self.host))
        return p == r or p.startswith(r + "/") or r.startswith(p + "/")

    @property
    def default_settings(self) -> Path | None:
        """charm.ini shipped with a local install (charm's settings when none are given)."""
        if self.host or self.package is None:
            return None
        p = Path(self.package) / "charm.ini"
        return p if p.is_file() else None

    def default_lines(self) -> list[str] | None:
        """The install's charm.ini, normalised (normalised_settings); None if unknown."""
        if self.settings_lines is not None:
            return list(self.settings_lines)
        p = self.default_settings
        return normalised_settings(p) if p else None

    @property
    def supported(self) -> bool | None:
        return is_supported(self.version)

    def label(self) -> str:
        ver = f"SimNIBS {self.version}" if self.version else "SimNIBS (version unknown)"
        where = f" · {self.module}" if self.module else ""
        on = f" · on {self.host}" if self.host else ""
        return f"{ver}{where} · {self.bin}{on}"

    def to_dict(self) -> dict:
        return {"bin": str(self.bin), "version": self.version, "package": str(self.package or ""),
                "host": self.host, "module": self.module, "pydeface": self.pydeface,
                "settings": list(self.settings_lines) if self.settings_lines is not None else None,
                "python": self.python, "ld_library_path": self.ld_library_path}

    @classmethod
    def from_dict(cls, d: dict) -> "Install":
        return cls(PurePosixPath(d["bin"]), d.get("version") or None,
                   PurePosixPath(d["package"]) if d.get("package") else None, d.get("host", ""),
                   d.get("module", ""), d.get("pydeface"),
                   tuple(d["settings"]) if d.get("settings") is not None else None,
                   d.get("python", "") or "", d.get("ld_library_path", "") or "")


def install_info(bin_path: str | Path) -> Install | None:
    """The install in this bin folder on this computer, or None if it holds no charm."""
    if not bin_path or not is_simnibs_bin(bin_path):
        return None
    pkg = package_dir(bin_path)
    pyd = any((Path(bin_path) / n).is_file() for n in ("pydeface", "pydeface.exe"))
    ld = pkg / "external" / "lib" / "linux" if pkg else None
    inst = Install(Path(bin_path), package_version(pkg), pkg, pydeface=pyd,
                   ld_library_path=str(ld) if ld and ld.is_dir() else "")
    return replace(inst, python=charm_python(bin_path, inst))


def charm_python(bin_path: str | Path, inst: Install | None = None) -> str:
    """The interpreter charm runs (as simnibs_install_info.m): the '#!' line of a conda console
    script, or the quoted python of an installer wrapper; else python in the environment."""
    charm = charm_file(bin_path)
    text = ""
    try:
        if charm is not None and charm.stat().st_size < 64_000:
            text = charm.read_text(errors="ignore")
    except OSError:  # best-effort: unreadable script: the environment's python is used below
        pass
    first = text.splitlines()[0] if text else ""
    if first.startswith("#!") and "/env " not in first:
        py = first[2:].strip().split(" ")[0]
        if "python" in py and Path(py).is_file():
            return py
    m = re.search(r'"([^"]*python[^"]*)"\s+-E', text)
    if m and Path(m.group(1)).is_file():
        return m.group(1)
    cands = [Path(bin_path) / "python"]
    root = inst.env_root if inst else None
    if root is not None:
        cands += [Path(root) / "bin" / "python", Path(root) / "python.exe"]
    return next((str(c) for c in cands if c.is_file()), "")


def _norm(path: str, remote: bool) -> str:
    if remote:
        return str(PurePosixPath(path)).rstrip("/")
    try:
        return Path(path).expanduser().resolve().as_posix().rstrip("/")
    except OSError:
        return Path(path).as_posix().rstrip("/")


# Config keys that follow startup.simnibs_bin_path while empty (load_parameters / PRESTUS):
# hpc.ld_library_path -> simnibs_ld_library_path.m; placement.plantus.env_path -> the
# install's Python (position_transducer_plantus.m)
LINKED_KEYS = ("hpc.ld_library_path", "placement.plantus.env_path")


def linked_mismatches(model_get, inst: Install | None) -> list[tuple[str, str]]:
    """Linked keys set explicitly to a folder of another SimNIBS install than the selected one."""
    out = []
    if inst is None:
        return out
    for key in LINKED_KEYS:
        v = model_get(key)
        if isinstance(v, str) and v and inst.owns(v) is False and \
                (key != "hpc.ld_library_path" or _looks_simnibs(v)):      # may be a system folder
            out.append((key, v))
    return out


def _looks_simnibs(path: str) -> bool:
    """A folder of some SimNIBS install (so a mismatch means another install, not a custom path)."""
    p = path.lower()
    return "simnibs" in p or "site-packages" in p


# ── installs on the cluster (found over SSH) ─────────────────────────────

_REMOTE: dict[str, Install] = {}                    # bin path on the host -> install
remote_host = lambda: None                          # noqa: E731  set by the main window: the SSH host
                                                    # runs are submitted to, or None (local runs)


_remote_generation = 0


def remote_generation() -> int:
    """Changes with every set_remote_installs (checks cache their answers on it)."""
    return _remote_generation


def set_remote_installs(host: str, installs: list[Install]) -> None:
    """Replace the installs known for one host (a new scan of it)."""
    global _remote_generation
    _remote_generation += 1
    for k in [k for k, v in _REMOTE.items() if v.host == host]:
        del _REMOTE[k]
    for i in installs:
        _REMOTE[str(i.bin)] = i


def remote_installs(host: str | None = None) -> list[Install]:
    found = [i for i in _REMOTE.values() if host is None or i.host == host]
    return sorted(found, key=_vkey, reverse=True)


def lookup(bin_path: str | Path | None) -> Install | None:
    """The install of startup.simnibs_bin_path: on this computer, else one found on the
    cluster that runs are submitted to."""
    if not bin_path:
        return None
    if Path(bin_path).is_dir():
        return install_info(bin_path)
    inst = _REMOTE.get(str(PurePosixPath(str(bin_path).rstrip("/"))))
    return inst if inst is not None and inst.host == (remote_host() or "") else None


# ── scanning ─────────────────────────────────────────────────────────────

def _conda_env_bins(home: Path) -> list[str]:
    pats = [str(home / d / "envs" / "*" / "bin")
            for d in (".conda", "miniconda3", "anaconda3", "miniforge3", "mambaforge")]
    out: list[str] = []
    for pat in pats:
        out += glob.glob(pat)
    txt = home / ".conda" / "environments.txt"           # conda's own list of environments
    if txt.is_file():
        for line in txt.read_text(errors="ignore").splitlines():
            if line.strip():
                out.append(str(Path(line.strip()) / ("Scripts" if os.name == "nt" else "bin")))
    return out


def candidate_bins(home: Path | None = None) -> list[str]:
    """Folders that may hold charm: PATH, conda environments, the installers' default
    locations and central installs (/opt on a cluster)."""
    home = home or Path.home()
    out: list[str] = []
    exe = shutil.which("charm")
    if exe:
        out.append(str(Path(exe).parent))
    out += _conda_env_bins(home)
    pats = [str(home / "SimNIBS*" / "bin"), str(home / "Applications" / "SimNIBS*" / "bin")]
    if sys.platform == "darwin":
        pats += ["/Applications/SimNIBS*/bin"]
    elif os.name == "nt":
        for env in ("LOCALAPPDATA", "PROGRAMFILES"):
            root = os.environ.get(env)
            if root:
                pats += [str(Path(root) / "SimNIBS*" / "bin")]
    else:
        pats += ["/opt/SimNIBS*/bin", "/opt/simnibs*/bin", "/opt/simnibs/*/bin", "/opt/simnibs/*/*/bin",
                 "/usr/local/SimNIBS*/bin"]
    for pat in pats:
        out += glob.glob(pat)
    return out


def _vkey(inst: Install) -> tuple:
    return tuple(int(x) for x in re.findall(r"\d+", inst.version or "")) or (-1,)


def find_installs(home: Path | None = None) -> list[Install]:
    """Every SimNIBS install found, newest version first (unknown versions last)."""
    seen: set[str] = set()
    found: list[Install] = []
    for b in candidate_bins(home):
        try:
            key = str(Path(b).resolve())
        except OSError:  # best-effort: a path that cannot be resolved is not an install
            continue
        if key in seen:
            continue
        seen.add(key)
        inst = install_info(b)
        if inst is not None:
            found.append(inst)
    return sorted(found, key=_vkey, reverse=True)


# ── an existing segmentation ─────────────────────────────────────────────

def m2m_version(m2m: Path) -> str | None:
    """SimNIBS version from the first line of charm_log.html (as simnibs_version.m reads it)."""
    log = Path(m2m) / "charm_log.html"
    if not log.is_file():
        return None
    try:
        with log.open(errors="ignore") as f:
            m = _LOG_VERSION.search(f.readline())
    except OSError:
        return None
    return m.group(1) if m else None


def normalised_settings(path: str | Path) -> list[str]:
    """Setting lines of a charm .ini without comments, blank lines and spacing (charm copies
    the file it used into the m2m folder verbatim, so the order is kept)."""
    out = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        line = re.sub(r"\s*=\s*", "=", " ".join(line.split("#", 1)[0].split()))
        if line:
            out.append(line)
    return out


def settings_equal(a: str | Path, b: str | Path) -> bool:
    """Same charm settings, ignoring comments, blank lines and spacing."""
    return normalised_settings(a) == normalised_settings(b)


@dataclass(frozen=True)
class SegmentationState:
    folder: Path
    exists: bool                    # final_tissues.nii.gz present
    version: str | None             # SimNIBS that made it (charm_log.html)
    settings: Path | None           # settings.ini charm used
    report: Path | None             # charm_report.html


def segmentation_state(m2m: str | Path) -> SegmentationState:
    m2m = Path(m2m)
    s, r = m2m / "settings.ini", m2m / "charm_report.html"
    return SegmentationState(m2m, (m2m / "final_tissues.nii.gz").is_file(), m2m_version(m2m),
                             s if s.is_file() else None, r if r.is_file() else None)


def requested_lines(settings_file: str | None, inst: Install | None) -> list[str] | None:
    """The settings a new charm run would use, normalised: segmentation.settings_file (when
    readable here), else the install's charm.ini; None when unknown."""
    if settings_file:
        p = Path(settings_file)
        return normalised_settings(p) if p.is_file() else None
    return inst.default_lines() if inst else None


def differences(state: SegmentationState, inst: Install | None, settings_file: str | None) -> list[str]:
    """How an existing segmentation differs from what the selected install and settings would
    make; empty when they match or cannot be compared."""
    out: list[str] = []
    if not state.exists:
        return out
    if inst is not None and inst.version and state.version and inst.version != state.version:
        out.append(f"made with SimNIBS {state.version}; the selected environment is {inst.version}")
    try:
        req = requested_lines(settings_file, inst)
        if req is not None and state.settings is not None and normalised_settings(state.settings) != req:
            out.append("made with other charm settings than "
                       + ("segmentation.settings_file" if settings_file else "the selected install's defaults"))
    except OSError:                    # the settings comparison is left out
        from .ui.errors import log_current
        log_current("SimNIBS: cannot read the charm settings to compare")
    return out
