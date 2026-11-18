"""Submit to a cluster over SSH, from a computer without sbatch/qsub.

The study folder must be visible on both sides, e.g. the project drive
mounted on a laptop (/Volumes/2425122.01) and the same drive on the cluster
(/project/2425122.01). Paths are translated by prefix (path_map): the study
folder, PRESTUS, and every config value under a mapped folder (path_overrides,
passed to load_parameters as overrides, so the config file stays as it is).
The Nexus then runs, on the login node,

    ssh <host> bash -lc '<setup> && cd <study dir> && matlab -batch "..."'

so PRESTUS submits exactly as it would on the cluster (prestus_pipeline_start
-> hpc_submit_job), and follows the jobs with squeue/sacct/qstat over SSH.

Commands run with BatchMode=yes, so ssh never prompts. They log in with a key, or
through a shared connection: *Connect* (connect_command) logs in once, asking for the
password in a Nexus dialog (SSH_ASKPASS, see askpass.py), and keeps an OpenSSH master
connection open for CONTROL_PERSIST; every later command reuses it (ControlPath).
"""

from __future__ import annotations

import os
import re
import shlex
import stat
import sys
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath

from . import launcher


CONTROL_PERSIST = "8h"        # how long the shared connection stays open once idle


def control_path() -> str:
    """Socket of the shared connection, one per user@host:port (%C is their hash)."""
    return str(Path.home() / ".ssh" / "prestus-nexus-%C")


def sharing_supported() -> bool:
    """OpenSSH on Windows has no connection sharing."""
    return os.name != "nt"


_WINDOWS_PATH = re.compile(r"^(?:[A-Za-z]:(?:[\\/]|$)|\\\\|//[^/])")


def _windows_style(text: str) -> bool:
    """A drive-letter (Z:\\, z:/) or UNC (\\\\server\\share, //server/share) path."""
    return bool(_WINDOWS_PATH.match(text))


def _under(path: str, prefix: str, fold: bool) -> bool:
    """path is prefix or lies below it (forward slashes); fold: ignore case (Windows)."""
    head = path[:len(prefix)]
    if fold:
        head, prefix = head.lower(), prefix.lower()
    return head == prefix and (len(path) == len(prefix) or path[len(prefix)] == "/")


def _local_prefixes(local: str) -> set[str]:
    """A mapped local folder as written and resolved, with forward slashes."""
    if _windows_style(local):
        out = {PureWindowsPath(local).as_posix().rstrip("/")}
        if os.name == "nt":              # a mapped drive may resolve to its UNC path
            out.add(Path(local).expanduser().resolve().as_posix().rstrip("/"))
        return out
    return {Path(local).expanduser().as_posix().rstrip("/"),
            Path(local).expanduser().resolve().as_posix().rstrip("/")}


class PathMapError(ValueError):
    """A local path is not under any mapped prefix."""


@dataclass
class RemoteHost:
    host: str                                   # user@login-node, or a ~/.ssh/config alias
    path_map: list[tuple[str, str]] = field(default_factory=list)   # (local prefix, remote prefix)
    setup: str = ""                             # shell commands run first, e.g. 'module load slurm matlab/R2023b'
    matlab: str = "matlab"                      # MATLAB command on the remote side
    prestus_root: str = ""                      # remote PRESTUS checkout; '' = mapped local root

    # ── paths ───────────────────────────────────────────────────────────
    def to_remote(self, path: str | Path) -> str:
        """Local path -> remote path through the longest matching prefix."""
        local = Path(path).resolve().as_posix()
        best: tuple[str, str] | None = None
        for lp, rp in self.path_map:
            lp_n = Path(lp).expanduser().resolve().as_posix().rstrip("/")
            if _under(local, lp_n, _windows_style(local) or _windows_style(lp_n)) and \
                    (best is None or len(lp_n) > len(best[0])):
                best = (lp_n, rp.rstrip("/"))
        if best is None:
            raise PathMapError(f"{local} is not under a mapped local folder "
                               f"({', '.join(lp for lp, _ in self.path_map) or 'none set'})")
        return best[1] + local[len(best[0]):]

    def map_value(self, text: str) -> str | None:
        """A config string under a mapped local folder, on the host; None if it is not.

        Prefix match on the text as written (it may be a pattern such as sub-%03d, so
        it is not resolved), also against the resolved local folder (/tmp vs /private/tmp).
        Windows paths (Z:\\data, \\\\server\\share) are compared with forward slashes and,
        as Windows does, without regard to case.
        """
        win = _windows_style(text)
        if win:
            text = PureWindowsPath(text).as_posix()
        best: tuple[str, str] | None = None
        for lp, rp in self.path_map:
            for prefix in _local_prefixes(lp):
                if prefix and _under(text, prefix, win or _windows_style(prefix)) and \
                        (best is None or len(prefix) > len(best[0])):
                    best = (prefix, rp.rstrip("/"))
        return None if best is None else best[1] + text[len(best[0]):]

    def to_local_or_none(self, remote_path: str) -> str | None:
        """Remote path -> local path through the longest matching remote prefix, or None."""
        best: tuple[str, str] | None = None
        for lp, rp in self.path_map:
            rp_n = rp.rstrip("/")
            if rp_n and (remote_path == rp_n or remote_path.startswith(rp_n + "/")) and \
                    (best is None or len(rp_n) > len(best[0])):
                best = (rp_n, Path(lp).expanduser().as_posix().rstrip("/"))
        return None if best is None else best[1] + remote_path[len(best[0]):]

    def remote_root(self, local_root: Path) -> str:
        return self.prestus_root.rstrip("/") or self.to_remote(local_root)

    # ── commands ────────────────────────────────────────────────────────
    def ssh(self, script: str) -> list[str]:
        """argv running script in a login shell on the host (so 'module' works)."""
        if self.setup.strip():
            script = f"{self.setup.strip()} && {script}"
        return ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", *self._shared(),
                self.host, "bash -lc " + shlex.quote(script)]

    def _shared(self) -> list[str]:
        """Reuse the connection opened by connect_command, if there is one."""
        return ["-o", f"ControlPath={control_path()}"] if sharing_supported() else []

    def wrap(self, argv: list[str]) -> list[str]:
        """A local argv (squeue, scancel, ...) run on the host instead."""
        return self.ssh(" ".join(shlex.quote(a) for a in argv))


def connect_command(remote: RemoteHost) -> list[str]:
    """ssh argv that logs in once (prompts through SSH_ASKPASS) and leaves a shared
    connection running in the background for CONTROL_PERSIST."""
    return ["ssh", "-o", "ControlMaster=yes", "-o", f"ControlPath={control_path()}",
            "-o", f"ControlPersist={CONTROL_PERSIST}", "-o", "ConnectTimeout=30",
            "-o", "ServerAliveInterval=60", remote.host, "true"]


def check_command(remote: RemoteHost) -> list[str]:
    """ssh argv that exits 0 when the shared connection is open."""
    return ["ssh", "-o", f"ControlPath={control_path()}", "-O", "check", remote.host]


def disconnect_command(remote: RemoteHost) -> list[str]:
    return ["ssh", "-o", f"ControlPath={control_path()}", "-O", "exit", remote.host]


def askpass_env(folder: Path | None = None) -> dict[str, str]:
    """Environment that makes ssh ask through the Nexus dialog (askpass.py).

    SSH_ASKPASS must be a single program, so a small wrapper script is written that
    starts this Python (or the frozen app) with --askpass. Needs OpenSSH >= 8.4 for
    SSH_ASKPASS_REQUIRE.
    """
    folder = folder or Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "prestus_nexus"
    folder.mkdir(parents=True, exist_ok=True)
    (Path.home() / ".ssh").mkdir(mode=0o700, exist_ok=True)
    if getattr(sys, "frozen", False):
        call = f'exec {shlex.quote(sys.executable)} --askpass "$@"'
    else:
        nexus_dir = Path(__file__).resolve().parents[1]
        call = (f'PYTHONPATH={shlex.quote(str(nexus_dir))}${{PYTHONPATH:+:$PYTHONPATH}} '
                f'exec {shlex.quote(sys.executable)} -m prestus_nexus --askpass "$@"')
    script = folder / "askpass.sh"
    script.write_text("#!/bin/sh\n" + call + "\n", encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return {"SSH_ASKPASS": str(script), "SSH_ASKPASS_REQUIRE": "force"}


def host_of(argv: list[str]) -> str:
    """The host in an ssh argv (the first word that is not an option or its value)."""
    i = 1
    while i < len(argv):
        if argv[i] in ("-o", "-O", "-p", "-i", "-l"):
            i += 2
        elif argv[i].startswith("-"):
            i += 1
        else:
            return argv[i]
    return ""


# Paths that can only be on this (macOS) computer; /home/... may well be a cluster path
LOCAL_ONLY = ("/Volumes/", "/Users/")


def path_overrides(model, remote: RemoteHost) -> tuple[list[tuple[str, str]], list[str]]:
    """Config values to translate for the host, as (dotted key, remote value) overrides
    for load_parameters, and keys whose value looks local but is under no mapped folder.

    Walks the merged config (toolbox + project defaults + study, as model.effective());
    transducers are keyed 'transducer.x' (one) or 'transducer(i).x' (several), as
    apply_config_overrides expects.
    """
    from .model import iter_leaves
    eff = model.effective()
    txs = eff.get("transducer")
    leaves = list(iter_leaves({k: v for k, v in eff.items() if k != "transducer"}))
    if isinstance(txs, list):
        for i, tx in enumerate(txs, 1):
            leaves += [(f"transducer({i}).{k}", v) for k, v in iter_leaves(tx)]
    elif isinstance(txs, dict):
        leaves += list(iter_leaves(txs, "transducer"))
    mapped, unmapped = [], []
    for key, value in leaves:
        if not isinstance(value, str) or not value:
            continue
        new = remote.map_value(value)
        if new is not None:
            mapped.append((key, new))
        elif value.startswith(LOCAL_ONLY) and not any(r and value.startswith(r) for _, r in remote.path_map):
            unmapped.append(key)
    return mapped, unmapped


def _overrides_cell(overrides: list[tuple[str, object]]) -> str:
    return "{" + "; ".join(f"{launcher._q(k)}, {launcher.matlab_value(v)}" for k, v in overrides) + "}"


def submit_command(remote: RemoteHost, study: Path, prestus_root: Path,
                   overrides: list[tuple[str, object]] | None = None, extra: str = "",
                   options=None, sequential_overrides: dict | None = None) -> list[str]:
    """ssh argv that loads the study config on the host and submits it with PRESTUS.

    Keeps the config's platform and hpc.* settings ('auto' detects the host's scheduler);
    overrides (from path_overrides) replace config values after merging, before
    load_parameters checks the paths. extra: MATLAB statements run on ``parameters``
    before submitting (e.g. launcher.PLACEMENT_ONLY). options (launcher.PipelineOptions)
    are passed to prestus_pipeline_start; each sequential config is loaded from its folder on
    the host with its own path overrides (sequential_overrides: resolved path -> overrides).
    """
    study = Path(study).resolve()
    rdir = PurePosixPath(remote.to_remote(study.parent))
    rroot = remote.remote_root(prestus_root)
    cell = f", {_overrides_cell(overrides)}" if overrides else ""
    seq_over = sequential_overrides or {}

    def load_remote(p: Path) -> str:
        ov = seq_over.get(str(p))
        c = f", {_overrides_cell(ov)}" if ov else ""
        return (f"load_parameters({launcher._q(p.name)}, "
                f"{launcher._q(str(PurePosixPath(remote.to_remote(p.parent))))}{c})")

    opts = launcher.options_lines(options, local=False, load_expr=load_remote)
    start = "prestus_pipeline_start(parameters, options)" if opts else "prestus_pipeline_start(parameters)"
    body = (f"addpath(genpath(fullfile({launcher._q(rroot)}, 'functions'))); "
            f"parameters = load_parameters({launcher._q(study.name)}, {launcher._q(str(rdir))}{cell}); "
            f"{extra + ' ' if extra else ''}{' '.join(opts) + ' ' if opts else ''}{start}")
    matlab = f"{remote.matlab} -batch {shlex.quote(launcher._wrapped(body))}"
    return remote.ssh(f"cd {shlex.quote(str(rdir))} && {matlab}")


def probe_command(remote: RemoteHost, prestus_root: Path) -> list[str]:
    """ssh argv that reports the scheduler, MATLAB and PRESTUS root on the host."""
    rroot = remote.remote_root(prestus_root)
    script = ("echo host: $(hostname); "
              "echo sbatch: $(command -v sbatch || echo -); "
              "echo qsub: $(command -v qsub || echo -); "
              f"echo matlab: $(command -v {shlex.quote(remote.matlab)} || echo -); "
              f"test -d {shlex.quote(rroot + '/functions')} "
              f"&& echo prestus: {shlex.quote(rroot)} || echo prestus: - {shlex.quote(rroot)}")
    return remote.ssh(script)


# Lists the SimNIBS installs on the host, as simnibs.find_installs does locally: charm on
# PATH, conda environments, ~/SimNIBS*, /opt/simnibs*, and the bin folders of the
# 'simnibs' environment modules (module show). Per install one 'simnibs' line
# (bin, version, package, pydeface 0/1, module, python, library folder) and its charm.ini as 'ini' lines,
# normalised like simnibs.normalised_settings. Runs nothing from SimNIBS.
SIMNIBS_SCAN = r'''
cands=""; mods=""
add() { cands="$cands
$1"; }
c=$(command -v charm 2>/dev/null) && add "$(dirname "$c")"
for d in "$HOME"/.conda/envs/*/bin "$HOME"/miniconda3/envs/*/bin "$HOME"/anaconda3/envs/*/bin \
         "$HOME"/miniforge3/envs/*/bin "$HOME"/mambaforge/envs/*/bin "$HOME"/SimNIBS*/bin \
         /opt/SimNIBS*/bin /opt/simnibs*/bin /opt/simnibs/*/bin /opt/simnibs/*/*/bin /usr/local/SimNIBS*/bin; do
  add "$d"
done
if [ -f "$HOME/.conda/environments.txt" ]; then
  while IFS= read -r e; do [ -n "$e" ] && add "$e/bin"; done < "$HOME/.conda/environments.txt"
fi
if type module >/dev/null 2>&1; then
  for m in $(module -t avail simnibs 2>&1 | grep -i '^simnibs' | sed 's/(default)//'); do
    for p in $(module show "$m" 2>&1 | grep -E '(^|[^A-Z_])PATH([^A-Z_]|$)' | grep -oE '/[^ "),:;]+'); do
      add "$p"; mods="$mods
$p	$m"
    done
  done
fi
seen="|"
echo "$cands" | while IFS= read -r b; do
  [ -n "$b" ] && [ -f "$b/charm" ] || continue
  r=$(cd "$b" 2>/dev/null && pwd -P) || continue
  case "$seen" in *"|$r|"*) continue;; esac
  seen="$seen$r|"
  pkg=$(grep -o "[^\"' ]*site-packages/simnibs/cli/charm.py" "$b/charm" 2>/dev/null | head -n 1)
  pkg=${pkg%/cli/charm.py}
  if [ -z "$pkg" ] || [ ! -d "$pkg" ]; then
    pkg=""
    for p in "$b"/../lib/python3*/site-packages/simnibs "$b"/../simnibs_env/lib/python3*/site-packages/simnibs \
             "$r"/../lib/python3*/site-packages/simnibs; do
      if [ -d "$p" ]; then pkg=$(cd "$p" && pwd -P); break; fi
    done
  fi
  ver=""
  if [ -n "$pkg" ]; then
    di=$(ls -d "$pkg"/../simnibs-*.dist-info 2>/dev/null | head -n 1)
    [ -n "$di" ] && ver=$(basename "$di" | sed -e 's/^simnibs-//' -e 's/\.dist-info$//')
    if [ -z "$ver" ] && [ -f "$pkg/_version.py" ]; then
      ver=$(sed -n "s/^__version__ *= *['\"]\([^'\"]*\)['\"].*/\1/p" "$pkg/_version.py" | head -n 1)
    fi
  fi
  pyd=0; [ -f "$b/pydeface" ] && pyd=1
  mod=$(echo "$mods" | awk -F'\t' -v p="$b" '$1==p {print $2; exit}')
  py=""; l1=$(head -n 1 "$b/charm" 2>/dev/null)
  case "$l1" in *'/env '*) ;; '#!'*python*) py=${l1#??}; py=${py# }; py=${py%% *};; esac
  [ -n "$py" ] && [ ! -f "$py" ] && py=""
  if [ -z "$py" ]; then
    py=$(grep -o '"[^"]*python[^"]*" -E' "$b/charm" 2>/dev/null | head -n 1 | sed -e 's/^"//' -e 's/" -E$//')
    [ -n "$py" ] && [ ! -f "$py" ] && py=""
  fi
  [ -z "$py" ] && [ -f "$b/python" ] && py="$b/python"
  ld=""; [ -n "$pkg" ] && [ -d "$pkg/external/lib/linux" ] && ld="$pkg/external/lib/linux"
  printf 'simnibs\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$b" "$ver" "$pkg" "$pyd" "$mod" "$py" "$ld"
  if [ -n "$pkg" ] && [ -f "$pkg/charm.ini" ]; then
    sed 's/#.*//' "$pkg/charm.ini" | tr -s ' \t' '  ' | sed -e 's/^ //' -e 's/ $//' -e 's/ *= */=/g' \
      | grep -v '^$' | while IFS= read -r l; do printf 'ini\t%s\t%s\n' "$b" "$l"; done
  fi
done
echo simnibs-scan-done
'''


def simnibs_scan_command(remote: RemoteHost) -> list[str]:
    """ssh argv that lists the SimNIBS installs on the host (parse_simnibs_scan)."""
    return remote.ssh(SIMNIBS_SCAN)


def parse_simnibs_scan(text: str, host: str):
    """Output of simnibs_scan_command -> [simnibs.Install] on host, or None when the scan
    did not finish (no connection, failed setup commands)."""
    from pathlib import PurePosixPath
    from .simnibs import Install
    if "simnibs-scan-done" not in text:
        return None
    rows: dict[str, dict] = {}
    ini: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if parts[0] == "simnibs" and len(parts) >= 5:
            rows[parts[1]] = {"version": parts[2] or None, "package": parts[3],
                              "pydeface": parts[4] == "1", "module": parts[5] if len(parts) > 5 else "",
                              "python": parts[6] if len(parts) > 6 else "",
                              "ld": parts[7] if len(parts) > 7 else ""}
        elif parts[0] == "ini" and len(parts) >= 3:
            ini.setdefault(parts[1], []).append("\t".join(parts[2:]))
    return [Install(PurePosixPath(b), r["version"], PurePosixPath(r["package"]) if r["package"] else None,
                    host=host, module=r["module"], pydeface=r["pydeface"],
                    settings_lines=tuple(ini[b]) if b in ini else None,
                    python=r["python"], ld_library_path=r["ld"])
            for b, r in rows.items()]


def cat_command(remote: RemoteHost, path: str) -> list[str]:
    """ssh argv that prints a file on the host (e.g. a cluster install's charm.ini)."""
    return remote.ssh(f"cat {shlex.quote(path)}")


def host_from_session(session) -> RemoteHost | None:
    """The SSH host runs are submitted to (Run tab, saved in the session), or None."""
    if not session.get("remote_enabled", False) or not str(session.get("remote_host", "") or "").strip():
        return None
    local, rem = session.get("remote_local_prefix", ""), session.get("remote_remote_prefix", "")
    return RemoteHost(host=str(session.get("remote_host")).strip(),
                      path_map=[(local, rem)] if local and rem else [],
                      setup=str(session.get("remote_setup", "") or ""),
                      matlab=str(session.get("remote_matlab", "") or "matlab"),
                      prestus_root=str(session.get("remote_prestus_root", "") or ""))


def parse_probe(text: str) -> dict[str, str]:
    """'key: value' lines from probe_command -> {key: value}; '-' means not found."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep and key.strip() in {"host", "sbatch", "qsub", "matlab", "prestus"}:
            out[key.strip()] = value.strip()
    return out


def scheduler_from_log(text: str) -> str | None:
    """'slurm'/'qsub' as prestus_pipeline_start reported it ('➤ auto-detected: SLURM')."""
    for line in text.splitlines():
        if "auto-detected:" in line or "deploying:" in line:
            word = line.rsplit(":", 1)[-1].strip().lower()
            if word in ("slurm", "qsub"):
                return word
    return None

