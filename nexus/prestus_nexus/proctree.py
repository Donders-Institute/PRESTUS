"""Stopping a started program together with everything it started.

``matlab`` is a launcher script on macOS and Linux and starts the MATLAB binary as a
child; killing only the launcher leaves MATLAB running (and holding a licence). On
POSIX, processes started with ``new_session`` lead their own process group, which
``kill_tree`` signals as a whole; children that left the group are found through the
parent-process table. On Windows, ``taskkill /T /F`` ends the tree.
"""
from __future__ import annotations

import os
import signal
import subprocess


def new_session(proc) -> None:
    """Have a QProcess start its program in a session (and process group) of its own
    (POSIX, Qt >= 6.6; elsewhere nothing changes and kill_tree walks the process table)."""
    if os.name == "nt":
        return
    try:
        from PyQt6.QtCore import QProcess
        params = QProcess.UnixProcessParameters()
        params.flags = QProcess.UnixProcessFlag.CreateNewSession
        proc.setUnixProcessParameters(params)
    except (AttributeError, ImportError, TypeError):
        pass


def _no_window() -> dict:
    return {"creationflags": getattr(subprocess, "CREATE_NO_WINDOW", 0)} if os.name == "nt" else {}


def descendants(pid: int) -> list[int]:
    """Process IDs below pid (POSIX, from ps); [] if they cannot be listed."""
    try:
        out = subprocess.run(["ps", "-A", "-o", "pid=", "-o", "ppid="], capture_output=True,
                             text=True, timeout=5, check=False).stdout
    except Exception:
        return []
    children: dict[int, list[int]] = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            children.setdefault(int(parts[1]), []).append(int(parts[0]))
    found, todo = [], [pid]
    while todo:
        for child in children.get(todo.pop(), []):
            if child not in found:
                found.append(child)
                todo.append(child)
    return found


def kill_tree(pid: int, sig: int | None = None) -> None:
    """Send sig (default: kill) to pid and all its descendants. Never raises."""
    if not pid or pid <= 0:
        return
    if os.name == "nt":
        try:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(pid)], capture_output=True,
                           timeout=15, check=False, **_no_window())
        except Exception:
            pass
        return
    sig = signal.SIGKILL if sig is None else sig
    tree = descendants(pid)
    try:
        if os.getpgid(pid) == pid:           # started with new_session: its own group
            os.killpg(pid, sig)
    except OSError:
        pass
    for p in [pid] + tree:
        try:
            os.kill(p, sig)
        except OSError:
            pass


def group_of(pid: int) -> int:
    """pid's process group when pid leads one (started with new_session), else 0."""
    if os.name == "nt" or not pid or pid <= 0:
        return 0
    try:
        return pid if os.getpgid(pid) == pid else 0
    except OSError:
        return 0


def group_alive(pgid: int) -> bool:
    """A member of process group pgid is still running (it outlives a killed leader)."""
    if os.name == "nt" or not pgid:
        return False
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True
    return True


def kill_group(pgid: int) -> None:
    if os.name == "nt" or not pgid:
        return
    try:
        os.killpg(pgid, signal.SIGKILL)
    except OSError:
        pass


def alive(pid: int) -> bool:
    if not pid or pid <= 0:
        return False
    if os.name == "nt":
        try:
            out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True,
                                 text=True, timeout=10, check=False, **_no_window()).stdout
        except Exception:
            return False
        return str(pid) in out
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True
    return True
