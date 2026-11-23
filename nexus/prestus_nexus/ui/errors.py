"""Errors the user should see instead of a crash.

``install_excepthook`` makes an uncaught exception in a Qt slot a dialog (short message,
traceback under "Show Details…") plus a log entry, and the Nexus keeps running.
``guarded`` wraps one user action (save, open, new, reload) and turns the errors such an
action can meet -- file system, YAML, invalid values -- into the same dialog.

    with guarded(self, "Save study") as g:
        self.model.save(path)
    if g.failed:
        return False
"""
from __future__ import annotations

import contextlib
import sys
import threading
import time
import traceback
from pathlib import Path

LOG_NAME = "nexus_errors.log"
LOG_MAX_BYTES = 1_000_000

_showing = threading.local()            # a dialog is open on this thread: do not stack another


def log_path() -> Path | None:
    """The error log in the Nexus cache folder, if that folder exists."""
    try:
        from ..localcopy import cache_root
        root = cache_root()
    except Exception:  # best-effort: the error log must never cause an error
        return None
    return root / LOG_NAME if root.is_dir() else None


def log_exception(exc_type, exc, tb, context: str = "") -> str:
    """Write the traceback to stderr and the error log; return it as text."""
    text = "".join(traceback.format_exception(exc_type, exc, tb))
    head = time.strftime("%Y-%m-%d %H:%M:%S") + (f" {context}" if context else "")
    try:
        sys.stderr.write(f"[PRESTUS Nexus] {head}\n{text}")
    except Exception:  # best-effort: stderr may be closed or None (pythonw)
        pass
    path = log_path()
    if path is not None:
        try:
            if path.is_file() and path.stat().st_size > LOG_MAX_BYTES:
                path.replace(path.with_suffix(".log.1"))
            with open(path, "a", encoding="utf-8") as fh:
                fh.write(f"── {head}\n{text}\n")
        except OSError:  # best-effort: the log file is optional; stderr has the traceback
            pass
    return text


_logged_once: set[tuple[str, str, str]] = set()


def log_current(context: str) -> None:
    """Log the exception being handled (call inside an except block) without a dialog,
    for errors the Nexus works around but that should not vanish. The same error in the
    same context is logged once per session, so a repeated poll does not flood the log."""
    exc_type, exc, tb = sys.exc_info()
    if exc_type is None:
        return
    key = (context, exc_type.__name__, str(exc))
    if key in _logged_once:
        return
    _logged_once.add(key)
    log_exception(exc_type, exc, tb, context)


def short_message(exc: BaseException) -> str:
    """One readable line for exc (OSError: reason and file)."""
    if isinstance(exc, OSError) and (exc.strerror or exc.filename):
        where = f": {exc.filename}" if exc.filename else ""
        return f"{exc.strerror or type(exc).__name__}{where}"
    text = str(exc).strip()
    first = text.splitlines()[0] if text else ""
    return f"{type(exc).__name__}: {first}" if first else type(exc).__name__


def show_error(parent, title: str, message: str, details: str = "") -> None:
    """A warning dialog with expandable details. Never raises; at most one at a time."""
    if getattr(_showing, "on", False):
        return
    try:
        from PyQt6.QtWidgets import QApplication, QMessageBox
        if QApplication.instance() is None:
            return
        _showing.on = True
        box = QMessageBox(parent)
        box.setIcon(QMessageBox.Icon.Warning)
        box.setWindowTitle(title)
        box.setText(message)
        if details:
            box.setDetailedText(details)
        box.exec()
    except Exception:  # best-effort: an error dialog must never raise
        pass
    finally:
        _showing.on = False


def install_excepthook(parent_for=lambda: None) -> None:
    """Uncaught exceptions: log them and show a dialog; the Nexus keeps running.
    parent_for() gives the window to show the dialog over (or None)."""
    previous = sys.excepthook

    def hook(exc_type, exc, tb):
        if issubclass(exc_type, KeyboardInterrupt):
            previous(exc_type, exc, tb)
            return
        try:
            details = log_exception(exc_type, exc, tb, "uncaught")
            parent = None
            try:
                parent = parent_for()
            except Exception:  # best-effort: no parent window: the dialog opens on its own
                pass
            where = log_path()
            note = f"\n\nDetails were written to {where}." if where else ""
            show_error(parent, "PRESTUS Nexus: unexpected error",
                       "Something went wrong; the Nexus keeps running, but save your study "
                       f"before you continue.\n\n{short_message(exc)}{note}", details)
        except Exception:  # best-effort: the hook must never fail; Python's own hook reports it
            previous(exc_type, exc, tb)

    sys.excepthook = hook

    def thread_hook(args):             # worker threads: log only (no dialogs off the GUI thread)
        if args.exc_type is SystemExit:
            return
        log_exception(args.exc_type, args.exc_value, args.exc_traceback,
                      f"in thread {getattr(args.thread, 'name', '?')}")

    threading.excepthook = thread_hook


def action_errors() -> tuple[type[BaseException], ...]:
    """Errors a file action is expected to meet."""
    errors: list[type[BaseException]] = [OSError, ValueError, UnicodeError]
    try:
        import yaml
        errors.append(yaml.YAMLError)
    except ImportError:  # best-effort: without PyYAML there are no YAML errors to catch
        pass
    return tuple(errors)


class guarded(contextlib.ContextDecorator):
    """Context manager (or decorator) for one user action: OSError, YAML errors and
    ValueError become a dialog titled "<action_label> failed" and are logged; other
    exceptions pass on (to the excepthook). After the block, .failed and .error say
    what happened."""

    def __init__(self, parent, action_label: str, errors: tuple | None = None):
        self.parent = parent
        self.action_label = action_label
        self.errors = errors or action_errors()
        self.failed = False
        self.error: BaseException | None = None

    def __enter__(self) -> "guarded":
        self.failed, self.error = False, None
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc_type is None or not issubclass(exc_type, self.errors):
            return False
        self.failed, self.error = True, exc
        details = log_exception(exc_type, exc, tb, self.action_label)
        show_error(self.parent, f"{self.action_label} failed",
                   f"{self.action_label} failed.\n\n{short_message(exc)}", details)
        return True
