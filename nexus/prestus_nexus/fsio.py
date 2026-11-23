"""Atomic file writes.

The study folder usually lives on a network share that can fill up. Writing a file in place
truncates it before the new content lands, so a failed write would leave an empty or partial
config. Every write here goes to a temporary file in the same folder first and replaces the
target only once the new content is complete; on failure the old file is untouched.
"""
from __future__ import annotations

import os
import shutil
import tempfile
from pathlib import Path


def _umask() -> int:
    mask = os.umask(0)                   # reading the umask means setting it; restore at once
    os.umask(mask)
    return mask


def _replace_from_temp(path: Path, fill) -> None:
    path = Path(path)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as fh:
            fill(fh)
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError:  # best-effort: some network file systems refuse fsync; the replace still holds
                pass
        try:
            if path.exists():
                shutil.copymode(path, tmp)
            else:                        # mkstemp makes 0600; a new file gets the usual mode
                os.chmod(tmp, 0o666 & ~_umask())
        except OSError:  # best-effort: mode bits only; the content is written below
            pass
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:  # best-effort: temp file cleanup; the original error is raised
            pass
        raise


def write_bytes_atomic(path: str | Path, data: bytes) -> None:
    """Write data to path, replacing it only once the write has succeeded."""
    _replace_from_temp(Path(path), lambda fh: fh.write(data))


def write_text_atomic(path: str | Path, text: str, encoding: str = "utf-8") -> None:
    """Write text to path, replacing it only once the write has succeeded."""
    write_bytes_atomic(path, text.encode(encoding))


def copy_atomic(src: str | Path, dst: str | Path) -> None:
    """Copy src to dst (contents and modification time), replacing dst only on success."""
    src, dst = Path(src), Path(dst)

    def fill(fh):
        with open(src, "rb") as fin:
            shutil.copyfileobj(fin, fh)

    _replace_from_temp(dst, fill)
    try:
        shutil.copystat(src, dst)
    except OSError:  # best-effort: timestamps only; the content is copied
        pass
