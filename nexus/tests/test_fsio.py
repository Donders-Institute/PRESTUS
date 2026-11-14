"""Atomic writes (fsio) and saving a study: a failed write leaves the old file intact."""

import os
import stat

import pytest

from prestus_nexus import fsio, yamlio
from prestus_nexus.model import ConfigModel


def test_round_trip(tmp_path):
    f = tmp_path / "a.txt"
    fsio.write_text_atomic(f, "één\n")
    assert f.read_text(encoding="utf-8") == "één\n"
    fsio.write_bytes_atomic(f, b"\x00\x01")
    assert f.read_bytes() == b"\x00\x01"
    g = tmp_path / "b.txt"
    fsio.copy_atomic(f, g)
    assert g.read_bytes() == b"\x00\x01"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.txt", "b.txt"]   # no temp files left


@pytest.mark.skipif(os.name == "nt", reason="POSIX permission bits")
def test_new_file_gets_the_usual_mode_and_existing_keeps_its_own(tmp_path):
    old = os.umask(0o022)
    try:
        f = tmp_path / "new.yaml"
        fsio.write_text_atomic(f, "x")
        assert stat.S_IMODE(f.stat().st_mode) == 0o644          # not mkstemp's 0600
        f.chmod(0o640)
        fsio.write_text_atomic(f, "y")
        assert stat.S_IMODE(f.stat().st_mode) == 0o640
    finally:
        os.umask(old)


def test_failed_write_keeps_the_old_file(tmp_path, monkeypatch):
    f = tmp_path / "config.yaml"
    f.write_text("old: 1\n", encoding="utf-8")

    def boom(fh):
        fh.write(b"partial")
        raise OSError(28, "No space left on device")
    with pytest.raises(OSError):
        fsio._replace_from_temp(f, boom)
    assert f.read_text(encoding="utf-8") == "old: 1\n"
    assert [p.name for p in tmp_path.iterdir()] == ["config.yaml"]


@pytest.mark.skipif(os.name == "nt" or (hasattr(os, "geteuid") and os.geteuid() == 0),
                    reason="needs POSIX permissions and a non-root user")
def test_save_to_read_only_folder_raises_and_keeps_old_file(schema, tmp_path):
    folder = tmp_path / "study"
    folder.mkdir()
    target = folder / "config.yaml"
    m = ConfigModel(schema)
    m.set("subject_id", 7)
    m.save(target)
    before = target.read_bytes()
    m.set("subject_id", 8)
    folder.chmod(0o555)
    try:
        with pytest.raises(OSError):
            m.save(target)
        assert target.read_bytes() == before                     # not truncated
    finally:
        folder.chmod(0o755)
    assert yamlio.load_file(target)["subject_id"] == 7
