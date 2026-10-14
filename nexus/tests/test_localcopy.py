"""Running the Nexus from a local copy of its code (network drives)."""

import os
import time

from prestus_nexus import localcopy


def _package(root):
    src = root / "prestus_nexus"
    (src / "ui").mkdir(parents=True)
    (src / "__init__.py").write_text("")
    (src / "app.py").write_text("def main(argv=None):\n    return 0\n")
    (src / "overlay.yaml").write_text("a: 1\n")
    (src / "ui" / "page.py").write_text("X = 1\n")
    (src / "notes.txt").write_text("not copied\n")
    return src


def test_mirror_copies_updates_and_removes(tmp_path):
    src = _package(tmp_path / "checkout")
    dest = localcopy.mirror(src, tmp_path / "cache") / "prestus_nexus"
    assert (dest / "app.py").is_file() and (dest / "ui" / "page.py").is_file()
    assert (dest / "overlay.yaml").is_file() and not (dest / "notes.txt").exists()
    (src / "ui" / "page.py").write_text("X = 22\n")
    later = time.time() + 5
    os.utime(src / "ui" / "page.py", (later, later))
    (src / "overlay.yaml").unlink()
    localcopy.mirror(src, tmp_path / "cache")
    assert (dest / "ui" / "page.py").read_text() == "X = 22\n"
    assert not (dest / "overlay.yaml").exists()


def test_local_main_opt_out_and_nesting(monkeypatch):
    monkeypatch.delenv(localcopy.SOURCE_ENV, raising=False)
    monkeypatch.setenv(localcopy.OPT_OUT_ENV, "0")
    assert localcopy.local_main() is None
    monkeypatch.delenv(localcopy.OPT_OUT_ENV)
    monkeypatch.setenv(localcopy.SOURCE_ENV, "/somewhere")    # already running from a copy
    assert localcopy.local_main() is None


def test_package_dir_follows_the_source(monkeypatch, tmp_path):
    from prestus_nexus import paths
    monkeypatch.setenv(localcopy.SOURCE_ENV, str(tmp_path / "nexus" / "prestus_nexus"))
    assert paths.package_dir() == tmp_path / "nexus" / "prestus_nexus"
