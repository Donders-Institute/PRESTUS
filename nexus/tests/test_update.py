"""Updating PRESTUS to a release (update.py): every starting state of a checkout, the
headless command line, and the Help ▸ Check for updates… dialog. The 'remote' is a local
git repository, so nothing touches the network."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from prestus_nexus import update

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")

GIT_ENV = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.org",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.org",
               GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_SYSTEM=os.devnull,
               # submodules of local test repositories (git blocks file:// ones by default)
               GIT_CONFIG_COUNT="1", GIT_CONFIG_KEY_0="protocol.file.allow",
               GIT_CONFIG_VALUE_0="always")


def git(cwd: Path, *args: str) -> str:
    out = subprocess.run(["git", "-C", str(cwd), *args], capture_output=True, text=True,
                         env=GIT_ENV, check=True)
    return out.stdout.strip()


def commit(repo: Path, version: str, env_text: str | None = None) -> str:
    (repo / "VERSION").write_text(version + "\n", encoding="utf-8")
    if env_text is not None:
        (repo / "nexus" / "environment.yml").write_text(env_text, encoding="utf-8")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", version)
    return git(repo, "rev-parse", "HEAD")


@pytest.fixture
def upstream(tmp_path) -> Path:
    """The release repository: v0.1.0 (lightweight tag), v0.2.0 (annotated, changes the
    Nexus environment), a pre-release tag that is not a release, and one later commit."""
    repo = tmp_path / "upstream"
    for sub in ("config", "functions", "nexus"):
        (repo / sub).mkdir(parents=True)
    (repo / "config" / "config_default.yaml").write_text("a: 1\n", encoding="utf-8")
    (repo / "functions" / "f.m").write_text("% f\n", encoding="utf-8")
    git(repo, "init", "-q", "-b", "main")
    commit(repo, "v0.1.0", env_text="deps: [a]\n")
    git(repo, "tag", "v0.1.0")
    commit(repo, "v0.2.0", env_text="deps: [a, b]\n")
    git(repo, "tag", "-a", "v0.2.0", "-m", "release")
    commit(repo, "v0.3.0-dev")
    git(repo, "tag", "v0.3.0-rc1")
    return repo


def clone_at(upstream: Path, dest: Path, ref: str, branch: str | None = None) -> Path:
    """A user's clone standing on `ref`: detached, or on a new branch `branch`."""
    subprocess.run(["git", "clone", "-q", str(upstream), str(dest)], env=GIT_ENV, check=True)
    git(dest, "checkout", "-q", *(["-b", branch] if branch else ["--detach"]), ref)
    # the clone must not know the tags beyond its commit yet: an update has to fetch them
    for tag in git(dest, "tag", "--no-merged", "HEAD").split():
        git(dest, "tag", "-d", tag)
    return dest


@pytest.fixture(autouse=True)
def _git_identity(monkeypatch):
    for k in ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL", "GIT_COMMITTER_NAME", "GIT_COMMITTER_EMAIL",
              "GIT_CONFIG_GLOBAL", "GIT_CONFIG_SYSTEM", "GIT_CONFIG_COUNT", "GIT_CONFIG_KEY_0",
              "GIT_CONFIG_VALUE_0"):
        monkeypatch.setenv(k, GIT_ENV[k])


# ── release tags ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("tag, key", [("v0.7.0", (0, 7, 0)), ("v10.2.33", (10, 2, 33)),
                                      ("v0.7.0-dev", None), ("v0.3.0-rc1", None),
                                      ("0.7.0", None), ("v0.7", None), ("nightly", None)])
def test_release_key(tag, key):
    assert update.release_key(tag) == key


def test_remote_releases_keep_only_releases_and_name_commits(upstream):
    tags = update.remote_releases(str(upstream))
    assert sorted(tags) == ["v0.1.0", "v0.2.0"]                    # no pre-release tag
    for tag, sha in tags.items():                                  # annotated tag: its commit
        assert sha == git(upstream, "rev-parse", f"{tag}^{{commit}}")
    assert update.latest_release(tags) == "v0.2.0"
    assert update.latest_release({"v0.9.0": "", "v0.10.0": ""}) == "v0.10.0"   # not by text


def test_unreachable_remote_is_reported(tmp_path, upstream):
    st = update.check(clone_at(upstream, tmp_path / "c", "v0.1.0"), str(tmp_path / "nowhere"))
    assert st.error and st.relation == "unknown" and not st.update_available


# ── check: every starting state ──────────────────────────────────────────────

def _diverged(upstream, dest):
    c = clone_at(upstream, dest, "v0.1.0", branch="mine")
    (c / "functions" / "mine.m").write_text("% mine\n", encoding="utf-8")
    git(c, "add", "-A")
    git(c, "commit", "-q", "-m", "local work")
    return c


STATES = {   # name: (how the clone stands, relation to v0.2.0, does apply() move it)
    "on_old_release": (lambda u, d: clone_at(u, d, "v0.1.0"), "behind", True),
    "branch_behind": (lambda u, d: clone_at(u, d, "v0.1.0", branch="main2"), "behind", True),
    "on_latest_release": (lambda u, d: clone_at(u, d, "v0.2.0"), "current", False),
    "branch_ahead": (lambda u, d: clone_at(u, d, "main", branch="dev"), "ahead", False),
    "branch_diverged": (_diverged, "diverged", False),
}


@pytest.mark.parametrize("name", STATES)
def test_check_and_apply_per_state(name, upstream, tmp_path):
    make, relation, moves = STATES[name]
    c = make(upstream, tmp_path / "clone")
    before = git(c, "rev-parse", "HEAD")
    branch = git(c, "branch", "--show-current")

    st = update.check(c, str(upstream))
    assert (st.relation, st.latest, st.clone, st.error) == (relation, "v0.2.0", True, "")
    assert st.update_available == (relation == "behind")
    assert st.branch == (branch or None)
    assert git(c, "rev-parse", "HEAD") == before                   # a check changes nothing

    res = update.apply(c, remote=str(upstream))
    assert res.changed == moves
    assert res.ok == (relation != "diverged")
    after = git(c, "rev-parse", "HEAD")
    if moves:
        assert after == git(upstream, "rev-parse", "v0.2.0^{commit}")
        assert (c / "VERSION").read_text().strip() == "v0.2.0"
        assert res.env_changed                                     # environment.yml differs
        assert any("install.sh" in m for m in res.messages)
        assert update.check(c, str(upstream)).relation == "current"
    else:
        assert after == before
    assert git(c, "branch", "--show-current") == branch            # never leaves the branch
    if name == "branch_diverged":
        assert (c / "functions" / "mine.m").is_file()              # local commits kept
        assert "checkout v0.2.0" in res.messages[-1] and "Ask whoever maintains" in res.messages[-1]


def test_local_changes_stop_the_update(upstream, tmp_path):
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    before = git(c, "rev-parse", "HEAD")
    (c / "functions" / "f.m").write_text("% my edit\n", encoding="utf-8")
    st = update.check(c, str(upstream))
    assert st.changed == ["functions/f.m"] and st.update_available and not st.can_update
    assert "block the update" in st.headline()
    res = update.apply(c, remote=str(upstream))
    assert not res.ok and not res.changed and "functions/f.m" in res.messages[0]
    assert git(c, "rev-parse", "HEAD") == before
    assert (c / "functions" / "f.m").read_text() == "% my edit\n"  # the edit survives


def test_untracked_files_do_not_stop_the_update(upstream, tmp_path):
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    (c / "my_study.yaml").write_text("x: 1\n", encoding="utf-8")
    assert update.check(c, str(upstream)).changed == []
    assert update.apply(c, remote=str(upstream)).changed
    assert (c / "my_study.yaml").is_file()


def test_explicit_target(upstream, tmp_path):
    c = clone_at(upstream, tmp_path / "clone", "v0.2.0")
    res = update.apply(c, "v0.1.0", str(upstream))                 # detached: going back is allowed
    assert res.ok and res.changed and (c / "VERSION").read_text().strip() == "v0.1.0"
    bad = update.apply(c, "v0.3.0-rc1", str(upstream))             # not a release
    assert not bad.ok and "not a PRESTUS release" in bad.messages[0]
    b = clone_at(upstream, tmp_path / "b", "v0.2.0", branch="main2")
    res = update.apply(b, "v0.1.0", str(upstream))                 # a branch is never moved back
    assert res.ok and not res.changed and (b / "VERSION").read_text().strip() == "v0.2.0"


def test_environment_unchanged_is_not_flagged(upstream, tmp_path):
    commit(upstream, "v0.2.1")
    git(upstream, "tag", "v0.2.1")
    c = clone_at(upstream, tmp_path / "clone", "v0.2.0")
    res = update.apply(c, remote=str(upstream))
    assert res.changed and res.target == "v0.2.1" and not res.env_changed
    assert not any("install.sh" in m for m in res.messages)


def test_installed_submodules_follow_the_release_and_others_stay_out(tmp_path):
    lib = tmp_path / "lib"
    lib.mkdir()
    git(lib, "init", "-q", "-b", "main")
    (lib / "lib.m").write_text("% 1\n", encoding="utf-8")
    git(lib, "add", "-A")
    git(lib, "commit", "-q", "-m", "1")
    up = tmp_path / "upstream"
    for sub in ("config", "functions", "nexus"):
        (up / sub).mkdir(parents=True)
    (up / "config" / "config_default.yaml").write_text("a: 1\n", encoding="utf-8")
    (up / "functions" / "f.m").write_text("% f\n", encoding="utf-8")
    git(up, "init", "-q", "-b", "main")
    for name in ("used", "optional"):
        git(up, "submodule", "add", "-q", str(lib), f"external/{name}")
    commit(up, "v0.1.0")
    git(up, "tag", "v0.1.0")
    (lib / "lib.m").write_text("% 2\n", encoding="utf-8")
    git(lib, "commit", "-q", "-am", "2")
    for name in ("used", "optional"):
        git(up / "external" / name, "pull", "-q", "origin", "main")
    commit(up, "v0.2.0")
    git(up, "tag", "v0.2.0")

    c = clone_at(up, tmp_path / "clone", "v0.1.0")
    git(c, "submodule", "update", "-q", "--init", "--", "external/used")
    assert update._installed_submodules(c) == ["external/used"]
    res = update.apply(c, remote=str(up))
    assert res.ok and res.changed
    assert (c / "external" / "used" / "lib.m").read_text() == "% 2\n"   # moved to the new pin
    assert not (c / "external" / "optional" / "lib.m").exists()          # not installed: left out


def test_downloaded_folder_is_not_a_clone(upstream, tmp_path):
    folder = tmp_path / "PRESTUS-0.1.0"
    shutil.copytree(upstream, folder, ignore=shutil.ignore_patterns(".git"))
    st = update.check(folder, str(upstream))
    assert not st.clone and st.latest == "v0.2.0" and "not installed with git" in st.error
    res = update.apply(folder, remote=str(upstream))
    assert not res.ok and "not installed with git" in res.messages[0]
    # also when the folder sits inside some other repository
    outer = tmp_path / "project"
    shutil.copytree(folder, outer / "tools" / "PRESTUS")
    git(outer, "init", "-q")
    assert not update.is_clone(outer / "tools" / "PRESTUS")


# ── headless command line ────────────────────────────────────────────────────

def test_cli_check_then_update(upstream, tmp_path, capsys):
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    common = ["--prestus-root", str(c), "--remote", str(upstream)]
    assert update.main(["--check", *common]) == update.EXIT_AVAILABLE
    out = capsys.readouterr().out
    assert "Latest release: v0.2.0" in out and "prestus-nexus --update" in out
    assert (c / "VERSION").read_text().strip() == "v0.1.0"
    assert update.main(common) == update.EXIT_OK
    assert "PRESTUS is now on v0.2.0" in capsys.readouterr().out
    assert update.main(["--check", *common]) == update.EXIT_OK
    assert update.main(common) == update.EXIT_OK                   # again: nothing to do


def test_cli_failures_exit_nonzero(upstream, tmp_path, capsys):
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    (c / "functions" / "f.m").write_text("% edit\n", encoding="utf-8")
    assert update.main(["--prestus-root", str(c), "--remote", str(upstream)]) == update.EXIT_FAILED
    assert update.main(["--prestus-root", str(tmp_path)]) == update.EXIT_FAILED   # no PRESTUS
    assert "not a PRESTUS folder" in capsys.readouterr().err


def _no_qt_env() -> dict:
    """A child whose imports of PyQt6 fail: proves the path is headless. Its console is
    cp1252 (a Windows default) on every platform, so the report must survive "▸" and "…"."""
    return dict(GIT_ENV, QT_QPA_PLATFORM="", DISPLAY="", PYTHONIOENCODING="cp1252")


BLOCK_QT = ("import sys\n"
            "class _B:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] == 'PyQt6': raise ImportError('Qt must not load')\n"
            "sys.meta_path.insert(0, _B())\n")


def test_entry_point_updates_without_qt(upstream, tmp_path):
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    code = BLOCK_QT + ("from prestus_nexus.app import main\n"
                       f"sys.exit(main(['--update', '--prestus-root', {str(c)!r}, "
                       f"'--remote', {str(upstream)!r}]))\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, encoding="cp1252",
                         env=_no_qt_env(), cwd=Path(update.__file__).parents[1])
    assert out.returncode == 0, out.stdout + out.stderr
    assert (c / "VERSION").read_text().strip() == "v0.2.0"
    assert "Update project defaults" in out.stdout      # the report got to the end


def test_runs_as_a_plain_script(upstream, tmp_path):
    """python3 nexus/prestus_nexus/update.py: no installed package, no environment."""
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    out = subprocess.run([sys.executable, "-S", update.__file__, "--check", "--prestus-root",
                          str(c), "--remote", str(upstream)],
                         capture_output=True, text=True, env=GIT_ENV, cwd=tmp_path)
    assert out.returncode == update.EXIT_AVAILABLE, out.stdout + out.stderr


def test_default_root_follows_the_local_copy_source(monkeypatch, tmp_path):
    monkeypatch.delenv("PRESTUS_ROOT", raising=False)
    monkeypatch.setenv("PRESTUS_NEXUS_SOURCE", str(tmp_path / "P" / "nexus" / "prestus_nexus"))
    assert update.default_root() == tmp_path / "P"


# ── dialog ───────────────────────────────────────────────────────────────────

@pytest.fixture
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _finish(dialog):
    dialog._worker.wait()
    from PyQt6.QtWidgets import QApplication
    QApplication.processEvents()


@pytest.mark.parametrize("name", STATES)
def test_dialog_offers_the_update_only_when_it_can_run(name, upstream, tmp_path, qapp):
    from prestus_nexus.ui.update_dialog import UpdateDialog
    make, relation, _moves = STATES[name]
    d = UpdateDialog(make(upstream, tmp_path / "clone"), remote=str(upstream))
    d.start_check()
    _finish(d)
    assert d.status.relation == relation
    assert d.update_btn.isEnabled() == (relation == "behind")
    assert "v0.2.0" in d.details.text()


def test_dialog_updates_after_confirmation(upstream, tmp_path, qapp, monkeypatch):
    from prestus_nexus.ui.update_dialog import UpdateDialog
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    d = UpdateDialog(c, remote=str(upstream))
    d.start_check()
    _finish(d)
    monkeypatch.setattr(d, "confirm", lambda: False)               # declined: nothing happens
    d.start_update()
    assert (c / "VERSION").read_text().strip() == "v0.1.0" and d.result_ is None
    monkeypatch.setattr(d, "confirm", lambda: True)
    d.start_update()
    _finish(d)
    assert d.result_.changed and (c / "VERSION").read_text().strip() == "v0.2.0"
    assert "restart the Nexus" in d.headline.text() and "install.sh" in d.log.toPlainText()
    assert not d.update_btn.isEnabled()


def test_dialog_blocks_the_update_on_local_changes(upstream, tmp_path, qapp):
    from prestus_nexus.ui.update_dialog import UpdateDialog
    c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
    (c / "functions" / "f.m").write_text("% edit\n", encoding="utf-8")
    d = UpdateDialog(c, remote=str(upstream))
    d.start_check()
    _finish(d)
    assert d.status.update_available and not d.update_btn.isEnabled()
    assert "functions/f.m" in d.details.text()


def test_main_window_wires_the_help_entry():
    src = (Path(update.__file__).parent / "ui" / "main_window.py").read_text(encoding="utf-8")
    assert "Check for updates…" in src and "def check_for_updates" in src


# ── wording: no message needs git knowledge to act on ────────────────────────

GIT_TERMS = ("commit", "stash", "revert", "tracked", "branch", "fast-forward", "checkout",
             "detached", "HEAD", "clone", "submodule", "fetch", "tag")


def _plain(text: str) -> str:
    """The part meant for every user: 'Technical details' lines are for maintainers."""
    return "\n".join(l for l in text.splitlines() if not l.startswith("Technical details"))


@pytest.mark.parametrize("name", [*STATES, "edited", "archive", "offline"])
def test_messages_use_no_git_terms(name, upstream, tmp_path):
    remote = str(upstream)
    if name in STATES:
        c = STATES[name][0](upstream, tmp_path / "clone")
    elif name == "archive":
        c = tmp_path / "PRESTUS"
        shutil.copytree(upstream, c, ignore=shutil.ignore_patterns(".git"))
    else:
        c = clone_at(upstream, tmp_path / "clone", "v0.1.0")
        if name == "edited":
            (c / "functions" / "f.m").write_text("% edit\n", encoding="utf-8")
        else:
            remote = str(tmp_path / "nowhere")
    st = update.check(c, remote)
    res = update.apply(c, remote=remote)
    text = _plain("\n".join([st.headline(), *st.lines()[1:], *res.messages,
                              update.RUNNING_JOBS_NOTE]))
    text = text.replace("installed with git", "").replace("program git", "")
    found = [t for t in GIT_TERMS if __import__("re").search(rf"\b{t}", text)]
    assert not found, f"{found} in:\n{text}"
    if st.technical():                                             # negative control: the
        assert "commit" in st.technical()                          # git state is still shown
