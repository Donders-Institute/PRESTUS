"""Checking for a new PRESTUS release and moving the checkout to it.

PRESTUS is installed as a git clone, so a release is a tag (``v0.7.0``) and updating means
moving the clone to that tag. The update never discards work:

* changed tracked files stop it (nothing is stashed or overwritten);
* on a branch it only fast-forwards (the release must contain every local commit);
* a clone that sits on a tag (detached HEAD) is switched to the new tag.

Afterwards the submodules that are installed follow the release's pins, and the result
says whether the Nexus environment changed (then the installer has to run again).

No Qt and no other module of this package: only the standard library. It therefore runs

* from the Nexus (*Help ▸ Check for updates…*, ``ui/update_dialog.py``),
* headless as ``prestus-nexus --update [--check] [--to v0.7.1]``, and
* as a plain script, ``python3 nexus/prestus_nexus/update.py --check``, where no Nexus
  environment exists (e.g. on a compute cluster).

Nothing contacts the network unless one of these is asked for.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

REPO_URL = "https://github.com/Donders-Institute/PRESTUS.git"
_RELEASE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")
# the files that define the Nexus environment: a change means the installer must run again
ENV_FILES = ("nexus/environment.yml", "nexus/pyproject.toml")
# exit codes of the command line
EXIT_OK, EXIT_FAILED, EXIT_AVAILABLE = 0, 1, 10

RUNNING_JOBS_NOTE = ("Simulations that are still running or waiting in a queue would switch "
                     "to the new version part-way through. Wait until they have finished.")
DOC_HINT = "See 'Updating PRESTUS' in the Nexus documentation (doc/doc_nexus.md)."


def _names(files: list[str], limit: int = 5) -> str:
    return ", ".join(files[:limit]) + (f", and {len(files) - limit} more" if len(files) > limit else "")


def edited_files_text(files: list[str]) -> str:
    """Why edited PRESTUS files stop an update, and what to do, without git terms."""
    n = len(files)
    return (f"{n} PRESTUS file{' was' if n == 1 else 's were'} edited in this folder "
            f"({_names(files)}). An update would overwrite {'it' if n == 1 else 'them'}, so "
            "PRESTUS is not updated while the edits are there. Copy your edits somewhere else "
            "if you need them, then restore the original files and check again.")


def release_key(tag: str) -> tuple[int, int, int] | None:
    """(major, minor, patch) of a release tag ``vX.Y.Z``; None for anything else."""
    m = _RELEASE.match(tag.strip())
    return tuple(int(x) for x in m.groups()) if m else None


def run_git(root: Path | None, *args: str, timeout: float = 600) -> tuple[int, str]:
    """Run git (in `root`); returns (exit code, output). 127: git is not installed.
    git never asks for a password here (it fails instead), so a headless call cannot hang."""
    cmd = ["git"] + (["-C", str(root)] if root is not None else []) + list(args)
    env = dict(os.environ, GIT_TERMINAL_PROMPT="0", LC_ALL="C")
    kwargs = {"creationflags": getattr(subprocess, "CREATE_NO_WINDOW", 0)} if os.name == "nt" else {}
    try:
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", timeout=timeout, env=env,
                             stdin=subprocess.DEVNULL, check=False, **kwargs)
    except FileNotFoundError:
        return 127, "git was not found"
    except subprocess.TimeoutExpired:
        return 124, f"git {args[0]} did not finish within {timeout:.0f} s"
    return out.returncode, out.stdout.strip()


def is_clone(root: Path) -> bool:
    """True when `root` is the top folder of a git clone (not a downloaded archive, and not
    a folder that merely sits inside some other repository)."""
    rc, top = run_git(root, "rev-parse", "--show-toplevel")
    if rc != 0 or not top:
        return False
    try:
        return os.path.samefile(top, root)
    except OSError:
        return False


def remote_releases(remote: str = REPO_URL, root: Path | None = None) -> dict[str, str]:
    """{release tag: commit hash} of `remote`. Raises RuntimeError when it cannot be read."""
    rc, out = run_git(root, "ls-remote", "--tags", remote, "v*", timeout=120)
    if rc != 0:
        if rc == 127:
            raise RuntimeError("The program git was not found on this computer. PRESTUS "
                               "updates need it: install git (https://git-scm.com), then "
                               "check again.")
        raise RuntimeError("The list of releases could not be read. Check the internet "
                           f"connection and try again.\nTechnical details: {remote}: {out}")
    tags: dict[str, str] = {}
    for line in out.splitlines():
        sha, _, ref = line.partition("\t")
        name = ref.removeprefix("refs/tags/")
        peeled = name.endswith("^{}")          # an annotated tag: this line names its commit
        name = name.removesuffix("^{}")
        if release_key(name) and (peeled or name not in tags):
            tags[name] = sha.strip()
    return tags


def latest_release(tags: dict[str, str]) -> str | None:
    return max(tags, key=release_key) if tags else None


@dataclass
class Status:
    """Where the checkout stands relative to the latest release."""
    root: Path
    clone: bool = False
    version: str = ""               # the VERSION file
    head: str = ""                  # short commit hash
    branch: str | None = None       # None: detached HEAD (e.g. on a release tag)
    changed: list[str] = field(default_factory=list)   # changed tracked files
    latest: str | None = None       # latest release tag of the remote
    relation: str = "unknown"       # current | behind | ahead | diverged | unknown
    error: str = ""

    @property
    def update_available(self) -> bool:
        return self.relation == "behind"

    @property
    def can_update(self) -> bool:
        """A newer release exists and nothing stands in the way of installing it."""
        return self.update_available and not self.changed

    def headline(self) -> str:
        """The outcome in one short sentence."""
        if self.error:
            return "The check could not be completed"
        if self.update_available and self.changed:
            return f"{self.latest} is available, but edited files block the update"
        return {"behind": f"{self.latest} is available",
                "current": "You have the latest release",
                "ahead": "You have a development version, newer than the latest release",
                "diverged": "This PRESTUS folder cannot be updated automatically",
                }.get(self.relation, "The check could not be completed")

    def lines(self) -> list[str]:
        """What was found and what to do, in plain sentences (command line and dialog)."""
        out = [f"PRESTUS folder: {self.root}",
               f"Installed version: {self.version or 'unknown'}"]
        if self.latest:
            out.append(f"Latest release: {self.latest}")
        if self.error:
            return out + [self.error]
        out.append({
            "current": "Nothing to update.",
            "ahead": f"This folder holds a development version that already includes "
                     f"{self.latest}. Nothing to update.",
            "behind": "A newer release than the installed version exists.",
            "diverged": f"This folder holds a version with changes of its own that are not "
                        f"part of {self.latest} (a development version, or changes someone "
                        "saved here). Updating could lose them, so nothing is changed. Ask "
                        "whoever maintains this PRESTUS folder. " + DOC_HINT,
        }.get(self.relation, "The installed version could not be compared with the release."))
        if self.changed and self.relation in ("behind", "diverged"):
            out.append(edited_files_text(self.changed))
        return out

    def technical(self) -> str:
        """The git state, for whoever maintains the folder ('' for a folder without git)."""
        if not self.head:
            return ""
        where = f"branch {self.branch}" if self.branch else "detached HEAD"
        return (f"git: commit {self.head}, {where}"
                + (f", {len(self.changed)} tracked file(s) modified" if self.changed else ""))


def read_version(root: Path) -> str:
    v = Path(root) / "VERSION"
    try:
        return v.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def changed_files(root: Path) -> list[str]:
    """Tracked files with uncommitted changes (untracked files and submodule states do not
    count: git itself refuses when an update would overwrite one)."""
    rc, out = run_git(root, "status", "--porcelain", "--untracked-files=no",
                      "--ignore-submodules=all")
    # 'XY path' per line; run_git strips the output, so the first line may have lost a
    # leading blank of its two status letters: split instead of cutting at a fixed column
    return [line.split(None, 1)[1] for line in out.splitlines()
            if len(line.split(None, 1)) == 2] if rc == 0 else []


def _relation(root: Path, commit: str) -> str:
    """How HEAD relates to `commit`. A commit this clone does not have yet is ahead of it
    for all we know: 'behind' (the update fetches it and then checks again)."""
    if run_git(root, "cat-file", "-e", f"{commit}^{{commit}}")[0] != 0:
        return "behind"
    rc, head = run_git(root, "rev-parse", "HEAD")
    if rc != 0:
        return "unknown"
    if head == commit:
        return "current"
    if run_git(root, "merge-base", "--is-ancestor", commit, "HEAD")[0] == 0:
        return "ahead"
    if run_git(root, "merge-base", "--is-ancestor", "HEAD", commit)[0] == 0:
        return "behind"
    return "diverged"


def check(root: Path, remote: str = REPO_URL) -> Status:
    """Compare the checkout with the latest release of `remote` (reads the remote's tag
    list; changes nothing)."""
    root = Path(root)
    st = Status(root=root, version=read_version(root))
    st.clone = is_clone(root)
    if st.clone:
        st.head = run_git(root, "rev-parse", "--short", "HEAD")[1]
        rc, branch = run_git(root, "symbolic-ref", "--quiet", "--short", "HEAD")
        st.branch = branch if rc == 0 and branch else None
        st.changed = changed_files(root)
    try:
        tags = remote_releases(remote, root if st.clone else None)
    except RuntimeError as err:
        st.error = str(err)
        return st
    st.latest = latest_release(tags)
    if st.latest is None:
        st.error = f"No release was found at {remote}."
    elif not st.clone:
        st.error = not_installed_with_git(remote, st.latest)
    else:
        st.relation = _relation(root, tags[st.latest])
    return st


def not_installed_with_git(remote: str, latest: str | None = None) -> str:
    return ("This PRESTUS folder was downloaded as an archive (not installed with git), so "
            f"it cannot be updated in place. Download {latest or 'the new release'} from "
            f"{remote.removesuffix('.git')}/releases and use that folder instead. " + DOC_HINT)


@dataclass
class Result:
    ok: bool
    changed: bool = False           # the checkout moved
    target: str = ""
    old: str = ""                   # commit before (short)
    new: str = ""                   # commit after (short)
    env_changed: bool = False       # the Nexus environment files differ: run the installer
    messages: list[str] = field(default_factory=list)


def _installed_submodules(root: Path) -> list[str]:
    """Paths of the submodules that are checked out (those not installed stay that way;
    prestus_setup fetches missing required ones)."""
    rc, out = run_git(root, "submodule", "status")
    paths = []
    for line in out.splitlines() if rc == 0 else []:
        m = re.match(r"^\s*([+U-]?)[0-9a-f]+\s+(\S+)", line)   # '-': not checked out
        if m and m.group(1) != "-":
            paths.append(m.group(2))
    return paths


def apply(root: Path, target: str | None = None, remote: str = REPO_URL,
          log: Callable[[str], None] = lambda _msg: None) -> Result:
    """Move the checkout to release `target` (default: the latest). See the module text for
    what it refuses. `log` receives one line per step."""
    root = Path(root)
    res = Result(ok=False, target=target or "")

    def stop(msg: str) -> Result:
        res.messages.append(msg)
        log(msg)
        return res

    if not is_clone(root):
        return stop(not_installed_with_git(remote))
    changed = changed_files(root)
    if changed:
        return stop(edited_files_text(changed) + " Nothing was changed.")
    try:
        tags = remote_releases(remote, root)
    except RuntimeError as err:
        return stop(str(err))
    target = target or latest_release(tags) or ""
    res.target = target
    if target not in tags:
        known = ", ".join(sorted(tags, key=release_key)[-5:]) or "none"
        return stop(f"'{target}' is not a PRESTUS release. The latest releases are: {known}.")

    log(f"Downloading {target}")
    rc, out = run_git(root, "fetch", "--no-tags", remote, f"refs/tags/{target}:refs/tags/{target}")
    if rc != 0:
        return stop(f"{target} could not be downloaded. Check the internet connection and "
                    f"try again; nothing was changed.\nTechnical details: {out}")
    old = run_git(root, "rev-parse", "HEAD")[1]
    new = run_git(root, "rev-parse", f"refs/tags/{target}^{{commit}}")[1]
    res.old, res.new = old[:8], new[:8]
    rc, branch = run_git(root, "symbolic-ref", "--quiet", "--short", "HEAD")
    branch = branch if rc == 0 and branch else None

    if old == new:
        res.ok = True
        return stop(f"PRESTUS is already {target}. Nothing to do.")
    if branch:
        if run_git(root, "merge-base", "--is-ancestor", new, old)[0] == 0:
            res.ok = True
            return stop(f"This folder holds a development version that already includes "
                        f"{target}. Nothing to do.")
        if run_git(root, "merge-base", "--is-ancestor", old, new)[0] != 0:
            return stop(f"This folder holds a version with changes of its own that are not "
                        f"part of {target}. Updating could lose them, so nothing was changed. "
                        "Ask whoever maintains this PRESTUS folder. " + DOC_HINT + "\n"
                        f"Technical details: branch {branch} cannot be fast-forwarded; to use "
                        f"the release as it is: git -C \"{root}\" checkout {target}")
        log(f"Installing {target}")
        rc, out = run_git(root, "merge", "--ff-only", f"refs/tags/{target}")
    else:
        log(f"Installing {target}")
        rc, out = run_git(root, "-c", "advice.detachedHead=false", "checkout", "--detach",
                          f"refs/tags/{target}")
    if rc != 0:
        return stop(f"{target} could not be installed; nothing was changed.\n"
                    f"Technical details: {out}")
    res.changed = True

    subs = _installed_submodules(root)
    if subs:
        log("Updating the bundled toolboxes: " + ", ".join(Path(p).name for p in subs))
        rc, out = run_git(root, "submodule", "update", "--init", "--", *subs)
        if rc != 0:
            res.messages.append("PRESTUS was updated, but its bundled toolboxes (such as "
                                "k-Wave) were not. Run prestus_setup in MATLAB to finish.\n"
                                f"Technical details: {out}")
    res.env_changed = run_git(root, "diff", "--quiet", old, new, "--", *ENV_FILES)[0] == 1
    res.ok = True
    res.messages.insert(0, f"PRESTUS is now on {target}.")
    steps = ["Close the Nexus and MATLAB and start them again, so they use the new version."]
    if res.env_changed:
        steps.append("Before starting the Nexus, run its installer once more: nexus/install.sh "
                     "(Windows: nexus\\install.bat). The new version needs other packages.")
    steps.append("Studies with their own copy of config_default.yaml keep their old default "
                 "settings. To move one to the new defaults, open it and choose "
                 "Tools ▸ Update project defaults…")
    res.messages.append("What to do now:")
    res.messages += [f"  {i}. {step}" for i, step in enumerate(steps, 1)]
    res.messages.append(f"Technical details: commit {res.new}, before {res.old}.")
    for msg in res.messages:
        log(msg)
    return res


def default_root() -> Path:
    """The PRESTUS folder this file belongs to (also when the Nexus runs from its local copy)."""
    source = os.environ.get("PRESTUS_NEXUS_SOURCE")
    here = Path(source) if source else Path(__file__).resolve().parent
    return Path(os.environ.get("PRESTUS_ROOT") or here.parents[1])


def main(argv: list[str] | None = None) -> int:
    """``prestus-nexus --update``: exit 0 = done or nothing to do, 1 = failed or refused;
    with --check, 10 = a newer release exists."""
    parser = argparse.ArgumentParser(
        prog="prestus-nexus --update",
        description="Update this PRESTUS folder to the latest release. Nothing is changed "
                    "when PRESTUS files were edited in the folder, or when it holds a version "
                    "with changes of its own. " + RUNNING_JOBS_NOTE)
    parser.add_argument("--check", action="store_true",
                        help="only report the installed and the latest release (exit 10 when "
                             "a newer one exists)")
    parser.add_argument("--to", metavar="TAG", help="release to move to (default: the latest)")
    parser.add_argument("--prestus-root", help="PRESTUS checkout (default: the one this file is in)")
    parser.add_argument("--remote", default=REPO_URL, help=f"repository to read (default: {REPO_URL})")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    root = Path(args.prestus_root).expanduser() if args.prestus_root else default_root()
    if not (root / "config" / "config_default.yaml").is_file():
        print(f"{root} is not a PRESTUS folder (no config/config_default.yaml).", file=sys.stderr)
        return EXIT_FAILED
    if args.check:
        st = check(root, args.remote)
        print("\n".join([st.headline() + "."] + st.lines() + [t for t in [st.technical()] if t]))
        if st.update_available:
            if st.can_update:
                print("To update: prestus-nexus --update\n" + RUNNING_JOBS_NOTE)
            return EXIT_AVAILABLE
        return EXIT_FAILED if st.error else EXIT_OK
    res = apply(root, args.to, args.remote, log=print)
    return EXIT_OK if res.ok else EXIT_FAILED


if __name__ == "__main__":
    sys.exit(main())
