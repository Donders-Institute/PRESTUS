#!/usr/bin/env python3
"""Generate tests/tests_function_map.md: which tests reach each PRESTUS source file.

Internal developer documentation (not part of the docs site). For every tracked
.m file under functions/ (prestus_pipeline.m included) the map lists:

  direct tests    test classes in tests/ whose code names the file's function;
  indirect tests  test classes that reach it through the static call graph;
  levels          the level tags (prestus_test_levels.m) of those tests;
  status          direct-unit / indirect-unit only / physics only (weekly CI) /
                  non-CI levels only / none.

For nexus/prestus_nexus/**/*.py it lists the pytest files in nexus/tests that
import the module (directly, or through imports inside prestus_nexus).

Heuristic (static, no MATLAB needed):
  * Comments ('%', '%{ ... %}', text after '...') and string literals are
    removed; every remaining identifier not preceded by '.' counts as a call
    when it is the name of a functions/ file (or a tests/ helper). A string
    literal that is exactly a function name counts too, but only in files
    that dispatch dynamically (feval, str2func, eval, evalc, parfeval, batch),
    e.g. the allowed-function list of nexus_bridge. Elsewhere such strings are
    labels (log_timer('start', 'prestus_pipeline')) and are ignored.
  * Names defined as local functions or methods in a file shadow files of the
    same name and are not counted for that file.
  * Call graph edges: file G -> file F when G names F. A test reaches F when it
    names G and G reaches F within MAX_DEPTH edges. tests/fixtures/*.m and
    the non-test .m files in tests/ are intermediate nodes.
  * A test class's levels are all level tags in its class- and method-level
    TestTags; the map does not attribute a call to a single method.

Reachability is an upper bound: a reached file may never execute (branches not
taken, early errors, options off), and calls by feval of a built string,
function handles stored in data, or cellfun over name lists are missed. It is
not measured coverage. Measured line coverage of the unit level comes from
tests/coverage_unit_job.m (Cobertura XML in tests/coverage/, gitignored); when
that file exists, a copy of the map with a "line coverage" column is written to
tests/coverage/tests_function_map_coverage.md. The tracked map never carries it,
so --check gives the same answer with or without a local coverage run.

Usage:
    python tests/docs/gen_test_map.py           # rewrite tests/tests_function_map.md
    python tests/docs/gen_test_map.py --check   # exit 1 if it is out of date
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "tests"
OUT = TESTS / "tests_function_map.md"
COVERAGE_XML = TESTS / "coverage" / "unit_coverage.xml"
COVERAGE_OUT = TESTS / "coverage" / "tests_function_map_coverage.md"
NEXUS_PKG = ROOT / "nexus" / "prestus_nexus"
NEXUS_TESTS = ROOT / "nexus" / "tests"

MAX_DEPTH = 8                 # call-graph edges from a test's direct calls
MAX_LISTED = 6                # test names shown per cell before "+n more"
CI_PUSH = {"unit"}            # .github/workflows/tests.yml: matlab-unit, every push
CI_WEEKLY = {"physics_axisym"}  # matlab-physics: weekly and on demand
STATUSES = ["direct-unit", "indirect-unit only", "physics only (weekly CI)",
            "non-CI levels only", "none"]
TOP_LEVEL = "(top level)"

IDENT = re.compile(r"(?<![\w.])([A-Za-z]\w*)")
FUNC_DEF = re.compile(r"^\s*function\b\s*(?:\[[^\]]*\]\s*=\s*|[A-Za-z]\w*\s*=\s*)?([A-Za-z]\w*)", re.M)
DYNAMIC = re.compile(r"\b(feval|str2func|evalc|eval|parfeval|batch)\s*\(")
TAGS = re.compile(r"TestTags\s*=\s*\{([^}]*)\}")
QUOTED = re.compile(r"""['"]([^'"]+)['"]""")


# --------------------------------------------------------------------------- files
def tracked(*pathspecs: str) -> list[Path] | None:
    """Tracked files (git ls-files), or None when git is not available."""
    try:
        res = subprocess.run(["git", "-C", str(ROOT), "ls-files", "-z", "--", *pathspecs],
                             capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return sorted(ROOT / p for p in res.stdout.decode("utf-8").split("\0") if p)


def list_files(pathspec: str, folder: Path, pattern: str) -> list[Path]:
    files = tracked(pathspec)
    if files is None:
        files = sorted(folder.glob(pattern))
    return [f for f in files if f.is_file() and "__pycache__" not in f.parts]


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- MATLAB parsing
def strip_matlab(text: str) -> tuple[str, list[str]]:
    """Code without comments and string literals, plus the string literals."""
    code: list[str] = []
    strings: list[str] = []
    block = 0
    for line in text.splitlines():
        s = line.strip()
        if s == "%{":
            block += 1
            continue
        if s == "%}" and block:
            block -= 1
            continue
        if block:
            continue
        out: list[str] = []
        i, n = 0, len(line)
        while i < n:
            c = line[i]
            if c == "%":
                break
            if c == "." and line.startswith("...", i):
                break
            if c == '"' or (c == "'" and not (i and re.match(r"[\w)\]}.']", line[i - 1]))):
                j, buf = i + 1, []
                while j < n:
                    if line[j] == c:
                        if j + 1 < n and line[j + 1] == c:
                            buf.append(c)
                            j += 2
                            continue
                        break
                    buf.append(line[j])
                    j += 1
                strings.append("".join(buf))
                out.append(" ")
                i = j + 1
                continue
            out.append(c)
            i += 1
        code.append("".join(out))
    return "\n".join(code), strings


def references(path: Path, names: set[str], text: str | None = None,
               defined: set[str] | None = None) -> set[str]:
    """Function names that `text` (default: the file) calls; see the module docstring."""
    code, strings = strip_matlab(read(path) if text is None else text)
    found = {m.group(1) for m in IDENT.finditer(code)} & names
    if DYNAMIC.search(code):
        found |= {s for s in strings if s in names}
    found -= set(FUNC_DEF.findall(code)) if defined is None else defined
    found.discard(path.stem)
    return found


def level_tags() -> set[str]:
    """Level tags: the tag lists on the right-hand sides in prestus_test_levels.m."""
    tags = {"standalone"}
    for line in read(TESTS / "prestus_test_levels.m").splitlines():
        code = line.split("%", 1)[0]
        if "=" in code:
            tags |= set(re.findall(r"'(\w+)'", code.split("=", 1)[1]))
    return tags


BLOCK = re.compile(r"^\s*(methods|properties|events|enumeration)\b")


def test_segments(path: Path, names: set[str], levels: set[str]) -> list[tuple[set[str], set[str]]]:
    """(levels, called names) per block of a test class.

    The file is cut at each methods/properties/events block header. A
    'methods (Test, ...)' block runs at the levels in its own TestTags plus the
    class TestTags (MATLAB merges them); setup, helper, and property blocks run
    with every test of the class and get the union of the class's levels.
    """
    text = read(path)
    defined = set(FUNC_DEF.findall(strip_matlab(text)[0]))
    chunks: list[list[str]] = [[]]
    for line in text.splitlines():
        if BLOCK.match(line):
            chunks.append([])
        chunks[-1].append(line)
    header = "\n".join(chunks[0])
    class_tags = set()
    for m in TAGS.finditer(header):
        class_tags |= set(QUOTED.findall(m.group(1)))
    segs: list[tuple[set[str] | None, set[str]]] = []
    for i, chunk in enumerate(chunks):
        body = "\n".join(chunk)
        refs = references(path, names, body, defined)
        head = chunk[0] if chunk else ""
        if i and re.match(r"^\s*methods\s*\(.*\bTest\b", head):
            tags = class_tags.copy()
            for m in TAGS.finditer(head):
                tags |= set(QUOTED.findall(m.group(1)))
            segs.append((tags & levels, refs))
        else:
            segs.append((None, refs))
    test_lv = set().union(*(lv for lv, _ in segs if lv is not None)) or (class_tags & levels)
    test_lv = test_lv or {"(no level tag)"}
    return [((lv if lv else test_lv) if lv is not None else test_lv, refs) for lv, refs in segs]


# --------------------------------------------------------------------------- MATLAB map
def matlab_map():
    sources = [p for p in list_files("functions", ROOT / "functions", "**/*.m") if p.suffix == ".m"]
    helpers = [p for p in list_files("tests", TESTS, "**/*.m")
               if p.suffix == ".m" and p.parent in (TESTS, TESTS / "fixtures")
               and not (p.parent == TESTS and p.name.startswith("test_"))]
    tests = [p for p in list_files("tests", TESTS, "*.m")
             if p.parent == TESTS and p.name.startswith("test_") and p.suffix == ".m"]

    node_paths = {p.stem: p for p in sources + helpers}
    names = set(node_paths)
    graph = {name: references(path, names) for name, path in node_paths.items()}

    levels = level_tags()
    test_info = {}
    # reached[name][test] = [depth, levels via direct calls, levels via any path]
    reached: dict[str, dict[str, list]] = {n: {} for n in names}
    # via[name] = shortest chain test -> ... -> name, preferring unit-level blocks
    via: dict[str, tuple[tuple[int, int], list[str]]] = {}
    for t in tests:
        segs = test_segments(t, names, levels)
        test_info[t.stem] = {"levels": sorted(set().union(*(lv for lv, _ in segs))),
                             "direct": set().union(*(r for _, r in segs))}
        for lv, refs in segs:
            dist = {n: 0 for n in refs}
            parent: dict[str, str] = {}
            queue = deque(sorted(refs))
            while queue:
                g = queue.popleft()
                if dist[g] >= MAX_DEPTH:
                    continue
                for f in sorted(graph[g]):
                    if f not in dist:
                        dist[f] = dist[g] + 1
                        parent[f] = g
                        queue.append(f)
            for n, d in dist.items():
                hit = reached[n].setdefault(t.stem, [d, set(), set()])
                hit[0] = min(hit[0], d)
                if d == 0:
                    hit[1] |= lv
                hit[2] |= lv
                if d > 0:
                    key = (0 if CI_PUSH & lv else 1, d)
                    if n not in via or key < via[n][0]:
                        chain = [n]
                        while chain[-1] in parent:
                            chain.append(parent[chain[-1]])
                        via[n] = (key, [t.stem] + chain[::-1])

    rows = []
    for src in sources:
        name = src.stem
        hits = reached[name]
        direct = sorted(t for t, h in hits.items() if h[0] == 0)
        indirect = sorted(t for t, h in hits.items() if h[0] > 0)
        lv = sorted(set().union(*(h[2] for h in hits.values())))
        unit_direct = any(CI_PUSH & h[1] for h in hits.values())
        unit_any = any(CI_PUSH & h[2] for h in hits.values())
        weekly = any(CI_WEEKLY & h[2] for h in hits.values())
        if unit_direct:
            status = STATUSES[0]
        elif unit_any:
            status = STATUSES[1]
        elif weekly:
            status = STATUSES[2]
        elif hits:
            status = STATUSES[3]
        else:
            status = STATUSES[4]
        rel = src.relative_to(ROOT / "functions")
        folder = rel.parts[0] if len(rel.parts) > 1 else TOP_LEVEL
        rows.append({"folder": folder, "rel": rel.as_posix(), "name": name, "status": status,
                     "direct": direct, "indirect": indirect, "levels": lv,
                     "via": via[name][1] if name in via and not unit_direct else []})
    return rows, test_info


# --------------------------------------------------------------------------- Nexus map
def py_imports(path: Path, modules: set[str]) -> set[str]:
    text = read(path)
    found: set[str] = set()
    for m in re.finditer(r"^\s*from\s+prestus_nexus((?:\.\w+)*)\s+import\s+(\([^)]*\)|[^\n]*)",
                         text, re.M):
        base = m.group(1).lstrip(".")
        if base in modules:
            found.add(base)
        for item in re.findall(r"(\w+)(?:\s+as\s+\w+)?", re.sub(r"#[^\n]*", "", m.group(2))):
            sub = f"{base}.{item}" if base else item
            if sub in modules:
                found.add(sub)
    in_pkg = NEXUS_PKG in path.parents
    for m in re.finditer(r"^\s*from\s+(\.+)((?:\w+\.?)*)\s+import\s+(\([^)]*\)|[^\n]*)",
                         text, re.M):
        if not in_pkg:
            continue
        pkg = path.relative_to(NEXUS_PKG).parent.parts
        up = len(m.group(1)) - 1
        base_parts = list(pkg[: len(pkg) - up] if up else pkg)
        if m.group(2):
            base_parts += m.group(2).strip(".").split(".")
        base = ".".join(base_parts)
        if base in modules:
            found.add(base)
        for item in re.findall(r"(\w+)(?:\s+as\s+\w+)?", re.sub(r"#[^\n]*", "", m.group(3))):
            sub = f"{base}.{item}" if base else item
            if sub in modules:
                found.add(sub)
    for m in re.finditer(r"\bprestus_nexus\.((?:\w+\.)*\w+)", text):
        parts = m.group(1).split(".")
        for k in range(len(parts), 0, -1):
            if ".".join(parts[:k]) in modules:
                found.add(".".join(parts[:k]))
                break
    # importing a module runs the __init__ of every package above it
    for mod in list(found):
        parts = mod.split(".")[:-1]
        for k in range(len(parts) + 1):
            init = ".".join(parts[:k] + ["__init__"])
            if init in modules:
                found.add(init)
    if "prestus_nexus" in text and re.search(r"^\s*(from|import)\s+prestus_nexus\b", text, re.M):
        found.add("__init__")
    found &= modules
    return found


def nexus_map():
    if not NEXUS_PKG.is_dir():
        return []
    srcs = [p for p in list_files("nexus/prestus_nexus", NEXUS_PKG, "**/*.py") if p.suffix == ".py"]
    mods = {".".join(p.relative_to(NEXUS_PKG).with_suffix("").parts): p for p in srcs}
    modules = set(mods)
    graph = {m: py_imports(p, modules) - {m} for m, p in mods.items()}
    tests = [p for p in list_files("nexus/tests", NEXUS_TESTS, "test_*.py")
             if p.parent == NEXUS_TESTS and p.name.startswith("test_") and p.suffix == ".py"]
    reached: dict[str, dict[str, int]] = {m: {} for m in modules}
    matlab_only = set()
    for t in tests:
        text = read(t)
        if re.search(r"^pytestmark\s*=.*pytest\.mark\.matlab", text, re.M):
            matlab_only.add(t.stem)
        direct = py_imports(t, modules)
        dist = {m: 0 for m in direct}
        queue = deque(direct)
        while queue:
            g = queue.popleft()
            if dist[g] >= MAX_DEPTH:
                continue
            for f in graph[g]:
                if f not in dist:
                    dist[f] = dist[g] + 1
                    queue.append(f)
        for m, d in dist.items():
            reached[m][t.stem] = d
    rows = []
    for m in sorted(modules):
        hits = reached[m]
        direct = sorted(t for t, d in hits.items() if d == 0)
        indirect = sorted(t for t, d in hits.items() if d > 0)
        ci_direct = [t for t in direct if t not in matlab_only]
        ci_any = [t for t in hits if t not in matlab_only]
        status = ("direct" if ci_direct else "indirect only" if ci_any
                  else "MATLAB-marked tests only" if hits else "none")
        rows.append({"module": m, "rel": mods[m].relative_to(ROOT).as_posix(), "status": status,
                     "direct": [t + (" (matlab)" if t in matlab_only else "") for t in direct],
                     "indirect": indirect})
    return rows


# --------------------------------------------------------------------------- coverage
def load_coverage() -> dict[str, float]:
    if not COVERAGE_XML.is_file():
        return {}
    cov: dict[str, float] = {}
    for cls in ET.parse(COVERAGE_XML).getroot().iter("class"):
        fname = cls.get("filename", "").replace("\\", "/")
        rate = cls.get("line-rate")
        if fname and rate is not None:
            cov[fname] = max(float(rate), cov.get(fname, 0.0))
    return cov


def coverage_of(rel: str, cov: dict[str, float]) -> str:
    for fname, rate in cov.items():
        if fname == rel or fname.endswith("/" + rel) or fname == "functions/" + rel:
            return f"{100 * rate:.0f}%"
    return "0% (not loaded)" if cov else ""


# --------------------------------------------------------------------------- render
def cell(items: list[str], limit: int = MAX_LISTED) -> str:
    if not items:
        return "-"
    shown = ", ".join(f"`{x}`" for x in items[:limit])
    return shown + (f" +{len(items) - limit} more" if len(items) > limit else "")


def render(rows, test_info, nexus_rows, cov: dict[str, float] | None) -> str:
    L: list[str] = []
    a = L.append
    a("# Test map: how each source file is tested")
    a("")
    a("<!-- Generated by tests/docs/gen_test_map.py; do not edit by hand. -->")
    a("")
    a("Generated by `tests/docs/gen_test_map.py`; do not edit. Regenerate with "
      "`python tests/docs/gen_test_map.py`; CI (docs job) runs it with `--check`. "
      "Internal developer documentation, not part of the docs site.")
    a("")
    a("**Static reachability is an upper bound, not measured coverage.** A test listed for a file "
      "names that file's function (direct) or names a function that names it within "
      f"{MAX_DEPTH} call-graph steps (indirect). Whether the line runs depends on the options the test "
      "sets. Calls through `feval` of a built string, function handles stored in data, or `cellfun` "
      "over a list of names are not seen; identifiers that only match a function name (a variable "
      "with the same name) are counted. Measured line coverage of the unit level: run "
      "`tests/coverage_unit_job.m`, then this generator (it writes "
      "`tests/coverage/tests_function_map_coverage.md`).")
    a("")
    a("For files without a direct `unit` test, the indirect column also gives one shortest call "
      "chain (from a `unit` test block when there is one), so the reader can judge whether the "
      "chain is real.")
    a("")
    a("Status, in order of precedence:")
    a("")
    a("- **direct-unit**: a `unit` test (CI on every push) calls it by name.")
    a("- **indirect-unit only**: only reached by `unit` tests through other functions.")
    a("- **physics only (weekly CI)**: reached only by `physics_axisym` tests (CI weekly and on demand).")
    a("- **non-CI levels only**: reached only by levels CI does not run (water, head, acoustic, "
      "thermal, physics_3d, itrusst, demo_*, slurm_live, telemetry_live, standalone; "
      "telemetry_live runs in CI only on schedule).")
    a("- **none**: no test reaches it.")
    a("")
    a("Levels are attributed per `methods (Test, TestTags = ...)` block (its own plus the class "
      "tags); setup, helper, and property blocks count for every level of the class. A level "
      "shown for a file is a level at which some test block reaches it.")
    a("")

    folders = sorted({r["folder"] for r in rows}, key=lambda f: (f == TOP_LEVEL, f))
    counts = {f: {s: 0 for s in STATUSES} for f in folders}
    for r in rows:
        counts[r["folder"]][r["status"]] += 1
    total = {s: sum(counts[f][s] for f in folders) for s in STATUSES}

    a("## Summary")
    a("")
    a("| Folder | Files | " + " | ".join(STATUSES) + " |")
    a("|---|---:|" + "---:|" * len(STATUSES))
    for f in folders:
        n = sum(counts[f].values())
        a(f"| [{f}](#{anchor(f)}) | {n} | " + " | ".join(str(counts[f][s]) for s in STATUSES) + " |")
    a(f"| **all** | **{len(rows)}** | " + " | ".join(f"**{total[s]}**" for s in STATUSES) + " |")
    a("")
    if nexus_rows:
        ns = ["direct", "indirect only", "MATLAB-marked tests only", "none"]
        nc = {s: sum(1 for r in nexus_rows if r["status"] == s) for s in ns}
        a(f"Nexus modules ({len(nexus_rows)}): " + ", ".join(f"{s} {nc[s]}" for s in ns) + ".")
        a("")

    show_cov = bool(cov)
    for f in folders:
        a(f"## {f}")
        a("")
        a(" | ".join(f"{s}: {counts[f][s]}" for s in STATUSES if counts[f][s]))
        a("")
        head = "| File | Status | Direct tests | Indirect tests | Levels |"
        sep = "|---|---|---|---|---|"
        if show_cov:
            head += " Line coverage (unit) |"
            sep += "---|"
        a(head)
        a(sep)
        for r in sorted((r for r in rows if r["folder"] == f), key=lambda r: r["rel"].lower()):
            ind = "-" if not r["indirect"] else f"{len(r['indirect'])}: " + cell(r["indirect"], 3)
            if r["via"]:
                ind += "<br>via " + " → ".join(f"`{x}`" for x in r["via"])
            line = (f"| `{r['rel']}` | {r['status']} | {cell(r['direct'])} | {ind} | "
                    f"{', '.join(r['levels']) or '-'} |")
            if show_cov:
                line += f" {coverage_of(r['rel'], cov)} |"
            a(line)
        a("")

    if nexus_rows:
        a("## Nexus (nexus/prestus_nexus)")
        a("")
        a("A pytest file counts as direct when it imports the module, indirect when it imports a "
          "module that imports it (within prestus_nexus). `(matlab)` marks files with "
          "`pytestmark = pytest.mark.matlab`, which CI skips.")
        a("")
        a("| Module | Status | Direct tests | Indirect tests |")
        a("|---|---|---|---|")
        for r in nexus_rows:
            ind = "-" if not r["indirect"] else f"{len(r['indirect'])}: " + cell(r["indirect"], 3)
            a(f"| `{r['module']}` | {r['status']} | {cell(r['direct'])} | {ind} |")
        a("")

    a("## Test classes and their levels")
    a("")
    a("| Test class | Levels | Functions named directly |")
    a("|---|---|---:|")
    src_names = {r["name"] for r in rows}
    for t in sorted(test_info):
        info = test_info[t]
        a(f"| `{t}` | {', '.join(info['levels'])} | {len(info['direct'] & src_names)} |")
    return "\n".join(L).rstrip("\n") + "\n"


def anchor(heading: str) -> str:
    return re.sub(r"[^\w\- ]", "", heading.lower()).strip().replace(" ", "-")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="fail if tests/tests_function_map.md is out of date")
    args = ap.parse_args()

    rows, test_info = matlab_map()
    nexus_rows = nexus_map()
    new = render(rows, test_info, nexus_rows, None)
    if args.check:
        old = OUT.read_text(encoding="utf-8") if OUT.is_file() else ""
        if new != old:
            print("tests/tests_function_map.md is out of date; run: python tests/docs/gen_test_map.py")
            return 1
        print("tests/tests_function_map.md is up to date")
        return 0
    OUT.write_text(new, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")
    cov = load_coverage()
    if cov:
        COVERAGE_OUT.write_text(render(rows, test_info, nexus_rows, cov), encoding="utf-8")
        print(f"wrote {COVERAGE_OUT.relative_to(ROOT)} (line coverage from {COVERAGE_XML.relative_to(ROOT)})")
    total = {s: sum(1 for r in rows if r["status"] == s) for s in STATUSES}
    print("; ".join(f"{s}: {n}" for s, n in total.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
