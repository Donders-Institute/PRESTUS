"""Links from the Nexus into the PRESTUS documentation (doc/*.md, published with mkdocs).

A doc reference is written as in the help texts: ``doc_x.md`` or ``doc_x.md#anchor``. It opens
the published page (``<site_url>/doc_x/#anchor``, mkdocs' default directory URLs) in the
system browser, or the local Markdown file when the user chose offline docs or mkdocs.yaml
names no site_url. overlay.yaml gives each section and tab its page (``doc:``) and may give
single keys their own (``docs:``); the Help menu and the card and help buttons use them.
"""

from __future__ import annotations

import fnmatch
import os
import re
import unicodedata
from functools import lru_cache
from pathlib import Path

from . import __version__, yamlio
from .paths import find_prestus_root, package_dir

OVERLAY_PATH = Path(__file__).with_name("overlay.yaml")
# a doc page mention in help text: doc_x.md, doc_x.md#anchor (README.md and CHANGELOG.md are not linked)
DOC_REF = re.compile(r"\b(doc_[A-Za-z0-9_-]+\.md)(#[A-Za-z0-9_-]+)?")
_offline: bool | None = None        # set_offline; None: $PRESTUS_DOCS_OFFLINE
# Version folder of the published docs: '' for the unversioned site that is live now. When
# the site moves to mike's versioned layout (.github/workflows/mkdocs.yml deploys 'main' with
# the alias 'latest'), set this to 'latest'. $PRESTUS_DOCS_VERSION overrides it.
DOCS_VERSION = ""


def docs_version() -> str:
    """The version folder of the published docs ('' = unversioned site)."""
    v = os.environ.get("PRESTUS_DOCS_VERSION")
    return (DOCS_VERSION if v is None else v).strip().strip("/")


# ── where the docs live ──────────────────────────────────────────────────
def prestus_root() -> Path:
    return find_prestus_root() or package_dir().parents[1]


def doc_dir() -> Path:
    return prestus_root() / "doc"


def _mkdocs_text() -> str:
    p = prestus_root() / "mkdocs.yaml"
    if not p.is_file():
        p = prestus_root() / "mkdocs.yml"
    return p.read_text(encoding="utf-8") if p.is_file() else ""


@lru_cache(maxsize=1)
def site_settings() -> tuple[str, bool]:
    """(site_url without trailing '/', use_directory_urls) from mkdocs.yaml. Read with a regex:
    mkdocs.yaml may carry !!python tags a safe YAML loader refuses."""
    text = _mkdocs_text()
    m = re.search(r"^site_url:\s*['\"]?([^'\"\s#]+)", text, re.M)
    d = re.search(r"^use_directory_urls:\s*(\w+)", text, re.M)
    return (m.group(1).rstrip("/") if m else ""), (d is None or d.group(1).lower() == "true")


@lru_cache(maxsize=1)
def nav_pages() -> frozenset[str]:
    """The Markdown files the mkdocs nav lists (only those are published with a menu entry)."""
    return frozenset(re.findall(r":\s*['\"]?([\w./-]+\.md)['\"]?\s*$", _mkdocs_text(), re.M))


def set_offline(offline: bool | None) -> None:
    """Open local Markdown files instead of the published pages (None: $PRESTUS_DOCS_OFFLINE)."""
    global _offline
    _offline = offline


def is_offline() -> bool:
    if _offline is not None:
        return _offline
    return os.environ.get("PRESTUS_DOCS_OFFLINE", "").strip().lower() in ("1", "true", "yes")


# ── references ───────────────────────────────────────────────────────────
def parse_ref(ref: str) -> tuple[str, str]:
    """'doc_x.md#anchor' -> ('doc_x.md', 'anchor'); a leading 'doc/' is dropped."""
    page, _, anchor = (ref or "").strip().partition("#")
    page = page.strip()
    if page.startswith("doc/"):
        page = page[4:]
    return page, anchor.strip()


def doc_url(ref: str, offline: bool | None = None) -> str:
    """URL of a doc reference: the published page, or a file:// URL of the local Markdown."""
    page, anchor = parse_ref(ref)
    site, dir_urls = site_settings()
    if offline is None:
        offline = is_offline()
    if offline or not site:
        return (doc_dir() / page).resolve().as_uri()
    stem = page[:-3] if page.endswith(".md") else page
    if stem in ("README", "index"):
        path = "/"
    else:
        path = f"/{stem}/" if dir_urls else f"/{stem}.html"
    version = docs_version()
    return site + (f"/{version}" if version else "") + path + (f"#{anchor}" if anchor else "")


def open_doc(ref: str) -> bool:
    """Open a doc reference in the system browser (or the local file's default app)."""
    from PyQt6.QtCore import QUrl
    from PyQt6.QtGui import QDesktopServices
    if not ref:
        return False
    return QDesktopServices.openUrl(QUrl(doc_url(ref)))


def linkify(escaped_html: str, colour: str | None = None) -> str:
    """Turn doc_x.md(#anchor) mentions in already-escaped HTML into links (href = the reference);
    open them with open_doc (QLabel.linkActivated)."""
    st = f" style='color:{colour}'" if colour else ""
    return DOC_REF.sub(lambda m: f"<a href='{m.group(0)}'{st}>{m.group(0)}</a>", escaped_html)


# ── anchors (mkdocs toc slugs) ───────────────────────────────────────────
_HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*#*\s*$")
_EXPLICIT = re.compile(r"""<a\s+(?:id|name)=["']([^"']+)["']""")
_ATTR_ID = re.compile(r"\s*\{[^}]*#([\w-]+)[^}]*\}\s*$")


def slugify(text: str) -> str:
    """Python-Markdown toc's default slug: inline markup dropped, ASCII only, lowercase,
    punctuation removed, runs of spaces and hyphens to one hyphen."""
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", text)      # links and images: their text
    text = re.sub(r"<[^>]+>", "", text)
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^\w\s-]", "", text).strip().lower()
    return re.sub(r"[-\s]+", "-", text)


@lru_cache(maxsize=None)
def _anchors_cached(path: str, mtime: float) -> tuple[tuple[str, str], ...]:
    out: list[tuple[str, str]] = []
    seen: dict[str, int] = {}
    fence = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        s = line.strip()
        f = re.match(r"^(`{3,}|~{3,})", s)
        if f:
            if fence is None:
                fence = f.group(1)[0]
            elif s.startswith(fence * 3):
                fence = None
            continue
        if fence is not None:
            continue
        out += [(a, "") for a in _EXPLICIT.findall(line)]
        m = _HEADING.match(line)
        if not m:
            continue
        title = m.group(2)
        a = _ATTR_ID.search(title)
        if a:
            slug, title = a.group(1), title[:a.start()]
        else:
            slug = slugify(title)
            if slug in seen:                # toc's unique(): slug, slug_1, slug_2, …
                seen[slug] += 1
                slug = f"{slug}_{seen[slug]}"
        seen.setdefault(slug, 0)
        out.append((slug, title.strip()))
    return tuple(out)


def anchors(page: str | Path) -> dict[str, str]:
    """{anchor: heading text} of a doc page (file name in doc/ or a path)."""
    p = Path(page)
    if not p.is_absolute():
        p = doc_dir() / parse_ref(str(page))[0]
    if not p.is_file():
        return {}
    return dict(_anchors_cached(str(p), p.stat().st_mtime))


def check_ref(ref: str) -> str | None:
    """Why a doc reference is broken (missing page, page not in the nav, unknown anchor), or None."""
    page, anchor = parse_ref(ref)
    if not (doc_dir() / page).is_file():
        return f"{page}: no such file in doc/"
    if page not in nav_pages():
        return f"{page}: not in the mkdocs.yaml nav"
    if anchor and anchor not in anchors(page):
        return f"{page}#{anchor}: no such heading"
    return None


# ── which page for what ──────────────────────────────────────────────────
@lru_cache(maxsize=1)
def _overlay() -> dict:
    return yamlio.load_file(OVERLAY_PATH) or {}


def section_doc(section: str) -> str:
    """doc: of an overlay section ('' when none)."""
    return str((_overlay().get("sections", {}) or {}).get(section, {}).get("doc") or "")


def tab_doc(tab: str, subtab: str | None = None) -> str:
    """doc: of an overlay tab, or of one of its sub-tabs (subtab_docs: {label: ref})."""
    for spec in _overlay().get("tabs", []) or []:
        if spec.get("key") == tab:
            if subtab:
                sub = (spec.get("subtab_docs") or {}).get(subtab)
                if sub:
                    return str(sub)
            return str(spec.get("doc") or "")
    return ""


def key_doc(path: str, schema=None) -> str:
    """Doc page for a config key: a 'docs:' override (fnmatch, the last match wins), else the
    doc of the section whose card shows the key."""
    hit = ""
    for pattern, ref in (_overlay().get("docs", {}) or {}).items():
        if fnmatch.fnmatchcase(path, pattern):
            hit = str(ref)
    if hit:
        return hit
    section = path.split(".", 1)[0] if "." in path else "general"
    node = schema.get(path) if schema is not None else None
    if node is not None:
        section = node.group or node.section
    return section_doc(section)


def all_refs() -> list[tuple[str, str]]:
    """(where, ref) of every doc reference in overlay.yaml: sections, tabs, sub-tabs, docs."""
    ov = _overlay()
    out = [(f"sections.{k}", v["doc"]) for k, v in (ov.get("sections", {}) or {}).items() if v.get("doc")]
    for t in ov.get("tabs", []) or []:
        if t.get("doc"):
            out.append((f"tabs.{t.get('key')}", t["doc"]))
        out += [(f"tabs.{t.get('key')}.{s}", r) for s, r in (t.get("subtab_docs") or {}).items()]
    out += [(f"docs.{k}", v) for k, v in (ov.get("docs", {}) or {}).items()]
    return out


# ── Help menu ────────────────────────────────────────────────────────────
HELP_ITEMS = [
    ("Nexus overview", "doc_nexus.md"),
    ("Tabs and pages", "doc_nexus_tabs.md"),
    ("Docs for this tab", None),            # the current tab's doc:
    (None, None),
    ("Parameters reference", "doc_parameters.md"),
    ("Troubleshooting", "doc_troubleshooting.md"),
]


def version_text() -> str:
    """'PRESTUS <VERSION file> · Nexus <package version>'."""
    v = prestus_root() / "VERSION"
    tb = v.read_text(encoding="utf-8").strip() if v.is_file() else ""
    return (f"PRESTUS {tb} · " if tb else "") + f"Nexus {__version__}"


def build_help_menu(menu, current_tab_key_getter, parent=None, session=None):
    """Fill `menu` (a QMenu) with the documentation entries. current_tab_key_getter() returns
    the current tab's key, or (key, sub-tab label); session (optional) keeps the offline choice.
    Returns the 'Docs for this tab' action (to give it a shortcut, e.g. F1)."""
    from PyQt6.QtGui import QAction

    if session is not None:
        set_offline(bool(session.get("docs_offline", False)) or None)

    def this_tab() -> str:
        cur = current_tab_key_getter() if current_tab_key_getter else None
        key, sub = (cur if isinstance(cur, tuple) else (cur, None))
        if not key or key == "quickstart":
            return "doc_nexus_tabs.md#quickstart-tab"
        return tab_doc(key, sub) or "doc_nexus_tabs.md"

    tab_action = None
    for text, ref in HELP_ITEMS:
        if text is None:
            menu.addSeparator()
            continue
        action = QAction(text, parent or menu)
        if ref is None:
            action.triggered.connect(lambda _=False: open_doc(this_tab()))
            action.setToolTip("Open the documentation of the tab shown")
            tab_action = action
        else:
            action.triggered.connect(lambda _=False, r=ref: open_doc(r))
        menu.addAction(action)
    menu.addSeparator()
    local = QAction("Use local documentation files", parent or menu, checkable=True)
    local.setToolTip("Open the Markdown files in doc/ instead of the published pages (offline)")
    local.setChecked(is_offline())

    def set_local(on: bool) -> None:
        set_offline(on)
        if session is not None:
            session.set("docs_offline", bool(on))
            session.save()
    local.toggled.connect(set_local)
    menu.addAction(local)
    menu.addSeparator()
    about = QAction("About PRESTUS Nexus", parent or menu)
    about.triggered.connect(lambda _=False: show_about(parent))
    menu.addAction(about)
    return tab_action


def show_about(parent=None) -> None:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QMessageBox
    site, _ = site_settings()
    link = f"<a href='{site}/doc_nexus/'>{site}/doc_nexus/</a>" if site else "doc/doc_nexus.md"
    box = QMessageBox(parent)
    box.setWindowTitle("About PRESTUS Nexus")
    box.setTextFormat(Qt.TextFormat.RichText)
    box.setText(f"<b>PRESTUS Nexus</b><br>{version_text()}<br><br>"
                "Builds and checks PRESTUS simulation configs before they run in MATLAB.<br>"
                f"Documentation: {link}")
    box.exec()
