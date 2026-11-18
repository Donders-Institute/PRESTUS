"""Doc links from the Nexus: every overlay doc reference resolves to a published page and heading."""

from __future__ import annotations

from pathlib import Path

import pytest

from prestus_nexus import docs
from prestus_nexus.schema import OVERLAY_PATH

ROOT = docs.prestus_root()


@pytest.fixture(scope="module")
def qapp():
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_overlay_doc_refs_resolve():
    refs = docs.all_refs()
    assert refs, "overlay.yaml has no doc: entries"
    bad = [f"{where}: {why}" for where, ref in refs if (why := docs.check_ref(ref))]
    assert not bad, "\n".join(bad)


def test_every_section_and_tab_has_a_doc(schema):
    tabbed = {s for t in schema.tabs for s in t["sections"]}
    missing = [s for s in schema.sections if s in tabbed and not docs.section_doc(s)]
    missing += [f"tab {t['key']}" for t in schema.tabs if not docs.tab_doc(t["key"])]
    assert not missing, f"no doc: for {missing}"


def test_key_overrides_match_keys(schema):
    import fnmatch
    paths = [n.path for n in schema.nodes()]
    unused = [p for p in (docs._overlay().get("docs") or {}) if not fnmatch.filter(paths, p)]
    assert not unused, f"docs: patterns matching no key: {unused}"


def test_doc_mentions_in_help_text_resolve():
    """doc_x.md(#anchor) mentions in help texts (config_default.yaml comments, overlay values)."""
    sources = {ROOT / "config" / "config_default.yaml": False, Path(OVERLAY_PATH): True}
    bad = []
    for path, skip_comment_lines in sources.items():
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if skip_comment_lines and line.lstrip().startswith("#"):
                continue                    # the overlay's own documentation (doc_x.md#anchor)
            for m in docs.DOC_REF.finditer(line):
                why = docs.check_ref(m.group(0))
                if why:
                    bad.append(f"{path.name}:{i}: {why}")
    assert not bad, "\n".join(bad)


def test_key_doc_prefers_override_then_section(schema):
    assert docs.key_doc("grid.mode", schema) == "doc_head-model.md#orientation"
    assert docs.key_doc("grid.resolution_mm", schema) == docs.section_doc("grid")
    # a key shown on another card (group:) follows that card
    node = schema.get("grid.async_layout")
    assert node is not None and node.group == "coupling"
    assert docs.key_doc("simulation.code_type", schema) == docs.section_doc("compute")


def test_tab_and_subtab_docs():
    assert docs.tab_doc("run") == "doc_nexus_tabs.md#run-tab"
    assert docs.tab_doc("placement", "Localite").startswith("doc_nexus_calibration-localite.md")
    assert docs.tab_doc("placement", "no such sub-tab") == docs.tab_doc("placement")
    assert docs.tab_doc("no such tab") == ""


def test_urls(monkeypatch):
    site, dir_urls = docs.site_settings()
    assert site.startswith("https://") and dir_urls
    assert docs.doc_url("doc_head-model.md#orientation", offline=False) == f"{site}/doc_head-model/#orientation"
    assert docs.doc_url("doc/doc_hpc.md", offline=False) == f"{site}/doc_hpc/"
    local = docs.doc_url("doc_hpc.md", offline=True)
    assert local.startswith("file:") and local.endswith("/doc/doc_hpc.md")
    monkeypatch.setattr(docs, "site_settings", lambda: ("", True))
    assert docs.doc_url("doc_hpc.md", offline=False).startswith("file:")


def test_slugify_matches_mkdocs_rules():
    assert docs.slugify("Manual placement (`placement.manual`)") == "manual-placement-placementmanual"
    assert docs.slugify("1. Fill In the QuickStart Tab") == "1-fill-in-the-quickstart-tab"
    assert docs.slugify("Thermal Dose (CEM43) Calculation") == "thermal-dose-cem43-calculation"
    assert docs.slugify("Solver Implementations (`code_type`)") == "solver-implementations-code_type"


def test_linkify():
    out = docs.linkify("see doc_outputs_reuse.md#overwrite-behaviour and README.md")
    assert "<a href='doc_outputs_reuse.md#overwrite-behaviour'>" in out and "README.md<" not in out


def test_help_popup_html_has_links(schema):
    from prestus_nexus.ui.tooltip import html_for
    node = schema.get("io.overwrite_files")
    plain = html_for(node, "Overwrite", "always")
    linked = html_for(node, "Overwrite", "always", links=True, doc=docs.key_doc(node.path, schema))
    assert "<a " not in plain
    assert "href='doc_outputs_reuse.md#overwrite-behaviour'" in linked and "Docs" in linked


def test_help_menu(qapp):
    from PyQt6.QtWidgets import QMenu
    menu = QMenu()
    tab_action = docs.build_help_menu(menu, lambda: ("placement", "Localite"))
    texts = [a.text() for a in menu.actions() if a.text()]
    assert texts[:3] == ["Nexus overview", "Tabs and pages", "Docs for this tab"]
    assert {"Parameters reference", "Troubleshooting", "About PRESTUS Nexus"} <= set(texts)
    assert tab_action is not None and tab_action.text() == "Docs for this tab"


def test_versioned_urls(monkeypatch):
    site, _ = docs.site_settings()
    monkeypatch.delenv("PRESTUS_DOCS_VERSION", raising=False)
    assert docs.doc_url("doc_hpc.md", offline=False) == f"{site}/doc_hpc/"      # unversioned now
    monkeypatch.setenv("PRESTUS_DOCS_VERSION", "latest")
    assert docs.doc_url("doc_head-model.md#orientation", offline=False) == \
        f"{site}/latest/doc_head-model/#orientation"
    monkeypatch.setattr(docs, "DOCS_VERSION", "dev")
    monkeypatch.delenv("PRESTUS_DOCS_VERSION")
    assert docs.doc_url("README.md", offline=False) == f"{site}/dev/"
