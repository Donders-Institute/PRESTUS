"""Choice icons and figures: overlay entries match keys, choices and drawings; all render."""

import os

import pytest

from prestus_nexus.ui import choice_icons


def test_choice_icon_entries_are_valid(schema):
    for section, rules, names in (("choice_icons", schema.choice_icon_rules, choice_icons.GLYPHS),
                                  ("choice_figures", schema.choice_figure_rules, choice_icons.FIGURES)):
        for path, per_choice in rules.items():
            node = schema.get(path)
            assert node is not None and node.choices, f"{section}: {path} is not a choice key"
            for choice, value in per_choice.items():
                assert choice in node.choices, f"{section}: {path} has no choice {choice!r}"
                for name in [value] if isinstance(value, str) else value:
                    assert name in names, f"{section}: {path}.{choice}: unknown drawing {name!r}"


def test_every_drawing_is_used(schema):
    used = {n for rules in schema.choice_icon_rules.values() for n in rules.values()}
    used |= {n for rules in schema.choice_figure_rules.values() for ns in rules.values() for n in ns}
    # tissue glyphs head the medium_properties table columns (ui/medium_table.py)
    used |= {choice_icons.tissue_glyph(t.name) for t in schema.get("medium_properties").children
             if t.is_group} - {None}
    assert set(choice_icons.GLYPHS) | set(choice_icons.FIGURES) <= used


def test_nodes_carry_icons(schema):
    node = schema.get("simulation.transducer_coupling")
    assert node.choice_icons == {"coherent": "coupling_coherent", "async": "coupling_async"}
    assert len(node.choice_figures["coherent"]) == 2
    assert schema.get("headmodel.smooth_method").choice_icons["off"] == "smooth_off"


def test_roles_become_attributes():
    body = choice_icons.qt('<path class="a" d="M0 0"/><circle class="af soft" r="1"/>')
    assert "class=" not in body
    assert 'stroke="accentColor"' in body and 'fill-opacity="0.35"' in body


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _inked(pix) -> int:
    img = pix.toImage()
    return sum(img.pixelColor(x, y).alpha() > 0
               for x in range(0, img.width(), 2) for y in range(0, img.height(), 2))


def test_all_drawings_render(app):
    from PyQt6.QtGui import QIcon
    from prestus_nexus.ui import icons
    for name in choice_icons.GLYPHS:
        pix = icons.icon(name, size=24).pixmap(24, 24, QIcon.Mode.Normal)
        assert _inked(pix) > 10, name
    for name in choice_icons.FIGURES:
        assert _inked(icons.figure(name, 160)) > 200, name
    assert "accentColor" not in icons.svg("coupling_coherent", "#000000", "#1f6fb2")
    assert "#1f6fb2" in icons.svg("coupling_coherent", "#000000", "#1f6fb2")


def test_choice_editor_shows_icons_and_figures(app, schema):
    from prestus_nexus.ui.fields import ChoiceEditor
    ed = ChoiceEditor(schema.get("simulation.transducer_coupling"))
    ed.set_value("coherent", "coherent")
    assert all(not ed.combo.itemIcon(i).isNull() for i in range(ed.combo.count()))
    figures = [ed.figures.itemAt(i).widget() for i in range(ed.figures.count())]
    assert len([w for w in figures if w is not None]) == 2
    ed.combo.setCurrentIndex(ed._items.index("async"))
    assert len([ed.figures.itemAt(i).widget() for i in range(ed.figures.count())
                if ed.figures.itemAt(i).widget() is not None]) == 2

    plain = ChoiceEditor(schema.get("simulation.precision"))     # no icons configured
    plain.set_value("single", "single")
    assert all(plain.combo.itemIcon(i).isNull() for i in range(plain.combo.count()))
    assert plain.figures.count() == 0


# ── Hover help (ui/tooltip.py) ─────────────────────────────────────────

def test_help_parse_reflows_and_lists_choices():
    from prestus_nexus.ui.tooltip import parse
    text = ("How they interact:\n'coherent'  — fields are superposed\n(models interference)\n"
            "'async'     — intensities\nare summed")
    assert parse(text, ["coherent", "async"]) == [
        ("p", "How they interact:"),
        ("items", [("coherent", "fields are superposed (models interference)", False),
                   ("async", "intensities are summed", False)])]


def test_help_parse_inline_lists_and_markers():
    from prestus_nexus.ui.tooltip import parse
    blocks = parse("Solver: 'kwave' (default) | 'jwave' (experimental)", ["kwave", "jwave"])
    assert blocks == [("p", "Solver:"), ("items", [("kwave", "", True), ("jwave", "experimental", False)])]
    # a quoted value that is not a choice stays in the paragraph
    assert parse("e.g.\n'80 mm' in a table", ["peak"])[0][0] == "p"


def test_tooltip_escapes_and_marks_default(app, schema):
    from prestus_nexus.ui.tooltip import html_for
    node = schema.get("grid.mode")
    tip = html_for(node, "Grid orientation", "focal_axis")
    assert "norm &lt; grid.focal_axis_min_norm" in tip           # '<' in help survives
    assert tip.count("(default)") == 1 and "data:image/svg+xml;utf8," in tip and "base64" not in tip


def test_medium_selector_has_icons_and_formatted_help(app, schema):
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.ui.medium_selector import MediumSelector
    sel = MediumSelector(ConfigModel(schema))
    assert all(not sel.combo.itemIcon(i).isNull() for i in range(sel.combo.count()))
    assert "simulation.medium" in sel.combo.toolTip() and "<table" in sel.combo.toolTip()


def test_placement_subtabs_carry_mode_icons(app, schema, tmp_path):
    from pathlib import Path
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import LOCALITE, PLANTUS, TUS_ENTRY, MainWindow
    w = MainWindow(ConfigModel(schema), Session(tmp_path / "s.yaml"), Path(__file__).resolve().parents[2])
    for key in (LOCALITE, PLANTUS, TUS_ENTRY):
        page = w.pages[key]
        bar = page.parentWidget()
        while bar is not None and not hasattr(bar, "indexOf"):
            bar = bar.parentWidget()
        tabs = bar.parentWidget() if bar is not None and not hasattr(bar, "tabIcon") else bar
        assert not tabs.tabIcon(tabs.indexOf(page)).isNull(), key
    w.close()
