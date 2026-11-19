"""Field labels with their units: no label repeats its unit's bracket ('(MNI) (MNI mm)')."""

import re

import pytest

pytest.importorskip("PyQt6")


def test_no_label_repeats_its_unit(schema):
    from prestus_nexus.ui.section import _label_html
    checked, bad = 0, []
    for node in schema.nodes():
        if not getattr(node, "unit", None):
            continue
        text = _label_html(node)
        checked += 1
        brackets = re.findall(r"\(([^()]*)\)", text)
        if len(brackets) >= 2 and brackets[-1].split()[0] == brackets[-2].strip():
            bad.append(text)
    assert checked and not bad, bad


def test_mni_labels_read_once(schema):
    from prestus_nexus.ui.section import _label_html
    nodes = {n.path: n for n in schema.nodes()}
    assert _label_html(nodes["placement.mni.focus_pos_mm"]).endswith("(MNI mm)")
    assert "(MNI) (" not in _label_html(nodes["analysis.target_mni_mm"])
