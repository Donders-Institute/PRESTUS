"""Schema drift: overlay.yaml must stay consistent with config_default.yaml and the code."""

import fnmatch
import functools
import re

import pytest

from prestus_nexus.schema import KINDS, extract_help


@functools.lru_cache(maxsize=1)
def _functions_source(prestus_root) -> str:
    return "\n".join(p.read_text(errors="ignore")
                     for p in (prestus_root / "functions").rglob("*.m"))


def test_every_default_key_is_in_schema(schema):
    def walk(d, prefix=""):
        for k, v in d.items():
            path = f"{prefix}.{k}" if prefix else k
            assert path in schema, path
            if isinstance(v, dict):
                walk(v, path)
    walk(schema.defaults)


def test_field_rules_match_existing_keys(schema):
    paths = [n.path for n in schema.nodes()]
    stale = [pat for pat in schema.field_rules if not fnmatch.filter(paths, pat)]
    assert not stale, f"overlay 'fields' patterns match no key: {stale}"


def test_extras_are_not_in_default(schema):
    in_default = [p for p in schema.extra_rules
                  if _in(schema.defaults, p.split("."))]
    assert not in_default, f"now in config_default.yaml, drop from overlay 'extra': {in_default}"


def _in(d, parts):
    for part in parts:
        if not isinstance(d, dict) or part not in d:
            return False
        d = d[part]
    return True


@pytest.mark.parametrize("section", ["extra", "deprecated"])
def test_overlay_keys_are_read_by_prestus(schema, prestus_root, section):
    source = _functions_source(prestus_root)
    keys = schema.extra_rules if section == "extra" else schema.deprecated
    unused = [p for p in keys if not re.search(r"\b" + re.escape(p.rsplit(".", 1)[-1]) + r"\b", source)]
    assert not unused, f"overlay '{section}' keys no longer referenced in functions/: {unused}"


def test_kinds_and_choices_are_consistent(schema):
    for node in schema.leaves():
        assert node.kind in KINDS, (node.path, node.kind)
        if node.choices and node.default not in ([], "", None):
            assert node.default in node.choices, (node.path, node.default, node.choices)


def test_help_text_from_comments(schema):
    assert "voxel" in schema.get("placement.manual.coordinate_space").help
    # continuation lines are included
    assert "cpp_gpu" in schema.get("simulation.code_type").help
    # block comment above a key without inline comment
    assert "flhm_3db" in schema.get("calibration.focal_reference").help


def test_extract_help_ignores_commented_keys():
    text = ("grid:\n"
            "  a: 1 # first\n"
            "  # b: 2  # commented out\n"
            "  #        more on b\n"
            "  c: 3\n")
    helps = extract_help(text)
    assert helps == {"grid.a": "first"}


def test_integer_defaults_are_numbers(schema):
    # MATLAB has no int/float distinction: timing.pd: 0 must accept 0.06
    assert schema.get("timing.pd").kind == "number"
    assert schema.get("subject_id").kind == "int"


def test_display_rules_match_existing_keys(schema):
    import fnmatch
    keys = [n.path for n in schema.nodes()]
    for pattern in schema.display_rules:
        assert any(fnmatch.fnmatchcase(k, pattern) for k in keys), f"display rule matches nothing: {pattern}"


def test_every_section_is_on_a_tab_and_groups_exist(schema):
    tabbed = {s for t in schema.tabs for s in t["sections"]}
    for section in schema.section_order():
        assert section in tabbed, f"section {section} is not on a tab"
    for node in schema.nodes():
        if node.group and not node.group.startswith("__"):
            assert node.group in tabbed, f"{node.path}: group {node.group} is not a section on a tab"


def test_labels_are_readable(schema):
    from prestus_nexus import labels
    for node in schema.leaves():
        text = labels.label(node)
        assert text and (node.label or "_" not in text), f"{node.path}: label {text!r}"
