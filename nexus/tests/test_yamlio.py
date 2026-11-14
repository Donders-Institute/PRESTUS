"""YAML I/O must agree with PRESTUS' MATLAB loader (SnakeYAML 1.30, YAML 1.1)."""

import math

import pytest

from prestus_nexus import yamlio


@pytest.mark.parametrize("text, expected", [
    ("2.17e-3", 0.00217),      # config_default.yaml style: no exponent sign
    ("1e5", 100000.0),
    ("1.5E+3", 1500.0),
    (".5", 0.5),
    ("30", 30),
    ("yes", True),              # YAML 1.1 boolean
    ("off", False),
    ("~", None),
    ("'04:00:00'", "04:00:00"),
])
def test_load_matches_snakeyaml(text, expected):
    assert yamlio.load(f"v: {text}")["v"] == expected


@pytest.mark.parametrize("value", [
    "yes", "no", "on", "04:00:00", "sub-%1$03d_T1w.nii*", "it's", "", "~", "1e5", "#x",
    0.00217, 1e-5, 1e20, 300e3, -0.0, 12, True, False, None, math.inf,
])
def test_scalar_round_trip(value):
    assert yamlio.load(yamlio.dump({"v": value}))["v"] == value


def test_strings_are_always_quoted():
    assert yamlio.format_scalar("yes") == "'yes'"
    assert yamlio.format_scalar("it's") == "'it''s'"


def test_floats_keep_a_dot():
    # SnakeYAML's float rule needs '.' or an exponent; keep both parsers in step
    assert yamlio.format_scalar(1e-5) == "1.0e-05"
    assert yamlio.format_scalar(1e20) == "1.0e+20"
    assert yamlio.format_scalar(300000.0) == "300000.0"


def test_nested_structures_round_trip():
    data = {
        "grid": {"default_dims": [144, 144, 400], "pml_size": "auto"},
        "transducer": [
            {"name": "left", "annular": {"elem_phase_deg": [[0, 10.5], [20, 30]]}},
            {"trans_pos": [1, 2, 3]},
        ],
        "startup": {"modules_to_load": ["matlab/R2023b", "fsl"]},
        "empty": {}, "empty_list": [],
    }
    assert yamlio.load(yamlio.dump(data)) == data


def test_header_is_commented():
    text = yamlio.dump({"a": 1}, header="line one\nline two")
    assert text.startswith("# line one\n# line two\n\na: 1")


def test_load_file_cached_rereads_a_changed_file(tmp_path, monkeypatch):
    import os
    path = tmp_path / "entry.yaml"
    path.write_text("a: 1\n")
    assert yamlio.load_file_cached(path) == {"a": 1}
    calls = []
    real = yamlio.load_file
    monkeypatch.setattr(yamlio, "load_file", lambda p: calls.append(p) or real(p))
    assert yamlio.load_file_cached(path) == {"a": 1}
    assert calls == []                                   # unchanged: not read again
    path.write_text("a: 22\n")
    st = path.stat()
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns + 10**9))
    assert yamlio.load_file_cached(path) == {"a": 22}
    assert calls == [path]
