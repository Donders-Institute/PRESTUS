"""Usage-statistics consent: the record telemetry_set_consent.m writes and track_usage.m reads."""

import json

from prestus_nexus import telemetry


def test_record_matches_matlab_jsonencode(tmp_path):
    f = tmp_path / ".prestus" / "telemetry.json"
    telemetry.set_consent(True, None, f, today="2026-11-14")
    assert f.read_text(encoding="utf-8") == \
        '{"opt_in":true,"decided_on":"2026-11-14","prestus_ver":"unknown"}\n'
    c = telemetry.read_consent(f)
    assert c.opted_in and c.decided_on == "2026-11-14" and c.prestus_ver == "unknown"


def test_round_trip_and_reset(tmp_path, prestus_root):
    f = tmp_path / "telemetry.json"
    assert telemetry.read_consent(f).state == "undecided"          # default: nothing sent
    telemetry.set_consent(False, prestus_root, f)
    c = telemetry.read_consent(f)
    assert c.state == "opted_out" and not c.opted_in and c.prestus_ver
    assert json.loads(f.read_text())["opt_in"] is False
    assert telemetry.reset_consent(f) and not f.exists()
    assert not telemetry.reset_consent(f)


def test_matlab_written_and_broken_records(tmp_path):
    f = tmp_path / "telemetry.json"
    f.write_text('{"opt_in":true,"decided_on":"2025-01-02","prestus_ver":"628edd0"}\n')
    assert telemetry.read_consent(f).opted_in
    f.write_text('{"opt_in": 1}')                 # isequal(1, true) holds in MATLAB
    assert telemetry.read_consent(f).opted_in
    f.write_text('{"opt_in": "true"}')            # a string is not consent
    assert telemetry.read_consent(f).state == "opted_out"
    f.write_text("{oops")
    c = telemetry.read_consent(f)
    assert c.state == "unreadable" and not c.opted_in


def test_default_location(monkeypatch, tmp_path):
    monkeypatch.setattr(telemetry.Path, "home", classmethod(lambda cls: tmp_path))
    assert telemetry.consent_file() == tmp_path / ".prestus" / "telemetry.json"


def test_env_override(monkeypatch):
    monkeypatch.setenv("PRESTUS_TELEMETRY", "OFF")
    assert telemetry.env_override() == "off"
    monkeypatch.setenv("PRESTUS_TELEMETRY", "maybe")
    assert telemetry.env_override() == ""
