"""Usage-statistics consent: the record telemetry_set_consent.m writes and track_usage.m reads."""

import json
import os

import pytest

from prestus_nexus import telemetry


@pytest.fixture(scope="module")
def app():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PyQt6")
    from PyQt6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_record_matches_matlab_jsonencode(tmp_path):
    f = tmp_path / ".prestus" / "telemetry.json"
    telemetry.set_consent(True, None, f, today="2026-11-14")
    assert f.read_text(encoding="utf-8") == \
        '{"opt_in":true,"decided_on":"2026-11-14","prestus_ver":"unknown","notice":2}\n'
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


def _consent(tmp_path, text):
    f = tmp_path / "telemetry.json"
    f.write_text(text)
    return telemetry.read_consent(f)


def test_notice_version_read_as_track_usage(tmp_path):
    assert _consent(tmp_path, '{"opt_in":true}').notice == 1                  # absent: version 1
    assert _consent(tmp_path, '{"opt_in":true,"notice":2}').notice == 2
    assert _consent(tmp_path, '{"opt_in":true,"notice":"2"}').notice == 1     # not numeric


def test_status_text_for_every_state():
    """Every recorded state and override: green exactly when records are sent."""
    cases = {   # (state, env): (words, colour)
        ("undecided", ""): ("Not sharing", "amber"), ("opted_in", ""): ("Sharing", "green"),
        ("opted_out", ""): ("Not sharing", "amber"), ("unreadable", ""): ("Not sharing", "amber"),
        ("undecided", "on"): ("Sharing", "green"), ("opted_in", "off"): ("Not sharing", "amber"),
        ("opted_out", "on"): ("Not sharing", "amber"),          # a recorded opt-out wins
    }
    for (state, env), (words, colour) in cases.items():
        text, tip, role = telemetry.status(telemetry.Consent(state, "2026-11-19", notice=2), env)
        assert text.startswith(words + " usage statistics") and role == colour, (state, env, text)
        assert "Click" in tip
    text, tip, role = telemetry.status(telemetry.Consent("opted_in", "2026-01-01", notice=1), "")
    assert "notice updated" in text and role == "green"
    # declining under an old notice is not outdated: nothing is sent either way
    assert not telemetry.outdated(telemetry.Consent("opted_out", notice=1))


def test_launch_prompt_rule():
    C = telemetry.Consent
    v = telemetry.NOTICE_VERSION
    assert telemetry.should_prompt(C("undecided"), "")
    assert telemetry.should_prompt(C("undecided"), "", asked_notice=v)       # until answered
    assert not telemetry.should_prompt(C("undecided"), "on")                 # env decides
    assert not telemetry.should_prompt(C("undecided"), "off")
    assert not telemetry.should_prompt(C("opted_in", notice=v), "")
    assert not telemetry.should_prompt(C("opted_out", notice=1), "")
    assert not telemetry.should_prompt(C("unreadable"), "")                  # a decision, as MATLAB
    assert telemetry.should_prompt(C("opted_in", notice=v - 1), "", asked_notice=v - 1)
    assert not telemetry.should_prompt(C("opted_in", notice=v - 1), "", asked_notice=v)  # once
    assert not telemetry.should_prompt(C("opted_in", notice=v - 1), "off")


def test_window_shows_decision_and_prompts_once(app, schema, tmp_path, monkeypatch):
    """The status bar names the recorded decision and follows a change; the launch
    prompt opens only when should_prompt says so, and an outdated Allow only once."""
    from prestus_nexus.model import ConfigModel
    from prestus_nexus.session import Session
    from prestus_nexus.ui.main_window import MainWindow
    monkeypatch.setattr(telemetry.Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.delenv("PRESTUS_TELEMETRY", raising=False)
    shown = []
    monkeypatch.setattr(telemetry, "show_consent_dialog",
                        lambda parent, root, at_launch=False: shown.append(at_launch))
    session = Session(tmp_path / "s.yaml")
    w = MainWindow(ConfigModel(schema), session, telemetry.Path(__file__).resolve().parents[2])
    assert "Not sharing" in w.usage_status.text()
    w.ask_usage_consent()
    assert shown == [True]
    telemetry.set_consent(True, None, today="2026-11-19")
    w.refresh_usage_status()
    assert "Sharing" in w.usage_status.text() and "2026-11-19" in w.usage_status.toolTip()
    w.ask_usage_consent()
    assert shown == [True]                                     # decided: no prompt
    telemetry.consent_file().write_text('{"opt_in":true,"decided_on":"2026-01-01"}\n')
    w.ask_usage_consent()
    w.ask_usage_consent()
    assert shown == [True, True]                               # outdated Allow: asked once
    w.usage_status.linkActivated.emit("usage")
    assert shown == [True, True, False]                        # click opens the dialog
    w.close()
