"""SSH_ASKPASS helper: ssh runs it for each prompt (password, one-time code, new host key)
and reads the answer from stdout. Started as `prestus-nexus --askpass <prompt>`; see
remote.askpass_env for the wrapper script ssh calls."""

from __future__ import annotations

import sys


def run(prompt: str) -> int:
    from PyQt6.QtWidgets import QApplication, QInputDialog, QLineEdit
    app = QApplication.instance() or QApplication(sys.argv[:1])  # noqa: F841 (keep a reference)
    prompt = prompt.strip() or "Password:"
    # a host-key question (yes/no/[fingerprint]) is answered in clear text
    echo = QLineEdit.EchoMode.Normal if "(yes/no" in prompt else QLineEdit.EchoMode.Password
    text, ok = QInputDialog.getText(None, "SSH login — PRESTUS Nexus", prompt, echo)
    if not ok:
        return 1
    sys.stdout.write(text + "\n")
    return 0
