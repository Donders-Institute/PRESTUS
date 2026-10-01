"""Entry point for the PRESTUS Python GUI.

DEPRECATED: superseded by the PRESTUS Nexus (nexus/, ``python -m
prestus_nexus``), which is schema-driven and round-trip tested against
load_parameters. Kept until the Nexus covers the results views.

Usage
-----
  python main.py                    # launch with defaults
  python main.py path/to/config.yaml  # pre-load a config file
  python -c "import main; main.prestus_gui(params)"  # from Python

Can also be called from MATLAB via:
  run_command({'python', 'functions/gui_py/main.py', 'config_study.yaml'}, 'Background', true);
"""

from __future__ import annotations
import os
import sys


def prestus_gui(init_params: dict | None = None):
    """Launch the PRESTUS GUI application."""
    print("Note: functions/gui_py is deprecated; use the PRESTUS Nexus "
          "(python -m prestus_nexus, see nexus/README.md).", file=sys.stderr)
    # Ensure the gui_py folder is on sys.path
    gui_dir = os.path.dirname(os.path.abspath(__file__))
    if gui_dir not in sys.path:
        sys.path.insert(0, gui_dir)

    from PyQt6.QtWidgets import QApplication
    import style

    app = QApplication.instance() or QApplication(sys.argv)
    app.setStyleSheet(style.QSS_BASE)

    from app import PrestusApp
    window = PrestusApp(init_params=init_params)
    window.show()

    # Only call exec() when we own the event loop
    if not QApplication.instance() or not QApplication.instance().property("_loop_running"):
        app.setProperty("_loop_running", True)
        sys.exit(app.exec())


if __name__ == "__main__":
    init_params = None

    if len(sys.argv) > 1:
        config_path = sys.argv[1]
        if os.path.isfile(config_path):
            try:
                import yaml
                with open(config_path) as f:
                    init_params = yaml.safe_load(f) or {}
                init_params["config_path"] = config_path
            except Exception as exc:
                print(f"Warning: could not load {config_path}: {exc}", file=sys.stderr)

    prestus_gui(init_params)
