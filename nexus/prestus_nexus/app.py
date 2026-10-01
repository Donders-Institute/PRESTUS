"""Entry point: ``prestus-nexus [config.yaml]`` or ``python -m prestus_nexus``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="prestus-nexus", description=__doc__)
    parser.add_argument("config", nargs="?", help="study config YAML to open directly")
    parser.add_argument("--prestus-root", help="PRESTUS checkout (default: auto-detect)")
    args = parser.parse_args(argv)

    from PyQt6.QtWidgets import QApplication, QMessageBox

    from .equipment import EquipmentLibrary
    from .model import ConfigModel
    from .paths import base_default_for, find_prestus_root, toolbox_default
    from .schema import Schema
    from .session import Session
    from .ui.main_window import MainWindow

    app = QApplication(sys.argv[:1])
    app.setApplicationName("PRESTUS Nexus")
    session = Session()

    root = find_prestus_root(args.prestus_root or session.get("prestus_root"))
    if root is None:
        QMessageBox.critical(None, "PRESTUS Nexus",
                             "Could not find PRESTUS (config/config_default.yaml).\n"
                             "Pass --prestus-root or set PRESTUS_ROOT.")
        return 2
    session.set("prestus_root", str(root))
    from .ui import style
    style.apply(app, root, dark=bool(session.get("dark_mode", False)))

    if args.config:
        path = Path(args.config).resolve()
        model = ConfigModel.from_file(Schema(base_default_for(path, root)), path,
                                      EquipmentLibrary(root))
        session.add_recent(path)
    else:           # untitled study; File > New study… / Open… / Open recent
        model = ConfigModel(Schema(toolbox_default(root)), EquipmentLibrary(root))

    window = MainWindow(model, session, root)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
