"""Entry point: ``prestus-nexus [config.yaml]`` or ``python -m prestus_nexus``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["--askpass"]:                 # ssh's password prompt (remote.askpass_env)
        from .askpass import run
        return run(" ".join(argv[1:]))
    if argv[:1] == ["--update"]:                  # headless: check for / move to a release
        from .update import main as update_main
        return update_main(argv[1:])
    from .localcopy import local_main     # from a network drive: start from a local copy
    local = local_main()
    if local is not None:
        return local(argv)
    parser = argparse.ArgumentParser(prog="prestus-nexus", description=__doc__)
    parser.add_argument("config", nargs="?", help="study config YAML to open directly")
    parser.add_argument("--prestus-root", help="PRESTUS checkout (default: auto-detect)")
    parser.epilog = ("Updating PRESTUS without the window: prestus-nexus --update [--check] "
                     "[--to TAG]  (see prestus-nexus --update --help)")
    args = parser.parse_args(argv)

    from PyQt6.QtCore import QCoreApplication, Qt
    from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox

    from .equipment import EquipmentLibrary
    from .model import ConfigModel
    from .paths import base_default_for, find_prestus_root, toolbox_default
    from .schema import Schema
    from .session import Session
    from .ui.main_window import MainWindow

    # Results ▸ Report imports QtWebEngine on first use; Qt allows that after the
    # QApplication exists only with shared OpenGL contexts (set before it is created).
    QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    app = QApplication(sys.argv[:1])
    app.setApplicationName("PRESTUS Nexus")
    from .ui import errors
    errors.install_excepthook(app.activeWindow)    # an uncaught error: dialog + log, no exit
    from .ui import selectable
    selectable.install(app)                        # all label text can be selected and copied
    from .bridge import shutdown_shared
    app.aboutToQuit.connect(shutdown_shared)       # end the MATLAB session with the app
    session = Session()

    root = find_prestus_root(args.prestus_root or session.get("prestus_root"))
    while root is None:
        # Standalone builds (and checkouts elsewhere) cannot find PRESTUS by their own
        # location: ask once; the choice is remembered in the session.
        QMessageBox.information(None, "PRESTUS Nexus",
                                "Select your PRESTUS folder (the one containing config/ and "
                                "functions/). You can also pass --prestus-root or set PRESTUS_ROOT.")
        chosen = QFileDialog.getExistingDirectory(None, "PRESTUS folder", str(Path.home()))
        if not chosen:
            return 2
        root = find_prestus_root(chosen)
        if root is None:
            QMessageBox.warning(None, "PRESTUS Nexus",
                                f"{chosen} does not contain config/config_default.yaml and functions/.")
    session.set("prestus_root", str(root))
    from .ui import style
    style.apply(app, root, dark=bool(session.get("dark_mode", False)))

    model = None
    if args.config:
        path = Path(args.config).resolve()
        with errors.guarded(None, f"Opening {path.name}") as opened:   # falls back to an untitled study
            model = ConfigModel.from_file(Schema(base_default_for(path, root),
                                                 toolbox_path=toolbox_default(root)),
                                          path, EquipmentLibrary(root))
        if not opened.failed:
            session.add_recent(path)
    if model is None:   # untitled study; File > New study… / Open… / Open recent
        model = ConfigModel(Schema(toolbox_default(root), toolbox_path=toolbox_default(root)),
                            EquipmentLibrary(root))

    window = MainWindow(model, session, root)
    window.show()
    from PyQt6.QtCore import QTimer
    QTimer.singleShot(0, window.ask_usage_consent)   # once the window is on screen
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
