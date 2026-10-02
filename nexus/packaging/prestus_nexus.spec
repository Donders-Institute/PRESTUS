# PyInstaller spec for a standalone PRESTUS Nexus app.
#
#   cd nexus
#   pip install -e ".[viewer]" pyinstaller
#   pyinstaller packaging/prestus_nexus.spec      # -> dist/PRESTUS Nexus(.app)
#
# The app does not bundle PRESTUS itself: on first start it asks for the PRESTUS
# folder (or use --prestus-root / PRESTUS_ROOT) and remembers it. MATLAB is
# found as usual (PRESTUS_MATLAB, PATH, standard install folders).
import sys
from pathlib import Path

here = Path(SPECPATH).resolve().parent          # nexus/
pkg = here / "prestus_nexus"

a = Analysis(
    [str(here / "packaging" / "run_nexus.py")],
    pathex=[str(here)],
    datas=[(str(pkg / "overlay.yaml"), "prestus_nexus"),
           (str(pkg / "run_presets.yaml"), "prestus_nexus")],
    hiddenimports=["nibabel", "matplotlib.backends.backend_qtagg"],
    excludes=["tkinter"],
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="PRESTUS Nexus",
          console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.datas, name="PRESTUS Nexus", upx=False)
if sys.platform == "darwin":
    app = BUNDLE(coll, name="PRESTUS Nexus.app", bundle_identifier="nl.donders.prestus.nexus")
