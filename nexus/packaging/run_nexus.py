"""Entry script for frozen builds (PyInstaller): the package's __main__ uses relative imports."""
import sys

from prestus_nexus.app import main

if __name__ == "__main__":
    sys.exit(main())
