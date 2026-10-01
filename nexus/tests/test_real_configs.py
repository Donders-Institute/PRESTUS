"""Every study config in $PRESTUS_REAL_CONFIG_DIR loads and round-trips."""

import os
from pathlib import Path

import pytest

from prestus_nexus import yamlio
from prestus_nexus.model import ConfigModel, values_equal

_dir = os.environ.get("PRESTUS_REAL_CONFIG_DIR", "")
CONFIG_DIR = Path(_dir) if _dir else None
FILES = sorted(p for p in CONFIG_DIR.glob("*.yaml") if p.name != "config_default.yaml") \
    if CONFIG_DIR and CONFIG_DIR.is_dir() else []


@pytest.mark.skipif(not FILES, reason="set PRESTUS_REAL_CONFIG_DIR to a folder of study configs")
@pytest.mark.parametrize("path", FILES, ids=lambda p: p.name)
def test_round_trip(schema, path):
    m = ConfigModel.from_file(schema, path)
    again = ConfigModel(schema)
    again.load_dict(yamlio.load(m.to_yaml()))
    assert values_equal(m.effective(), again.effective())
