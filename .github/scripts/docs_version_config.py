#!/usr/bin/env python3
"""Prepare a docs config for the versioned site (called by build_docs.sh).

    docs_version_config.py <mkdocs.yaml> <version> [banner]

Appends <version>/ to site_url, so the sitemap lists the pages of this version
(the version selector maps the open page to the same page in the chosen
version through the sitemaps; without it, switching lands in the wrong
version). Adds extra.version (provider mike, default main) when missing, so the
theme shows the version selector; sets extra.version_banner to [banner]; and points
theme.custom_dir at _version_overrides when the config has no overrides of its
own. Plain text edits, because the config may hold tags (!!python/name) that a
YAML loader rejects.
"""

import re
import sys

VERSION = "  version:\n    provider: mike\n    default: main"


def prepare(text: str, version: str, banner: str = "") -> str:
    text = re.sub(r"(?m)^site_url:\s*['\"]?([^'\"\s#]+?)/?['\"]?\s*(#.*)?$",
                  lambda m: "site_url: " + m.group(1) + "/" + version + "/", text, count=1)
    if not re.search(r"(?m)^extra:\s*$", text):
        text = text.rstrip("\n") + "\n\nextra:\n" + VERSION + "\n"
    elif "provider: mike" not in text:
        text = re.sub(r"(?m)^extra:\s*$", lambda m: "extra:\n" + VERSION, text, count=1)
    if banner:
        quoted = "'" + banner.replace("'", "''") + "'"
        text = re.sub(r"(?m)^extra:\s*$", lambda m: "extra:\n  version_banner: " + quoted, text, count=1)
    if not re.search(r"(?m)^  custom_dir:", text):
        text = re.sub(r"(?m)^theme:\s*$", lambda m: "theme:\n  custom_dir: _version_overrides", text, count=1)
    return text


if __name__ == "__main__":
    cfg, version = sys.argv[1], sys.argv[2]
    banner = sys.argv[3] if len(sys.argv) > 3 else ""
    with open(cfg) as f:
        text = f.read()
    with open(cfg, "w") as f:
        f.write(prepare(text, version, banner))
