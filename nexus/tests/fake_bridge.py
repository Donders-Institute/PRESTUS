"""A stand-in for functions/gui/nexus_bridge.m that speaks the same file protocol, so the
Nexus side of the bridge can be tested without MATLAB.

  python fake_bridge.py <bridge_dir> [--fail-start]

nexus_ping echoes its args; args.action == 'fail' answers with an error; args.progress
writes a progress message and waits a little; args.canned returns that value as result.
"""

import json
import os
import sys
import time
from pathlib import Path


def write(path: Path, data) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    os.replace(tmp, path)


def main() -> None:
    d = Path(sys.argv[1])
    if "--fail-start" in sys.argv:
        print("Error: no licence")
        sys.exit(3)
    for sub in ("req", "resp", "progress"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    write(d / "ready.json", {"pid": os.getpid(), "matlab": "fake"})
    print("[nexus_bridge] ready", flush=True)
    while not (d / "stop").exists():
        for f in sorted((d / "req").glob("*.json")):
            req = json.loads(f.read_text(encoding="utf-8"))
            f.unlink()
            args = req.get("args") or {}
            rid = req["id"]
            if args.get("progress"):
                (d / "progress" / f"{rid}.txt").write_text(args["progress"], encoding="utf-8")
                time.sleep(0.4)
                (d / "progress" / f"{rid}.txt").unlink(missing_ok=True)
            if args.get("action") == "fail":
                resp = {"id": rid, "ok": False, "error": "it failed", "identifier": "fake:fail",
                        "report": "trace", "log": "printed"}
            elif "canned" in args:
                resp = {"id": rid, "ok": True, "result": args["canned"], "log": ""}
            else:
                resp = {"id": rid, "ok": True, "result": {"echo": args, "fn": req["fn"]}, "log": "hi\n"}
            write(d / "resp" / f"{rid}.json", resp)
        time.sleep(0.02)
    print("[nexus_bridge] stop requested", flush=True)


if __name__ == "__main__":
    main()
