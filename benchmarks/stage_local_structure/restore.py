"""Check a freshly compiled release against the start, recording actual Git and ELF identities."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
START = "666f49d2570c047470d3c281df30d18c5ebd3ce1"
RESTORED = "1c3141238511ba98e13232ffff8f2ba51efef3d1"


def tree(sha):
    return subprocess.check_output(
        ["git", "-C", ROOT, "rev-parse", sha + ":src"], text=True
    ).strip()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("baseline", "restored", "runs", "output"):
        p.add_argument(arg, type=Path)
    args = p.parse_args()
    if tree(START) != tree(RESTORED):
        raise ValueError("restored source differs from anchor")
    verifier = HERE.parent / "search_trial_panel/restore.py"
    subprocess.run(
        [
            sys.executable,
            verifier,
            args.baseline,
            args.restored,
            args.runs,
            args.output,
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    path = args.output / "verification.json"
    rec = json.loads(path.read_text())
    rec.update(source_commit=RESTORED, source_tree=tree(RESTORED), round_start=START)
    path.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2))
