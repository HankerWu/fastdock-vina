"""Recover the 65 candidate command lines from the exact runner and input generator."""

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "geometry_comparison/upstream/four_versions"))
import run_panel


class Captured(Exception):
    pass


def main(work, binary):
    result = {}

    def capture(argv, log, timeout):
        result[log.stem] = dict(
            argv=list(map(str, argv)), timeout_s=timeout, address_limit=run_panel.LIMIT
        )
        raise Captured

    run_panel.execute = capture
    for spec in run_panel.panel(work / "fixtures"):
        try:
            run_panel.dock(("candidate", binary.resolve(), work / "raw", spec))
        except Captured:
            pass
    if len(result) != 65:
        raise ValueError("incomplete reconstructed argv")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    main(**vars(parser.parse_args()))
