"""Replay one prespecified seed per receptor/compound with the fork's default AD4 mode."""

import argparse
import itertools
import json
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from run_official import COMPOUNDS, RECEPTORS, SEEDS, run_one, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "ligands", "archive", "binary", "report"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    work, lig, archive, binary = [
        getattr(a, k).resolve() for k in ("work", "ligands", "archive", "binary")
    ]
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        (out / "pools").mkdir()
        (out / "maps").symlink_to(archive / "maps", target_is_directory=True)
        with ThreadPoolExecutor(3) as exe:
            jobs = itertools.product(RECEPTORS, COMPOUNDS, SEEDS[:1], ["ad4"])
            got = list(exe.map(lambda j: run_one(work, lig, out, binary, j), jobs))
        pairs = []
        for row in got:
            stem = f"{row['compound']}_{row['receiver']}_{row['seed']}_ad4"
            official = archive / "pools" / (stem + ".pdbqt")
            exact = sha(official) == row["sha256"]
            pairs.append(
                dict(
                    cell=stem,
                    byte_identical=exact,
                    official_sha256=sha(official),
                    fork_sha256=row["sha256"],
                )
            )
        if not all(r["byte_identical"] for r in pairs):
            raise ValueError("AD4 pose file differs between binaries")
    a.report.write_text(
        json.dumps({"fork_binary_sha256": sha(binary), "pools": pairs}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
