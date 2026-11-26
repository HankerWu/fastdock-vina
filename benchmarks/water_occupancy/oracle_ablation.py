"""After prediction, remove surely displaced external sites identified by each test crystal."""

import argparse
import itertools
import json
from pathlib import Path

from evaluate import TARGET_OF, iterate, parse_pool, read, score_components, write


def run(args):
    out = args.out
    audit = read(out / "query_site_audit.csv")
    key = lambda r: (
        r["panel"],
        r["rec"],
        r["lig"],
        str(r["seed"]),
        r["arm"],
        r["candidate"],
    )
    saved = {key(r): r for r in read(out / "poses.csv")}
    cells = []
    for name, rec, lig, seed, arm, _, models in iterate(args):
        target = TARGET_OF[name]
        sites = json.loads((out / f"{rec}_sites.json").read_text())["sites"]
        removed = [
            sites[int(r["site_index"])]
            for r in audit
            if r["target"] == target
            and r["rec"] == rec
            and r["query"] == lig
            and r["crystal_displaces"] == "1"
        ]
        scored = []
        for mid, path, index, _, rms, sign in models:
            record = saved[name, rec, lig, str(seed), arm, mid]
            subtotal = score_components(removed, *parse_pool(path)[index])["joint"]
            scored.append(
                {
                    "candidate": mid,
                    "oracle_score": float(record["joint"]) - subtotal,
                    "rmsd": rms,
                    "signed_hit": int(rms <= 2 and sign is not False),
                }
            )
        chosen = min(scored, key=lambda r: r["oracle_score"])
        cells.append(
            {
                "panel": name,
                "rec": rec,
                "lig": lig,
                "seed": seed,
                "arm": arm,
                "removed": len(removed),
                **chosen,
            }
        )
    write(out / "oracle_cells.csv", cells)
    for g, rows in itertools.groupby(
        sorted(cells, key=lambda r: (r["panel"], r["arm"])),
        lambda r: (r["panel"], r["arm"]),
    ):
        rows = list(rows)
        print(
            g,
            sum(r["rmsd"] <= 2 for r in rows),
            sum(r["signed_hit"] for r in rows),
            len(rows),
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("brd_pools", "ca_pools", "state_pools", "out"):
        p.add_argument("--" + key.replace("_", "-"), required=True, type=Path)
    run(p.parse_args())
