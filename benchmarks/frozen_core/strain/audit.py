"""Cross-check replacement rankings independently from the raw PDBQT INTER field."""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def audit(work, affinity):
    pools = defaultdict(list)
    with (work / "modes.csv").open() as f:
        for r in csv.DictReader(f):
            q = r["compound" if affinity else "query"]
            group = r["cap" if affinity else "panel"]
            pools[
                group, q, r["receptor" if affinity else "receiver"], r["seed"]
            ].append(r)
    maxerr = nchanges = classchanges = nmodes = 0
    for (group, q, r, s), modes in pools.items():
        if affinity:
            path = work / group / "pools" / f"{r}_{q}_{s}_armijo.pdbqt"
        else:
            arm = "native" if group.startswith("matched_") else "armijo"
            folder = work / group
            if group in ("high", "low"):
                folder = folder / "poses"
            path = folder / f"{q}-{r}-{s}-{arm}.pdbqt"
        inter = list(
            map(
                float,
                re.findall(
                    r"^REMARK INTER:\s+(-?[\d.]+)", path.read_text(), re.MULTILINE
                ),
            )
        )
        alts = []
        for val, m in zip(inter, modes, strict=True):
            scale = float(m["flexibility_scale"])
            new = scale * (val + float(m["fixed_mmff"]) - float(m["free_mmff"]))
            maxerr = max(maxerr, abs(new - float(m["full_unit"])))
            alts.append(new)
        nmodes += len(modes)
        a = min(range(len(modes)), key=alts.__getitem__)
        b = min(range(len(modes)), key=lambda i: float(modes[i]["full_unit"]))
        nchanges += a != b
        if not affinity:
            classchanges += (float(modes[a]["rmsd_A"]) <= 2) != (
                float(modes[b]["rmsd_A"]) <= 2
            )
    if maxerr > 0.002:
        raise ValueError("Replacement inconsistent with archived INTER energy")
    return {
        "pools": len(pools),
        "modes": nmodes,
        "max_abs_score_difference": maxerr,
        "selected_rank_differences": nchanges,
        "pose_hit_differences": None if affinity else classchanges,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pose", type=Path)
    p.add_argument("affinity", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    out = {"pose": audit(args.pose, False), "affinity": audit(args.affinity, True)}
    args.output.write_text(json.dumps(out, indent=2) + "\n")
