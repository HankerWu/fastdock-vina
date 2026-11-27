"""Cross-score fixed coordinates under the two graphs, on their generating receptor."""

import argparse
import csv
import hashlib
import json
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay
from summarize import read, save


def compare(work, run, probe, dest):
    start = time.monotonic()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    rows, pools = [], []
    for drug in replay.geom.REFERENCES:
        ligand, bonds, _, _ = replay.parameters(work, drug)
        for receiver in ("far_both", "both"):
            with tempfile.TemporaryDirectory(dir=dest) as directory:
                paths, poses = [], []
                for cohort in ("historical", "fresh"):
                    folder = run / cohort
                    modes = {
                        (r["arm"], r["drug"], r["seed"], int(r["rank"])): r
                        for r in read(folder / "modes.csv")
                    }
                    for cell in read(folder / "cells.csv"):
                        source, seed = cell["arm"], cell["seed"]
                        if drug != cell["drug"] or source not in (
                            receiver,
                            receiver + "_graph",
                        ):
                            continue
                        pool = folder / "blind" / f"{source}_{drug}_{seed}.pdbqt"
                        if replay.digest(pool) != cell["pool_sha256"]:
                            raise ValueError(f"Unexpected pool hash: {pool}")
                        blocks = [
                            b for b in pool.read_text().split("ENDMDL") if b.strip()
                        ]
                        if len(blocks) != int(cell["n"]):
                            raise ValueError("Wrong mode count")
                        for rank, block in enumerate(blocks, 1):
                            text = (
                                "\n".join(
                                    l
                                    for l in block.strip().splitlines()
                                    if not l.startswith("MODEL")
                                )
                                + "\n"
                            )
                            path = Path(directory) / f"{len(paths)}.pdbqt"
                            path.write_text(text)
                            paths.append(str(path.resolve()))
                            xyz = "\n".join(
                                l[6:11] + l[30:54]
                                for l in text.splitlines()
                                if l.startswith(("HETATM", "ATOM  "))
                            )
                            poses.append(
                                dict(
                                    cohort=cohort,
                                    drug=drug,
                                    receiver=receiver,
                                    seed=seed,
                                    source=source,
                                    rank=rank,
                                    coord_sha=hashlib.sha256(xyz.encode()).hexdigest(),
                                    rmsd_A=modes[source, drug, seed, rank]["rmsd_A"],
                                )
                            )
                for graph, args in (("inferred", []), ("complete", [bonds])):
                    result = subprocess.run(
                        [
                            probe,
                            run / (receiver + ".pdbqt"),
                            ligand,
                            *map(str, center),
                            *args,
                        ],
                        input="\n".join(paths) + "\n",
                        capture_output=True,
                        text=True,
                        check=True,
                    )
                    for pose, energies in zip(
                        poses, csv.reader(result.stdout.splitlines()), strict=True
                    ):
                        for name, value in zip(
                            ("score", "inter", "intra"), energies, strict=True
                        ):
                            pose[graph + "_" + name] = float(value)
                rows.extend(poses)
                for seed in sorted({r["seed"] for r in poses}):
                    for source in (receiver, receiver + "_graph", "union"):
                        group = [
                            r
                            for r in poses
                            if r["seed"] == seed and source in (r["source"], "union")
                        ]
                        record = dict(
                            cohort=group[0]["cohort"],
                            drug=drug,
                            receiver=receiver,
                            seed=seed,
                            source=source,
                            n=len(group),
                        )
                        winners = {}
                        for graph in ("inferred", "complete"):
                            m = min(r[graph + "_score"] for r in group)
                            tied = [r for r in group if r[graph + "_score"] == m]
                            winners[graph] = {r["coord_sha"] for r in tied}
                            record[graph + "_winner_ids"] = "|".join(
                                r["source"] + "#" + str(r["rank"]) for r in tied
                            )
                            record[graph + "_winner_min_A"] = min(
                                float(r["rmsd_A"]) for r in tied
                            )
                            record[graph + "_winner_max_A"] = max(
                                float(r["rmsd_A"]) for r in tied
                            )
                            record[graph + "_min_score"] = m
                        pools.append(
                            dict(
                                **record,
                                same_winners=int(
                                    winners["inferred"] == winners["complete"]
                                ),
                            )
                        )
    save(dest / "fixed_graph_scores.csv", rows)
    save(dest / "fixed_graph_pools.csv", pools)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(
        wall_s=time.monotonic() - start,
        child_cpu_s=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
        n_poses=len(rows),
        scores=2 * len(rows),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "run", "probe", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    print(
        json.dumps(
            compare(*(x.resolve() for x in (a.work, a.run, a.probe, a.out))), indent=2
        )
    )
