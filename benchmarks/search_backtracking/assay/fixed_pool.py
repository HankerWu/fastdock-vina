"""Score frozen parent poses and measure close heavy-atom contacts in every first pose."""

import argparse
import hashlib
import re
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np
from analyze import read, write
from scipy.spatial import cKDTree


def xyz(lines):
    return np.array(
        [
            [float(s[k : k + 8]) for k in (30, 38, 46)]
            for s in lines
            if s.startswith(("ATOM  ", "HETATM"))
            and s.split()[-1] not in ("H", "HD", "HS")
        ]
    )


def main(work, bins, output):
    pre = work / "prepared"
    trees = {
        p: cKDTree(xyz((pre / (p + ".pdbqt")).read_text().splitlines()))
        for p in ("1A9U", "1M17")
    }
    scores, contacts = [], []
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "pose.pdbqt"
        for study in ("experiments", "cost_check"):
            for row in read(work / study / "runs.csv"):
                rec, cid, key = row["receptor"], row["compound"], row["key"]
                models = (
                    (work / study / "raw" / (key + ".pdbqt"))
                    .read_text()
                    .split("ENDMDL")[:-1]
                )
                assert len(models) == int(row["modes"]) and models
                near = trees[rec].query(xyz(models[0].splitlines()))[0]
                contacts.append(
                    {
                        "study": study,
                        "key": key,
                        "min_heavy_A": round(float(min(near)), 5),
                        "atoms_within_2A": int((near < 2).sum()),
                        "first_score": row["first_score"],
                    }
                )
                if (
                    study != "experiments"
                    or row["cap"] != "12000"
                    or row["arm"] != "parent"
                ):
                    continue
                for rank, model in enumerate(models, 1):
                    path.write_text(
                        "\n".join(
                            s
                            for s in model.splitlines()
                            if s and not s.startswith("MODEL")
                        )
                        + "\n"
                    )
                    args = [
                        "--score_only",
                        "--receptor",
                        pre / (rec + ".pdbqt"),
                        "--ligand",
                        pre / (cid + ".pdbqt"),
                        "--config",
                        pre / (rec + ".box"),
                        "--unbound_energy",
                        "0",
                        "--pose_coordinates",
                        path,
                        "--cpu",
                        "1",
                        "--seed",
                        row["seed"],
                    ]
                    parts = []
                    start = time.perf_counter()
                    for binary in bins:
                        log = subprocess.check_output([binary, *args], text=True)
                        parts.append(
                            [
                                s
                                for s in log.splitlines()
                                if "(kcal/mol)" in s
                                or "Ligand - " in s
                                or "Flex   - " in s
                            ]
                        )
                    if len(parts[0]) != 10 or any(v != parts[0] for v in parts[1:]):
                        raise ValueError((key, rank, parts))
                    vals = [float(re.search(r":\s+([+.\d-]+)", s)[1]) for s in parts[0]]
                    scores.append(
                        {
                            "key": key,
                            "rank": rank,
                            "arms": len(bins),
                            "pose_sha256": hashlib.sha256(
                                path.read_bytes()
                            ).hexdigest(),
                            "terms": " ".join(map(str, vals)),
                            "diagnostic_seconds": time.perf_counter() - start,
                        }
                    )
    write(output / "fixed_scores.csv", scores)
    write(output / "clashes.csv", contacts)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("binaries", type=Path, nargs=3)
    a = p.parse_args()
    main(a.work, a.binaries, a.output)
