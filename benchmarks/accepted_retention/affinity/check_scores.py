"""Compare two real fixed-coordinate scores for every model in both parent pools."""

import argparse
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path

from run import RECS, SEED, cells, read, write


def check(work, binaries, cpu):
    os.sched_setaffinity(0, {cpu})
    out = []
    with tempfile.TemporaryDirectory() as tmp:
        pose = Path(tmp) / "mode.pdbqt"
        for panel in RECS:
            setups = {(s[1], s[2]): s for s in cells(work, panel)}
            for row in read(work / (panel + "_fixed.csv")):
                if row["arm"] != "parent":
                    continue
                _, rec, c, protein, box, lig = setups[row["receptor"], row["compound"]]
                models = (
                    (work / "raw" / (row["key"] + ".pdbqt")).read_text().split("ENDMDL")
                )
                if len(models) != 1 + int(row["modes"]):
                    raise ValueError(f"Incomplete fixed parent pool {row['key']}")
                for i, model in enumerate(models[:-1], 1):
                    pose.write_text(
                        "\n".join(
                            l
                            for l in model.splitlines()
                            if l and not l.startswith("MODEL")
                        )
                        + "\n"
                    )
                    common = [
                        "--score_only",
                        "--receptor",
                        protein,
                        "--ligand",
                        lig,
                        "--config",
                        box,
                        "--unbound_energy",
                        "0",
                        "--seed",
                        str(SEED),
                        "--cpu",
                        "1",
                        "--pose_coordinates",
                        pose,
                    ]
                    parts = []
                    for binary in binaries:
                        proc = subprocess.run(
                            [binary, *common],
                            capture_output=True,
                            text=True,
                            check=True,
                        )
                        components = [
                            l
                            for l in proc.stdout.splitlines()
                            if "(kcal/mol)" in l or "Ligand - " in l or "Flex   - " in l
                        ]
                        if len(components) != 10:
                            raise ValueError((rec, c, i, components))
                        parts.append(components)
                    if parts[0] != parts[1]:
                        raise ValueError(("different native scores", rec, c, i, parts))
                    out.append(
                        {
                            "panel": panel,
                            "receptor": rec,
                            "compound": c,
                            "rank": i,
                            "equal_components": len(parts[0]),
                            "pose_sha256": hashlib.sha256(
                                pose.read_bytes()
                            ).hexdigest(),
                            "decomposition_sha256": hashlib.sha256(
                                "\n".join(parts[0]).encode()
                            ).hexdigest(),
                            "score_line": parts[0][0].strip(),
                        }
                    )
    write(work / "fixed_scores.csv", out)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("parent", type=Path)
    p.add_argument("candidate", type=Path)
    p.add_argument("--cpu", type=int, default=6)
    args = p.parse_args()
    check(args.work, [args.parent, args.candidate], args.cpu)
