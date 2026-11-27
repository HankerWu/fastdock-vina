"""Replay broken and intact rapamycin poses plus the crystallographic BACE1 start."""

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from run import cases, read, write
from geometry import coord_blocks


def main(a):
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb, cid, protein, original, box in cases(a.work, "macro"):
        if pdb not in ("1FKB", "5QCT"):
            continue
        initial = next(coord_blocks(original))
        ends = [i for i, (ty, v) in initial.items() if ty == "CG0"]

        def bond_length(pose):
            return float(np.linalg.norm(pose[ends[0]][1] - pose[ends[1]][1]))

        def table(path):
            return subprocess.check_output([a.probe, "table", "XS", "1", path])

        expected = table(original)
        if pdb == "5QCT":
            starts = [("5QCT_crystal", original)]
        else:
            meta = json.loads((a.work / "extra/1FKB/preparation.json").read_text())
            starts = []
            for row in read(a.tail / "runs.csv"):
                if row["pdb"] != pdb or row["arm"] != "parent" or row["cap"] != "3000":
                    continue
                text = (a.tail / "raw" / (row["key"] + ".pdbqt")).read_text()
                for rank, block in enumerate(
                    re.findall(r"^MODEL [^\n]*\n(.*?)^ENDMDL", text, re.M | re.S), 1
                ):
                    key = row["key"] + f"_{rank}"
                    path = raw / (key + ".input.pdbqt")
                    path.write_text(block)
                    starts.append((key, path))
        for key, path in starts:
            actual = table(path)
            if expected != actual:
                raise ValueError("Reparsing changed pair topology: " + key)
            for arm, binary in (("parent", a.parent), ("tail", a.candidate)):
                prefix = raw / (key + "_" + arm)
                cmd = [
                    binary,
                    "--local_only",
                    "--refinement_optimizer",
                    "armijo",
                    "--ligand",
                    path,
                    "--receptor",
                    protein,
                    "--config",
                    box,
                    "--cpu",
                    "1",
                    "--seed",
                    "740317",
                    "--out",
                    prefix.with_suffix(".pdbqt"),
                ]
                log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
                prefix.with_suffix(".log").write_text(log)
                before = next(coord_blocks(path))
                after = next(coord_blocks(prefix.with_suffix(".pdbqt")))
                keyfields = dict(
                    key=key,
                    arm=arm,
                    bond_before=bond_length(before),
                    bond_after=bond_length(after),
                    gradient_evals=int(re.search(r"gradient_evals=(\d+)", log)[1]),
                    pair_table_sha256=hashlib.sha256(actual).hexdigest(),
                )
                if pdb == "1FKB":

                    def rmsd(pose):
                        x = np.array([pose[i][1] for i in meta["serial"]])
                        return float(
                            np.sqrt(
                                ((x - np.asarray(meta["mapped"])) ** 2)
                                .sum(axis=2)
                                .mean(axis=1)
                            ).min()
                        )

                    keyfields.update(rmsd_before=rmsd(before), rmsd_after=rmsd(after))
                else:
                    # The prepared BACE input agrees with the registered crystal to 0.00051 A.
                    ix = [
                        i
                        for i, (ty, v) in before.items()
                        if ty not in ("HD", "H", "G0")
                    ]
                    keyfields.update(
                        rmsd_before=0.0,
                        rmsd_after=float(
                            np.sqrt(
                                np.mean(
                                    [
                                        np.sum((after[i][1] - before[i][1]) ** 2)
                                        for i in ix
                                    ]
                                )
                            )
                        ),
                    )
                rows.append(keyfields)
    write(a.out / "local.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "tail", "out", "parent", "candidate", "probe"):
        p.add_argument(name, type=Path)
    main(p.parse_args())
