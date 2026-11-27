"""Replay the shipping binary against frozen ordinary-ring and Kd docking outputs."""

import argparse
import hashlib
import subprocess
import tempfile
from pathlib import Path
from run import HERE, cases, read, write, REFINE, TRAJECTORY


def main(work, binary, out):
    cells = {(p, c): (r, l, b) for p, c, r, l, b in cases(work, "kd")}
    result = []
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "out.pdbqt"
        for r in read(HERE / "kd/runs.csv"):
            if r["arm"] != "candidate":
                continue
            rec, lig, box = cells[r["pdb"], r["compound"]]
            log = subprocess.check_output(
                [
                    binary,
                    "--receptor",
                    rec,
                    "--ligand",
                    lig,
                    "--config",
                    box,
                    "--seed",
                    r["seed"],
                    "--cpu",
                    "1",
                    "--exhaustiveness",
                    "1",
                    "--max_evals",
                    r["cap"],
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--out",
                    path,
                    "--verbosity",
                    "0",
                    "--search_counters",
                ],
                stderr=subprocess.STDOUT,
                text=True,
            )
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            t, e = TRAJECTORY.search(log), REFINE.search(log)
            fields = (
                "hunt_calls",
                "hunt_evals",
                "second_calls",
                "second_evals",
                "explicit_calls",
                "explicit_evals",
            )
            if sha != r["pool_sha256"] or [r[k] for k in fields] != [
                *t.groups(),
                *e.groups(),
            ]:
                raise ValueError(r["key"])
            result.append(dict(key=r["key"], pool_sha256=sha, counters_match=True))
    write(out, result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "binary", "out"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    main(a.work, a.binary, a.out)
