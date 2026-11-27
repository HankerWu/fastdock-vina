"""Exploratory local relaxation of the entire 61231 AT9283 pool pair on the audited graph."""

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay
from pose_probe import local_refinement, score_pose
from summarize import read, save


def run(work, blind, pro, dest, fixed):
    dest.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    drug, seed = "at9283", "61231"
    ligand, bonds, serials, refs = replay.parameters(work, drug)
    centre = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    cmd = [pro, "--receptor", blind / "both.pdbqt", "--cpu", "1", "--seed", seed]
    for axis, value in zip("xyz", centre, strict=True):
        cmd.extend(["--center_" + axis, str(value), "--size_" + axis, "26"])
    score = [*cmd, "--ligand", ligand, "--ligand_bonds", bonds]
    records = {
        (r["source"], int(r["rank"])): r
        for r in read(fixed)
        if r["drug"] == drug and r["receiver"] == "both" and r["seed"] == seed
    }
    out = []
    for source in ("both", "both_graph"):
        pool = blind / "historical/blind" / f"{source}_{drug}_{seed}.pdbqt"
        blocks = [b for b in pool.read_text().split("ENDMDL") if b.strip()]
        for rank, block in enumerate(blocks, 1):
            initial = dest / f"{source}_{rank}.pdbqt"
            initial.write_text(
                "\n".join(
                    l for l in block.strip().splitlines() if not l.startswith("MODEL")
                )
                + "\n"
            )
            raw = score_pose(score, initial)
            if abs(raw - float(records[source, rank]["complete_score"])) > 0.00051:
                raise ValueError(
                    "Independent CLI fixed-coordinate score disagrees with probe"
                )
            x = replay.geom.mode_coords(initial, serials)[0]
            for opt, options in [
                ("native", ()),
                ("armijo", ("--refinement_optimizer", "armijo", "--max_evals", "250")),
            ]:
                final = dest / f"{source}_{rank}_{opt}.pdbqt"
                local_refinement(
                    cmd, initial, bonds, final, final.with_suffix(".log"), options
                )
                s = score_pose(score, final)
                y = replay.geom.mode_coords(final, serials)[0]
                out.append(
                    dict(
                        source=source,
                        rank=rank,
                        optimizer=opt,
                        initial_A=replay.mapped_rmsd(x, refs),
                        relaxed_A=replay.mapped_rmsd(y, refs),
                        moved_same_atom_rms_A=np.sqrt(
                            np.mean(np.sum((y - x) ** 2, axis=1))
                        ),
                        input_score=raw,
                        output_score=s,
                    )
                )
    save(dest / "relax.csv", out)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    print(
        json.dumps(
            dict(
                n_starts=len(records),
                n_relaxations=len(out),
                binary_sha256=replay.digest(pro),
                wall_s=time.perf_counter() - start,
                cli_cpu_s=after.ru_utime
                + after.ru_stime
                - before.ru_utime
                - before.ru_stime,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("work", "blind", "pro", "out", "fixed_scores"):
        p.add_argument(k, type=Path)
    a = p.parse_args()
    run(*(v.resolve() for v in (a.work, a.blind, a.pro, a.out, a.fixed_scores)))
