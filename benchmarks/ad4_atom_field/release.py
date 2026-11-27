"""Rerun every joint search after the release-domain check, counting empty pools as failures."""

import argparse
import csv
import os
from pathlib import Path
import re
import subprocess
import time
from run_panel import write


def run(work, root, binary, key, ntraj):
    q, r, seed, cap = key
    stem = root / f"{q}-{r}-{seed}-{cap}-{ntraj}"
    pdbqt = stem.with_suffix(".pdbqt")
    # The caller may reuse an output path. A failed search must not leave stale models.
    pdbqt.write_text("OLD RESULT\n")
    with stem.with_suffix(".log").open("w") as out:
        t = time.perf_counter()
        p = subprocess.Popen(
            [
                binary,
                "--ligand",
                work / q / "ligand.pdbqt",
                "--scoring",
                "ad4",
                "--maps",
                work / "maps" / r / "receptor",
                "--seed",
                seed,
                "--max_evals",
                cap,
                "--num_modes",
                "5",
                "--energy_range",
                "100",
                "--cpu",
                "1",
                "--search_counters",
                "--exhaustiveness",
                str(ntraj),
                "--out",
                pdbqt,
            ],
            stdout=out,
            stderr=subprocess.STDOUT,
        )
        _, status, usage = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
        elapsed = time.perf_counter() - t
    if p.returncode:
        raise RuntimeError(stem)
    txt = pdbqt.read_text()
    assert "OLD RESULT" not in txt
    log = stem.with_suffix(".log").read_text()
    counts = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", log)
    assert len(counts) == ntraj
    return (
        [x.split("\n", 1)[1] for x in txt.split("ENDMDL\n") if x.strip()],
        sum(int(t) for row in counts for t in row),
        elapsed,
        usage.ru_utime + usage.ru_stime,
        usage.ru_maxrss,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "panel", "released", "pre", "post"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    a.released.mkdir(exist_ok=True, parents=True)
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    rows = []
    data = [
        r
        for f in sorted(a.panel.glob("*_cells.csv"))
        for r in csv.DictReader(f.open())
        if r["arm"] == "joint"
    ]
    for row in data:
        key = tuple(row[n] for n in ("query", "receiver", "seed", "cap"))
        prev = (a.panel / "runs" / ("-".join(key) + "-joint.pdbqt")).read_text()
        prev = [x.split("\n", 1)[1] for x in prev.split("ENDMDL\n") if x.strip()]
        blocks, evals, elapsed, cput, rss = run(a.work, a.released, a.post, key, 1)
        kept = [prev.index(x) + 1 for x in blocks if x in prev]
        rows.append(
            dict(
                query=key[0],
                receiver=key[1],
                seed=key[2],
                cap=key[3],
                before=len(prev),
                after=len(blocks),
                preserved_ranks=",".join(map(str, kept)),
                identical_coordinates_scores=int(len(kept) == len(blocks)),
                top1_hit=row["top1"] if blocks and kept and kept[0] == 1 else 0,
                evals=evals,
                unchanged_evals=int(evals == int(row["evals"])),
                wall_s=elapsed,
                cpu_s=cput,
                peak_rss_KiB=rss,
            )
        )
    write(a.released / "release.csv", rows)
    # A known infeasible trajectory can coexist with a valid one. Replay both
    # maps across independent trajectories, retaining every eligible survivor.
    multi = []
    for receiver in ("1FKG", "1FKH"):
        key = (receiver, receiver, "30553", "3000")
        for arm, binary in (("before", a.pre), ("after", a.post)):
            root = a.released / arm
            root.mkdir(exist_ok=True)
            blocks, ev, elapsed, cput, rss = run(a.work, root, binary, key, 3)
            multi.append(
                dict(
                    receiver=receiver,
                    arm=arm,
                    models=len(blocks),
                    evals=ev,
                    scores=",".join(
                        re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", x)[1]
                        for x in blocks
                    ),
                )
            )
    write(a.released / "multiple.csv", multi)


if __name__ == "__main__":
    main()
