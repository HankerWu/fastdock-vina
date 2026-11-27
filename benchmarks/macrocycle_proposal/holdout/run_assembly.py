"""Run the previously unseen MerTK pocket in its registered two-chain environment."""

import argparse
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "external"))
from panel import command, counts, sha, write
from external_score import main as score


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for field in (
        "work",
        "out",
        "upstream",
        "upstream_counts",
        "parent",
        "reset",
        "ccd",
        "lock",
    ):
        p.add_argument(field, type=Path)
    p.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    d = a.work / "5K0X"
    assembly = json.loads((d / "assembly.json").read_text())
    assert (
        set(assembly["chains"]) == {"A", "B"}
        and assembly["nearest_omitted_nonwater"][0] >= 8
    )
    meta = json.loads((d / "preparation.json").read_text())
    box = d / "box.txt"
    box.write_text(
        "".join(
            f"center_{k} = {v}\nsize_{k} = 26\n" for k, v in zip("xyz", meta["center"])
        )
    )
    case = ("5K0X", "K0X", d / "receptor.pdbqt", d / "ligand.pdbqt", box)
    cells = list(itertools.product((458921, 740317), ((3000, 2000), (12000, 9000))))
    raw = a.out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    rows = {}
    timings = []
    for block in range(4):
        order = cells.copy()
        random.Random(3281 + block // 2).shuffle(order)
        for i, (seed, (cap, lower)) in enumerate(order):
            arms = [
                (arm, cap) for arm in ("parent", "upstream", "ccd", "reset", "lock")
            ] + [("lock", lower)]
            j = (i + block // 2) % len(arms)
            arms = arms[j:] + arms[:j]
            if block % 2:
                arms.reverse()
            for arm, budget in arms:
                key = f"5K0X_{seed}_{budget}_{arm}"
                out = raw / (key + ".pdbqt")
                cmd = command(
                    case, seed, budget, getattr(a, arm), out, arm != "upstream"
                )
                start = time.perf_counter_ns()
                log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
                elapsed = time.perf_counter_ns() - start
                (raw / (key + ".log")).write_text(log)
                if arm == "upstream":
                    with tempfile.TemporaryDirectory() as tmp:
                        replay = Path(tmp) / "counted.pdbqt"
                        log = subprocess.check_output(
                            command(
                                case, seed, budget, a.upstream_counts, replay, False
                            ),
                            text=True,
                            stderr=subprocess.STDOUT,
                        )
                        if out.read_bytes() != replay.read_bytes():
                            raise ValueError(key)
                        (raw / (key + ".counts")).write_text(log)
                values = counts(log)
                digest = sha(out)
                if block:
                    if digest != rows[key]["sha256"] or any(
                        v != rows[key][k] for k, v in values.items()
                    ):
                        raise ValueError(key + " replay differs")
                else:
                    rows[key] = dict(
                        key=key,
                        pdb="5K0X",
                        arm=arm,
                        seed=seed,
                        cap=budget,
                        elapsed_ns=elapsed,
                        **values,
                        joint=sum(v for k, v in values.items() if k.endswith("_evals")),
                        modes=out.read_text().count("ENDMDL"),
                        sha256=digest,
                    )
                timings.append(dict(block=block, key=key, elapsed_ns=elapsed))
    write(a.out / "runs.csv", list(rows.values()))
    write(a.out / "timing.csv", timings)
    score(a.work, a.out)


if __name__ == "__main__":
    main()
