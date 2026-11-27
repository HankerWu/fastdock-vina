"""Interleave five engines, record child CPU and peak RSS, and verify replay consistency."""

import argparse
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
MACRO = HERE.parent / "macrocycle_proposal"
sys.path.insert(0, str(MACRO))
sys.path.insert(0, str(MACRO / "external"))
from panel import command, counts, sha, write
from external_score import main as evaluate


def timed(cmd, path):
    start = time.perf_counter_ns()
    with path.open("wb") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    elapsed = time.perf_counter_ns() - start
    if proc.returncode:
        raise RuntimeError(path.read_text())
    return elapsed, usage.ru_utime + usage.ru_stime, usage.ru_maxrss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in [
        "inputs",
        "output",
        "upstream",
        "upstream_counts",
        "parent",
        "lock",
        "null",
        "concerted",
    ]:
        parser.add_argument(name, type=Path)
    parser.add_argument("--cpu", type=int, default=min(os.sched_getaffinity(0)))
    parser.add_argument(
        "--lower",
        action="store_true",
        help="Use 9000 for concerted and 12000 for parent and upstream",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Original and Pro at 12000 vs three restricted engines at 9000 on fresh seeds",
    )
    parser.add_argument("--ids", nargs="+", default=["5T9W", "7MBO", "5K0X"])
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    raw = args.output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    engines = {
        n: getattr(args, n).resolve()
        for n in ("upstream", "parent", "lock", "null", "concerted")
    }
    if args.lower:
        engines = {
            k: v for k, v in engines.items() if k in ("parent", "upstream", "concerted")
        }
    seeds = (421831, 692441, 886117) if args.fresh else (191009, 458921, 740317)
    cells = list(
        itertools.product(
            args.ids, seeds, (12000,) if args.lower or args.fresh else (3000, 12000)
        )
    )
    random.Random(9049).shuffle(cells)
    timing, rows = [], {}
    for rep in range(args.repeats):
        for index, (pdb, seed, nominal) in enumerate(cells):
            directory = args.inputs / pdb
            center = json.loads((directory / "preparation.json").read_text())["center"]
            box = args.output / (pdb + ".box")
            box.write_text(
                "".join(
                    f"center_{d} = {v}\nsize_{d} = 26\n" for d, v in zip("xyz", center)
                )
            )
            case = (
                pdb,
                "",
                directory / "receptor.pdbqt",
                directory / "ligand.pdbqt",
                box,
            )
            order = list(engines)
            offset = index % len(order)
            order = order[offset:] + order[:offset]
            if rep % 2:
                order.reverse()
            for name in order:
                reduced = (args.lower and name == "concerted") or (
                    args.fresh and name in ("lock", "null", "concerted")
                )
                cap = 9000 if reduced else nominal
                key = f"{pdb}_{seed}_{cap}_{name}"
                dest = raw / (key + ".pdbqt")
                log = raw / (key + ".log")
                elapsed, cpu, rss = timed(
                    command(case, seed, cap, engines[name], dest, name != "upstream"),
                    log,
                )
                if name == "upstream":
                    replay, counted = raw / (key + ".tmp"), raw / (key + ".counts")
                    if not rep:
                        timed(
                            command(
                                case, seed, cap, args.upstream_counts, replay, False
                            ),
                            counted,
                        )
                        assert replay.read_bytes() == dest.read_bytes()
                        replay.unlink()
                    log = counted
                work = counts(log.read_text())
                digest = sha(dest)
                if rep:
                    assert rows[key]["sha256"] == digest and all(
                        rows[key][k] == v for k, v in work.items()
                    )
                else:
                    rows[key] = dict(
                        key=key,
                        pdb=pdb,
                        arm=name,
                        seed=seed,
                        cap=cap,
                        **work,
                        joint=sum(v for k, v in work.items() if k.endswith("_evals")),
                        modes=dest.read_text().count("ENDMDL"),
                        sha256=digest,
                    )
                timing.append(
                    dict(
                        key=key, repeat=rep, wall_ns=elapsed, cpu_s=cpu, max_rss_KiB=rss
                    )
                )
    write(args.output / "runs.csv", list(rows.values()))
    write(args.output / "timings.csv", timing)
    evaluate(args.inputs, args.output)
    (args.output / "binaries.json").write_text(
        json.dumps({k: sha(v) for k, v in engines.items()}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
