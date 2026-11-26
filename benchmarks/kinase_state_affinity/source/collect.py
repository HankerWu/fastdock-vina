"""Collect saved top modes and a target-independent prepared-ligand reference."""

import argparse
import csv
import hashlib
import itertools
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

from freeze import RECEPTORS, SEEDS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "implementation"))
from ensemble_affinity import pair_sensitivity


def write(path, rows):
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def score(line, tag):
    return float(
        next(l for l in line.splitlines() if l.startswith(tag))
        .split(":")[1]
        .strip()
        .split()[0]
    )


def collect(work, binary, arm="native", seeds=SEEDS):
    results = work / ("results" if arm == "native" else f"{arm}-results")
    results.mkdir(exist_ok=True)
    refs = work / "references"
    refs.mkdir(exist_ok=True)
    with (work / "panel.csv").open() as f:
        panel = list(csv.DictReader(f))
    reference = {}
    for cid in sorted({r["cid"] for r in panel}):
        runs = []
        for e0 in (0, 100):
            logfile = refs / f"{cid}-{e0}.log"
            with logfile.open("w") as log:
                subprocess.run(
                    [
                        binary,
                        "--score_only",
                        "--autobox",
                        "--unbound_energy",
                        str(e0),
                        "--cpu",
                        "1",
                        "--receptor",
                        work / "prepared/1HCK.pdbqt",
                        "--ligand",
                        work / f"prepared/{cid}.pdbqt",
                    ],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            runs.append(logfile.read_text())
        intra = score(runs[0], "(2) Final Total Internal Energy")
        slope = (
            score(runs[0], "Estimated Free Energy of Binding")
            - score(runs[1], "Estimated Free Energy of Binding")
        ) / 100
        assert 0 < slope <= 1 and intra == score(
            runs[1], "(2) Final Total Internal Energy"
        )
        reference[cid] = (intra, slope)
    cells = []
    samples = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for r in panel:
        target = r["target"]
        cid = r["cid"]
        intra_ref, slope = reference[cid]
        for state, _ in RECEPTORS[target]:
            for seed in seeds:
                path = (
                    work
                    / ("pools" if arm == "native" else f"{arm}-pools")
                    / f"{state}_{cid}_{seed}.pdbqt"
                )
                top = path.read_text().split("ENDMDL")[0]
                native = score(top, "REMARK VINA RESULT")
                unbound = score(top, "REMARK UNBOUND")
                total = score(top, "REMARK INTER + INTRA")
                assert abs(native - slope * (total - unbound)) < 0.004
                fixed = slope * (total - intra_ref)
                samples[target][cid][state].append((native, fixed))
                cells.append(
                    dict(
                        target=target,
                        ligand=cid,
                        state=state,
                        replicate=seed,
                        native=native,
                        total=total,
                        pool_unbound=unbound,
                        prepared_unbound=intra_ref,
                        normalization=slope,
                        energy=fixed,
                        pool_sha256=digest(path),
                    )
                )
    write(results / "cells.csv", cells)
    pairs = []
    for target, drugs in samples.items():
        lookup = {r["cid"]: r for r in panel if r["target"] == target}
        states = [state for state, _ in RECEPTORS[target]]
        for a, b in itertools.combinations(drugs, 2):
            ka, kb = (float(lookup[c]["kd_nm"]) for c in (a, b))
            if ka > kb:
                a, b = b, a
                ka, kb = kb, ka
            if kb < 3 * ka:
                continue
            record = dict(
                target=target,
                strong=lookup[a]["compound"],
                weak=lookup[b]["compound"],
                strong_cid=a,
                weak_cid=b,
                kd_strong_nm=ka,
                kd_weak_nm=kb,
                state_0=states[0],
                state_1=states[1],
            )
            for i, arm in enumerate(("native", "fixed")):
                va, vb = (
                    [median(obs[i] for obs in drugs[c][s]) for s in states]
                    for c in (a, b)
                )
                g = pair_sensitivity(va, vb)
                record.update({f"{arm}_{k}": v for k, v in g.items()})
                record.update(
                    {
                        f"{arm}_state0": vb[0] - va[0],
                        f"{arm}_state1": vb[1] - va[1],
                        f"{arm}_mean": sum(vb) / 2 - sum(va) / 2,
                    }
                )
            pairs.append(record)
    write(results / "pairs.csv", pairs)
    return cells, pairs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--arm", choices=("native", "mixed", "soft"), default="native")
    parser.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    args = parser.parse_args()
    collect(args.work, args.binary.resolve(), args.arm, args.seeds)
