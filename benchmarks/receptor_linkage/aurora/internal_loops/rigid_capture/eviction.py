"""Locate loss of the near-crystal t1/s1 witness inside the original online container."""

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as audit


def main(work, binary, original_raw, out):
    out.mkdir(parents=True, exist_ok=False)
    spec = original_raw / "fit_t1_s1.input"
    pose, events, xyz = out / "pool.pdbqt", out / "events.txt", out / "stages.xyz"
    box = [
        a
        for d, v in zip(
            "xyz",
            json.loads((work / "prepared/2X6E.json").read_text())["center"],
            strict=True,
        )
        for a in (f"--center_{d}", str(v), f"--size_{d}", "26")
    ]
    args = [
        binary,
        "--receptor",
        work / "template1OL7/trans_model/released/receptor.pdbqt",
        "--ligand",
        work / "prepared/at9283.pdbqt",
        "--ligand_bonds",
        work / "raw/at9283_chemical.bonds",
        *box,
        "--cpu",
        "1",
        "--exhaustiveness",
        "2",
        "--seed",
        "61231",
        "--num_modes",
        "9",
        "--energy_range",
        "100",
        "--max_evals",
        "25000",
        "--pose_selection",
        "late",
        "--search_counters",
        "--out",
        pose,
    ]
    tic = time.perf_counter()
    with (out / "search.log").open("w") as f:
        subprocess.run(
            args,
            env=os.environ
            | {
                k: str(v)
                for k, v in {
                    "FD_RIGID_HISTORY": events,
                    "FD_RIGID_COORDS": xyz,
                    "FD_RIGID_LOG": out / "injection.log",
                    "FD_RIGID_INTERVENTION": spec,
                }.items()
            },
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True,
        )
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    timing = dict(wall_s=time.perf_counter() - tic, cpu_s=cpu.ru_stime + cpu.ru_utime)
    assert pose.read_bytes() == (original_raw / "fit_t1_s1.pdbqt").read_bytes()
    query, order = audit.frame.read_query(work, "at9283")
    crystal, _ = audit.refmol(work, "at9283")
    maps = crystal.GetSubstructMatches(query, uniquify=False, maxMatches=10000)
    assert 0 < len(maps) < 10000
    with xyz.open() as f:
        codes = list(map(int, next(f)[1:].split()))
        ref = crystal.GetConformer().GetPositions()[np.array(maps)][
            :, [order.index(s) for s in codes]
        ]
        offers = {}
        for line in f:
            tag, tr, k, *vals = line.split()
            if tag == "offered":
                assert tr == "1"
                offers[int(k)] = np.array(vals, float).reshape(-1, 3)
    rms = lambda x, y: float(
        np.sqrt(np.mean(np.sum((x - y) ** 2, axis=-1), axis=-1)).min()
    )
    errors = {k: rms(x, ref) for k, x in offers.items()}
    previous = []
    energies = {}
    drops = []
    reported = []
    for ln in events.read_text().splitlines():
        i, energy, capacity, threshold, *ids = ln.split()
        i, energy, capacity, threshold = (
            int(i),
            float(energy),
            int(capacity),
            float(threshold),
        )
        resident = list(map(int, ids))
        energies[i] = energy
        geom = offers[i]
        # Independent replay of sorted direct-index RMSD exclusion and capacity.
        copy = list(previous)
        if len(copy) < capacity or energy < energies[copy[-1]]:
            first = next(
                (j for j, k in enumerate(copy) if energy < energies[k]), len(copy)
            )
            if all(rms(geom, offers[k]) >= threshold for k in copy[:first]):
                copy = (
                    copy[:first]
                    + [i]
                    + [k for k in copy[first:] if rms(geom, offers[k]) >= threshold]
                )[:capacity]
        # A repeated exact conformation can carry an older native identity.
        assert len(copy) == len(resident)
        for a, b in zip(copy, resident, strict=True):
            assert energies[a] == energies[b] and np.array_equal(offers[a], offers[b])
        prev_hits = [k for k in previous if errors[k] <= 2]
        hits = [k for k in resident if errors[k] <= 2]
        if prev_hits and not hits:
            drops.append(
                dict(
                    step=i,
                    previous=previous,
                    previous_near_ids=prev_hits,
                    current=resident,
                    capacity=capacity,
                    exclusion_A=threshold,
                    arrival_A=errors[i],
                    arrival_score=energy,
                    lost=[
                        dict(
                            step=k,
                            rmsd_A=errors[k],
                            score=energies[k],
                            distance_to_arrival_no_sym_A=rms(offers[k], geom),
                        )
                        for k in prev_hits
                    ],
                    smallest_remaining_rmsd_A=min(errors[k] for k in resident),
                )
            )
        if hits:
            reported.append(i)
        previous = resident
    import re

    n = sum(
        map(
            int,
            re.findall(
                r"\b(?:hunt|second|explicit)_evals=(\d+)",
                (out / "search.log").read_text(),
            ),
        )
    )
    doc = dict(
        offers_verified=len(energies),
        source_states=len(offers),
        loss_transitions=drops,
        near_container_first_last_offer_step=[min(reported), max(reported)],
        final_min_A=min(errors[k] for k in previous),
        normal_evals=n,
        additional_objective_evals=0,
        binary_sha256=audit.digest(binary),
        search=timing,
    )
    (out / "eviction.json").write_text(json.dumps(doc, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "binary", "original_raw", "out"):
        p.add_argument(key, type=Path)
    main(*(v.resolve() for v in vars(p.parse_args()).values()))
