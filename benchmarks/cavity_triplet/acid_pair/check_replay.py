"""Audit chemical state mappings, all released poses, and optional field-query counts."""

import argparse
import csv
import json
from pathlib import Path
import re
import sys
import tarfile

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from run_acid import BENCH, CASES, SEED, prepare, invoke, table, digest

sys.path.insert(0, str(BENCH / "confirmation"))
from run import reference, measure


def clean_check(p, out, field, official_field):
    with tarfile.open(BENCH / "human_panel/prepared.tar.xz") as archive:
        frames = {
            r: json.load(archive.extractfile(r + "/preparation.json"))
            for r in set(CASES) | set(CASES.values())
        }
    states = json.loads((BENCH / "acid_states/states.json").read_text())
    modes = {
        (r["query"], r["state"], r["arm"], int(r["rank"])): r
        for r in csv.DictReader((out / "modes.csv").open())
    }
    maxerr = 0.0
    stats = []
    for q, r in CASES.items():
        coords = []
        for state in ("HA", "A"):
            tag = q + "-" + state
            lns = (p / "inputs" / (tag + ".pdbqt")).read_text().splitlines()
            serial = states[q][state]["serial"]
            raw = {
                int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
                for l in lns
                if l.startswith(("ATOM  ", "HETATM"))
            }
            coords.append(np.array([raw[a] for a in serial]))
            refs = reference(frames[q], frames[r])
            oxy = states[q][state]["oxygen"]
            swapped = refs.copy()
            swapped[:, oxy, :] = refs[:, oxy[::-1], :]
            refs_swap = np.concatenate([refs, swapped])
            for arm in (
                "official",
                "official-double",
                "base",
                "high",
                "matched",
                "sham",
            ):
                file = p / "raw" / (tag + "_" + arm + ".pdbqt")
                fixed = measure(file, serial, refs)
                sym = measure(file, serial, refs_swap)
                for i, (fixed_value, sym_value) in enumerate(
                    zip(fixed, sym, strict=True), 1
                ):
                    row = modes.pop((q, state, arm, i))
                    maxerr = max(
                        maxerr,
                        abs(fixed_value - float(row["fixed_O_A"])),
                        abs(sym_value - float(row["exchange_O_A"])),
                    )
        smi = next(
            l.split()[2]
            for l in (p / "inputs" / (q + "-HA.pdbqt")).read_text().splitlines()
            if l.startswith("REMARK SMILES ")
        )
        scaffold = Chem.MolToSmiles(
            MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smi))
        )
        assert np.array_equal(coords[0], coords[1])
        stats.append(
            dict(
                query=q,
                target=frames[q]["target"],
                receiver=r,
                scaffold=scaffold,
                n_heavy=len(coords[0]),
                max_prepared_heavy_move_A=np.max(abs(coords[0] - coords[1])),
            )
        )
    assert not modes and maxerr < 1e-8
    table(out / "scaffolds.csv", stats)
    (out / "independent.json").write_text(
        json.dumps(
            dict(
                maximum_abs_released_RMSD_delta_A=maxerr,
                unique_scaffolds=len({r["scaffold"] for r in stats}),
            ),
            indent=2,
        )
        + "\n"
    )
    if field:
        assert official_field is not None
        fields = p / "field_check"
        (fields / "raw").mkdir(parents=True, exist_ok=True)
        costs = []
        results = []
        for tag, cell in prepare(p).items():
            for arm in (
                "base",
                "high",
                "matched",
                "sham",
                "official",
                "official-double",
            ):
                stem = tag + "_" + arm
                official = arm.startswith("official")
                env = {}
                if arm in ("matched", "sham"):
                    env = dict(
                        FD_RIGID_INTERVENTION=p / "raw" / (stem + ".input"),
                        FD_RIGID_LOG=fields / "state.log",
                    )
                orig = p / "raw" / (stem + ".pdbqt")
                check = invoke(
                    fields,
                    costs,
                    stem,
                    official_field if official else field,
                    cell,
                    SEED,
                    50000 if arm in ("high", "official-double") else 25000,
                    official,
                    env,
                )
                assert check.read_bytes() == orig.read_bytes()
                log = check.with_suffix(".log").read_text()
                if official:
                    h = dict(
                        (int(d), int(n))
                        for d, n in re.findall(
                            r"GRID_CHANNEL derivative=(\d+) n=(\d+)", log
                        )
                    )
                else:
                    h = dict(
                        (int(d), int(n))
                        for d, n in re.findall(
                            r"FIELD_COUNT combined=0 derivative=(\d+) positions=(\d+)",
                            log,
                        )
                    )
                assert len(h) == 2
                results.append(
                    dict(
                        stem=stem,
                        sha256=digest(check),
                        scalar_grid=h[0],
                        gradient_grid=h[1],
                        wall_s=costs[-1]["wall_s"],
                        cpu_s=costs[-1]["cpu_s"],
                    )
                )
        table(out / "field_calls.csv", results)


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("work", type=Path)
    a.add_argument("results", type=Path)
    a.add_argument("--field-counter", type=Path)
    a.add_argument("--upstream-field", type=Path)
    n = a.parse_args()
    clean_check(
        n.work.resolve(),
        n.results,
        n.field_counter.resolve() if n.field_counter else None,
        n.upstream_field.resolve() if n.upstream_field else None,
    )
