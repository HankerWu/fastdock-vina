"""Replay the candidate on the six historical hunt libraries, or audit its compressed returns."""

import argparse
import csv
import hashlib
import json
import lzma
import os
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bfgs_source_replay"))
from run_replay import ARMS, CASES, digest, run_panel


def rows(content):
    return [line.split("\t") for line in content.decode().splitlines()]


def replays(baseline, output, binary):
    fixtures = baseline / "fixtures/reciprocal"
    totals, receipts = [], []
    for query, receiver, seed in sorted(CASES):
        case = f"{query}-{receiver}-{seed}"
        source = baseline / case
        q, r = [
            json.loads((fixtures / key / "preparation.json").read_text())
            for key in (query, receiver)
        ]
        reference = run_panel.replay.reciprocal.reference(q, r)
        records = (fixtures / query / "ligand.pdbqt").read_text().splitlines()
        serials = [
            int(line[6:11])
            for line in records
            if line.startswith(("ATOM  ", "HETATM"))
            and line[77:].strip() not in {"H", "HD", "HS"}
        ]
        order = [serials.index(i) for i in q["serial"]]

        def errors(table):
            xyz = (
                np.array(
                    [
                        np.fromstring(row[8], dtype=int, sep=" ").reshape(-1, 3)[order]
                        for row in table
                    ]
                )
                / 1000
            )
            distances = np.sum((xyz[:, None] - reference[None]) ** 2, axis=-1).mean(
                axis=-1
            )
            return np.sqrt(distances.min(axis=1))

        for driver in ARMS:
            stem = f"{case}-{driver}"
            states = source / (driver + ".inputs")
            target = output / (stem + ".shadow.xz")
            if binary:
                shadow = output / (stem + ".shadow")
                os.environ.update(FD_BFGS_SHADOW=str(states), FD_BFGS_TRACE=str(shadow))
                cmd = [
                    str(binary),
                    "--ligand",
                    str(fixtures / query / "ligand.pdbqt"),
                    "--receptor",
                    str(fixtures / receiver / "receptor.pdbqt"),
                    "--cpu",
                    "1",
                    "--seed",
                    str(seed),
                    "--max_evals",
                    "12000",
                    "--exhaustiveness",
                    "3",
                    "--num_modes",
                    "9",
                    "--energy_range",
                    "100",
                    "--spacing",
                    "0.375",
                    "--out",
                    str(output / (stem + ".unused")),
                ]
                for axis, value in zip("xyz", r["center"]):
                    cmd += ["--size_" + axis, "24", "--center_" + axis, str(value)]
                stats = run_panel.execute(cmd, output / (stem + ".log"), 180)
                if stats["status"] != "ok":
                    raise ValueError((case, driver, stats))
                target.write_bytes(lzma.compress(shadow.read_bytes(), preset=3))
                receipts.append(
                    dict(
                        case=case,
                        driver=driver,
                        cmd=cmd,
                        states_sha256=digest(states),
                        raw_sha256=digest(shadow),
                        compressed_sha256=digest(target),
                        binary_sha256=digest(binary),
                        **stats,
                    )
                )
                shadow.unlink()
            else:
                old = next(
                    row
                    for row in json.loads((output / "costs.json").read_text())
                    if row["case"] == case and row["driver"] == driver
                )
                assert (
                    digest(states) == old["states_sha256"]
                    and digest(target) == old["compressed_sha256"]
                )
            data = lzma.decompress(target.read_bytes())
            if not binary:
                assert hashlib.sha256(data).hexdigest() == old["raw_sha256"]
            candidate = rows(data)
            rmsd = errors(candidate)
            for name, commit in (("stable", "ed0d42d"), ("common_scalar", "0bb827e")):
                control = rows(
                    (source / (driver + "-" + commit + ".shadow")).read_bytes()
                )
                assert (
                    len(candidate)
                    == len(control)
                    == len(states.read_text().splitlines())
                )
                assert all(
                    a[1] == b[1] == str(i) and a[9] == b[9] and float(a[10]) < 0
                    for i, (a, b) in enumerate(zip(candidate, control))
                )
                dr = rmsd - errors(control)
                de = np.array(
                    [float(a[2]) - float(b[2]) for a, b in zip(candidate, control)]
                )

                def count(records, field):
                    return sum(int(row[field]) for row in records)

                totals.append(
                    dict(
                        case=case,
                        driver=driver,
                        control=name,
                        proposals=len(control),
                        control_calls=count(control, 4),
                        candidate_calls=count(candidate, 4),
                        control_exhausted=count(control, 5),
                        candidate_exhausted=count(candidate, 5),
                        candidate_nondescent=count(candidate, 6),
                        lower_E=int((de < -1e-7).sum()),
                        higher_E=int((de > 1e-7).sum()),
                        lower_RMSD=int((dr < -0.001).sum()),
                        higher_RMSD=int((dr > 0.001).sum()),
                        median_delta_E=float(np.median(de)),
                        median_delta_RMSD_A=float(np.median(dr)),
                        control_best_A=float(errors(control).min()),
                        candidate_best_A=float(rmsd.min()),
                        candidate_within2=int((rmsd <= 2).sum()),
                    )
                )
    if binary:
        (output / "costs.json").write_text(json.dumps(receipts, indent=2) + "\n")
        run_panel.csvfile(output / "locked.csv", totals)
    else:
        with (output / "locked.csv").open() as f:
            expected = list(csv.DictReader(f))
        assert [{k: str(v) for k, v in row.items()} for row in totals] == expected
    print(
        f"Verified {sum(r['proposals'] for r in totals if r['control'] == 'stable')} paired inputs against both controls"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("baseline", type=Path, help="restored bfgs_source_replay archive")
    p.add_argument("output", type=Path)
    p.add_argument("--binary", type=Path)
    a = p.parse_args()
    if a.binary:
        a.output.mkdir(parents=True, exist_ok=False)
    replays(
        a.baseline.resolve(),
        a.output.resolve(),
        a.binary.resolve() if a.binary else None,
    )
