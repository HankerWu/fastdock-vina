"""Repeat 60 actual single points with the original Fortran xTB program and audit pairing."""

import argparse
import math
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import median

from rdkit import Chem
from data import RT, read, write
from quantum import KCAL, binding_site, ligand, merge, read_models, rotations


def legacy(task, exe):
    key, geom, solvent = task
    table = Chem.GetPeriodicTable()
    with tempfile.TemporaryDirectory() as tmp:
        xyz = Path(tmp) / "input.xyz"
        xyz.write_text(
            f"{len(geom.xyz)}\nfixed coordinates in angstrom\n"
            + "".join(
                f"{table.GetElementSymbol(int(z))} {x:.12f} {y:.12f} {w:.12f}\n"
                for z, (x, y, w) in zip(geom.numbers, geom.xyz, strict=True)
            )
        )
        flags = ["--alpb", solvent] if solvent else []
        env = {
            **os.environ,
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
        log = subprocess.run(
            [
                exe,
                xyz,
                "--gfn",
                "2",
                "--chrg",
                str(geom.charge),
                "--acc",
                "0.05",
                "--norestart",
                *flags,
            ],
            cwd=tmp,
            env=env,
            check=True,
            text=True,
            capture_output=True,
        ).stdout
        match = re.search(r"\| TOTAL ENERGY\s+([-0-9.]+) Eh", log)
        if not match:
            raise ValueError(key, log[-1024:])
        return key, float(match[1]) * KCAL


def original_jobs(folder):
    for rec in ("1BTY", "3PTB"):
        ps = rotations(
            binding_site(folder / (rec + ".pdbqt")), ("triad", "delta", "water")
        )
        for i, pocket in enumerate(ps):
            yield (rec, "", "pocket", i), pocket, "water"
        for name in ("H", "NH2", "OCH3"):
            path = folder / f"{rec}_{name}_271828_shell.pdbqt"
            l = ligand(read_models(path.read_text().splitlines())[1])
            yield (rec, name, "ligand_vacuum", 0), l, None
            yield (rec, name, "ligand_water", 0), l, "water"
            for i, pocket in enumerate(ps):
                yield (rec, name, "complex", i), merge(pocket, l), "water"


def audit_legacy(prepared, results, exe):
    with ThreadPoolExecutor(3) as workers:
        energies = dict(workers.map(lambda j: legacy(j, exe), original_jobs(prepared)))
    reference = {}
    for row in read(results / "energies.csv"):
        if (
            row["fragment"],
            row["histidine"],
            row["solvent"],
            row["seed"],
            row["arm"],
        ) != ("triad", "delta", "water", "271828", "shell") or row["ligand"] not in (
            "H",
            "NH2",
            "OCH3",
        ):
            continue
        rec, name, i = row["receptor"], row["ligand"], int(row["angle_deg"]) // 60
        for term, label, idx in (
            ("complex", name, i),
            ("pocket", "", i),
            ("ligand_water", name, 0),
            ("ligand_vacuum", name, 0),
        ):
            reference[rec, label, term, idx] = float(row[term])
    assert energies.keys() == reference.keys() and len(energies) == 60
    rows = [
        dict(
            receptor=r,
            ligand=l,
            term=t,
            angle_deg=60 * i,
            tblite=reference[k],
            legacy=energies[k],
            delta=energies[k] - reference[k],
        )
        for k in energies
        for r, l, t, i in (k,)
    ]
    write(results / "independent_quantum.csv", rows)
    if any(abs(row["delta"]) > 0.02 for row in rows):
        raise ValueError(
            "Independent xTB total-energy difference exceeds 0.02 kcal/mol"
        )
    print(
        "Independent legacy single points",
        len(rows),
        "max |E difference|",
        max(abs(row["delta"]) for row in rows),
    )


def audit_edges(prepared, results):
    degrees = defaultdict(dict)
    for r in read(results / "energies.csv"):
        if (r["fragment"], r["histidine"], r["solvent"]) != ("triad", "delta", "water"):
            continue
        key = (r["receptor"], r["seed"], r["arm"], r["ligand"])
        i = int(r["angle_deg"])
        assert i not in degrees[key]
        degrees[key][i] = r
    assert len(degrees) == 80
    recomputed = {}
    for k, rot in degrees.items():
        assert sorted(rot) == list(range(0, 360, 60))

        # Direct partition ratio around each state's most favorable angle.
        # Subtracting a minimum keeps even strongly clashing geometries finite.
        def log_partition(term):
            energies = [float(row[term]) for row in rot.values()]
            low = min(energies)
            return -low / RT + math.log(
                sum(math.exp(-(e - low) / RT) for e in energies)
            )

        zero = rot[0]
        recomputed[k] = -RT * (
            log_partition("complex") - log_partition("pocket")
        ) - float(zero["ligand_water"])
        rec, seed, arm, l = k
        lines = (prepared / f"{rec}_{l}_{seed}_{arm}.pdbqt").read_text().splitlines()
        assert int(lines[0].split()[1]) == 1 and float(lines[1].split()[3]) == float(
            zero["native"]
        )
    count = 0
    for r in read(results / "summary.csv"):
        gaps = []
        for (rec, seed, arm, name), Gstrong in recomputed.items():
            if (
                arm == r["arm"]
                and (r["receptor"] in (rec, "both"))
                and name == r["strong"]
            ):
                gaps.append(recomputed[rec, seed, arm, r["weak"]] - Gstrong)
        assert (
            len(gaps) == int(r["matched_cells"])
            and abs(median(gaps) - float(r["triad_delta_water"])) < 1e-8
        )
        count += 1
    print("Independent angular sums and matched gap rows", count)


def audit_resolution(results):
    # Compare fine and coarse independently generated grids at their common angles.
    key = lambda row: (
        row["receptor"],
        row["ligand"],
        row["seed"],
        float(row["angle_deg"]),
    )
    coarse = {
        key(r): r
        for r in read(results / "energies.csv")
        if r["fragment"] == "triad"
        and r["solvent"] == "water"
        and r["histidine"] == "delta"
        and r["arm"] == "shell"
        and r["ligand"] in ("NH2", "OCH3")
    }
    fine = {key(r): r for r in read(results / "angular_raw.csv")}
    assert len(coarse) == 96 and len(fine) == 768
    for k, x in coarse.items():
        for name in ("pocket", "complex"):
            assert abs(float(fine[k][name]) - float(x[name])) < 0.002
    stored = read(results / "angular.csv")
    assert len(stored) == 64
    for row in stored:
        rec, lig, seed = row["receptor"], row["ligand"], row["seed"]
        spacing = float(row["step_deg"])
        values = [
            r
            for (rrec, rlig, rseed, phi), r in fine.items()
            if (rrec, rlig, rseed) == (rec, lig, seed) and phi % spacing == 0
        ]
        Z = lambda term: sum(
            math.exp(-(float(r[term]) - min(float(v[term]) for v in values)) / RT)
            for r in values
        )
        E = lambda term: min(float(r[term]) for r in values) - RT * math.log(Z(term))
        fixed = float(coarse[rec, lig, seed, 0]["ligand_water"])
        assert abs((E("complex") - E("pocket") - fixed) - float(row["transfer"])) < 1e-8
    print(
        "Fine-grid/cross-run common energies and independent nested partitions: 96 shared angles, 64 reductions"
    )


def neutron(path, results):
    """Read deposited alternative identities and H/D occupancies independently of the prepared protonation."""
    chosen = {
        ("HIS", 57): {"HD1", "DD1", "HE2", "DE2"},
        ("SER", 195): {"OG", "DG", "HG"},
    }
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  ") or line[21] != "A":
            continue
        res = (line[17:20], int(line[22:26]))
        name = line[12:16].strip()
        if name in chosen.get(res, ()):
            rows.append(
                dict(
                    residue=res[0],
                    number=res[1],
                    atom=name,
                    alternate=line[16].strip(),
                    occupancy=float(line[54:60]),
                )
            )
    if not rows:
        raise ValueError("No matching 5MO0 atoms")
    write(results / "neutron_sites.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prepared", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--xtb", type=Path)
    p.add_argument("--neutron", type=Path, help="optional 5MO0 coordinate file")
    a = p.parse_args()
    if a.neutron:
        neutron(a.neutron, a.results)
    if a.xtb:
        audit_legacy(a.prepared, a.results, a.xtb.resolve())
    audit_edges(a.prepared, a.results)
    audit_resolution(a.results)
