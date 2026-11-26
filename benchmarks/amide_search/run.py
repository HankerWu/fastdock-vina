"""Equal-budget tests of a continuous amide rotor versus its two planar wells."""

import argparse
import csv
import itertools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from amide_rotamers import amides, rotate, lock_amide
from score_physical_pool import mapped_ligand_bonds

SEEDS = (314159, 271828)
FXA = ("CHEMBL172831", "CHEMBL151985", "CHEMBL355796")


def sources(base):
    for panel, folder, receptors in (
        ("ca", "fdv-second", ("6WQ4", "6WQ8")),
        ("brd", "fdv-brd4", ("3MXF", "4OGI", "4OGJ")),
    ):
        prep = base / folder / "prepared"
        for original in sorted(prep.glob("*_ligand.pdbqt")):
            lig = original.stem.split("_")[0]
            if not amides(original.read_text()):
                continue
            states = {"orig": original}
            if panel == "ca":
                states["mirror"] = (
                    base / "fdv-map-states-declared/inputs" / f"{lig}_reflected.pdbqt"
                )
            for rec in receptors:
                anchor = (
                    base / "fdv-map-states-declared/inputs" / f"{rec}_{lig}.anchor"
                    if panel == "ca"
                    else base / "fdv-map-brd4" / f"{rec}_{lig}.multi"
                )
                yield (
                    panel,
                    rec,
                    lig,
                    states,
                    prep / f"{rec}_receptor.pdbqt",
                    prep / f"{rec}_box.txt",
                    anchor,
                )
    with (ROOT / "benchmarks/affinity_partition/fxa_series.csv").open() as source:
        compounds = [
            r["ligand"] for r in csv.DictReader(source) if r["included"] == "1"
        ]
    for rec, d in (("1FJS", "fdv-fxa1"), ("2P16", "fdv-fxa")):
        for lig in compounds:
            prep = base / d
            repair = lig == "CHEMBL369480"
            yield (
                "fxa",
                rec,
                lig,
                {
                    "orig": (base / "fdv-ionpair/repair" if repair else prep)
                    / f"{lig}.pdbqt"
                },
                prep / f"{rec}.pdbqt",
                prep / "box.txt",
                (base / ("fdv-ionpair/repair" if repair else "fdv-ionpair/fxa"))
                / f"{rec}_{lig}.shell",
            )


def build(base, out):
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    tasks = []
    for panel, rec, lig, states, receptor, box, anchor in sources(base):
        for state, p in states.items():
            orig = p.read_text()
            site = amides(orig)[0] if panel != "fxa" or lig in FXA else None
            prefix = f"{panel}_{lig}_{state}"
            bonds = inputs / (prefix + ".bonds")
            bonds.write_text(
                "".join(f"{i} {j}\n" for i, j in mapped_ligand_bonds(orig))
            )
            for form in ("free", "near", "far") if site else ("free",):
                ligand = inputs / (prefix + "_" + form + ".pdbqt")
                ligand.write_text(
                    orig
                    if form == "free"
                    else lock_amide(
                        rotate(orig, site, alternate=form == "far", planar=True), site
                    )
                )
                for seed, slot in itertools.product(SEEDS, (0, 1)):
                    stem = f"{panel}_{rec}_{lig}_{state}_{seed}_{slot}_{form}"
                    tasks.append(
                        (
                            out / (stem + ".pdbqt"),
                            [receptor, ligand, box, anchor, bonds, seed + 1009 * slot],
                        )
                    )
    return tasks


def dock(binary, item, *, search_amide=None, rank_amide=None):
    dest, (rec, lig, box, anchor, bonds, seed) = item
    args = [
        binary,
        "--receptor",
        rec,
        "--ligand",
        lig,
        "--config",
        box,
        "--search_anchor",
        anchor,
        "--ligand_bonds",
        bonds,
        "--cpu",
        "1",
        "--exhaustiveness",
        "1",
        "--max_evals",
        "50000",
        "--num_modes",
        "9",
        "--pose_selection",
        "late",
        "--energy_range",
        "100",
        "--seed",
        str(seed),
        "--search_counters",
        "--out",
        dest,
    ]
    if search_amide is not None:
        args.extend(["--search_amide", search_amide])
    if rank_amide is not None:
        args.extend(["--rank_amide", rank_amide])
    with dest.with_suffix(".log").open("w") as f:
        subprocess.run(args, stdout=f, stderr=subprocess.STDOUT, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "base", "out"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    jobs = build(args.base, args.out)
    with ThreadPoolExecutor(args.jobs) as workers:
        list(workers.map(lambda j: dock(args.binary.resolve(), j), jobs))
    print(len(jobs), "fixed-budget one-trajectory searches")
