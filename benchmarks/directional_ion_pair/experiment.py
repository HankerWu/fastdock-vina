"""Use the receptor carboxylate to locate an aryl primary amidinium for search."""

import argparse
import csv
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from make_ion_pair_anchors import anchors, carboxylate

SEEDS = (271828, 314159)
ARMS = ("native", "point", "pair", "shell")


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def aspsite(path, chain):
    return carboxylate(path, chain, 189)


def prepare(work, rec, receptor, chain, input_ligands):
    acid = aspsite(receptor, chain)
    for ligand in input_ligands:
        for arm in ("point", "pair", "shell"):
            (work / f"{rec}_{ligand.stem}.{arm}").write_text(anchors(ligand, acid, arm))


def run(binary, output, inputs, cases, jobs):
    def one(case):
        rec, lig, seed, arm = case
        ligp, recp, config = inputs[rec, lig]
        base = output / f"{rec}_{lig}_{seed}_{arm}"
        command = [
            binary,
            "--ligand",
            ligp,
            "--receptor",
            recp,
            "--config",
            config,
            "--seed",
            str(seed),
            "--exhaustiveness",
            "4",
            "--cpu",
            "1",
            "--max_evals",
            "50000",
            "--num_modes",
            "9",
            "--pose_selection",
            "late",
            "--energy_range",
            "100",
            "--search_counters",
            "--out",
            base.with_suffix(".pdbqt"),
        ]
        if arm != "native":
            command.extend(["--search_anchor", output / f"{rec}_{lig}.{arm}"])
        with base.with_suffix(".log").open("w") as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(one, cases))


def fxa(a):
    manifest = table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
    compounds = [r["ligand"] for r in manifest if int(r["included"])]
    inputs = {}
    for rec, folder in (("1FJS", a.fjs), ("2P16", a.p16)):
        ligs = [folder / f"{lig}.pdbqt" for lig in compounds]
        prepare(a.out, rec, folder / f"{rec}.pdb", "A", ligs)
        for lig, file in zip(compounds, ligs, strict=True):
            inputs[rec, lig] = file, folder / f"{rec}.pdbqt", folder / "box.txt"
    run(
        a.binary.resolve(),
        a.out,
        inputs,
        itertools.product(("1FJS", "2P16"), compounds, a.seeds, a.arms),
        a.jobs,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("out", "fjs", "p16", "binary"):
        p.add_argument(key, type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--arms", nargs="+", choices=ARMS, default=["pair", "point"])
    p.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    a = p.parse_args()
    a.out.mkdir(exist_ok=True, parents=True)
    fxa(a)
