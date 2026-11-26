"""Independent unit-weight MMFF94s annotation of original native winners."""

import argparse
import csv
import hashlib
import itertools
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools"))
import score_physical_pool as physical


def run_one(task):
    work, receptors, binary, comp, receiver, seed = task
    lig = work / (comp + ".pdbqt")
    pool = work / "pools" / f"{comp}-{receiver}-{seed}.pdbqt"
    annotation = physical.score_pool(
        binary,
        receptors / receiver / "receptor.pdbqt",
        pool,
        prepared_ligand=lig,
        mmff_affinity=True,
    )
    chem = next(
        r
        for r in json.loads((work / "chemistry.json").read_text())
        if r["compound"] == comp
    )
    if hashlib.sha256(lig.read_bytes()).hexdigest() != chem["pdbqt_sha256"]:
        raise ValueError("Prepared chemistry checksum mismatch")
    if (
        Chem.MolToSmiles(Chem.MolFromSmiles(chem["smiles"]))
        != annotation["strain"]["identity"]
    ):
        raise ValueError("Prepared graph differs from independently frozen identity")
    aff = annotation["mmff_affinity"]
    top = aff["native_rank_one"]
    first, *_ = physical.read_poses(pool)
    if first.model != top["model"] or top["native_rank"] != 1:
        raise ValueError("Native first model expected")
    bound = annotation["strain"]["pose_kcal_mol"][first.model]
    free = aff["reference"]["energy_kcal_mol"]
    inter = float(
        re.search(r"^REMARK INTER:\s+(-?[\d.]+)", pool.read_text(), re.MULTILINE)[1]
    )
    direct = aff["vina_flexibility_scale"] * (inter + bound - free)
    if abs(direct - top["affinity_score"]) > 0.002:
        raise ValueError("Direct INTER expression failed")
    return {
        "compound": comp,
        "receptor": receiver,
        "seed": seed,
        "native_score": first.native_score,
        "unit_score": top["affinity_score"],
        "saved_inter": inter,
        "saved_intra": first.native_intra,
        "saved_unbound": first.native_unbound,
        "bound_mmff": bound,
        "free_mmff": free,
        "reference_conformers": aff["reference"]["conformers_converged"],
        "scale": aff["vina_flexibility_scale"],
        "direct_unit_score": direct,
        "pool_min_unit_score": aff["minimum_across_emitted_models"]["affinity_score"],
        "selected_model": aff["minimum_across_emitted_models"]["model"],
        "pool_sha256": hashlib.sha256(pool.read_bytes()).hexdigest(),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("receptors", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    tasks = [
        (a.work, a.receptors, a.binary.resolve(), c, r, s)
        for c, r, s in itertools.product(
            ("niraparib", "olaparib", "talazoparib", "veliparib"),
            ("4UND", "7KK4"),
            (98207, 98213),
        )
    ]
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = list(pool.map(run_one, tasks))
    with a.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
