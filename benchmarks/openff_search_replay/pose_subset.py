"""Audit and reproduce five human-protein cases from the earlier frozen PoseBusters list."""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem

CASES = ("7WPW_F15", "8EX2_Q2Q", "5SB2_1K2", "7VWF_K55", "7NR8_UOE")


def prepare(data, work):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "paired_crossdocking"))
    from run import coords, receptor

    for case in CASES:
        source = data / "posebusters_benchmark_set" / case
        lines = (source / (case + "_protein.pdb")).read_text().splitlines()
        ref = Chem.SDMolSupplier(str(source / (case + "_ligand.sdf")), removeHs=False)[
            0
        ]
        xyz = Chem.RemoveHs(ref).GetConformer().GetPositions()
        excluded = [
            l
            for l in lines
            if l.startswith(("HETATM", "ATOM  "))
            and l[76:78].strip() not in ("H", "D")
            and l[17:20] != "HOH"
            and (l.startswith("HETATM") or l[21] != "A")
        ]
        closest = min(
            (
                (float(np.min(np.linalg.norm(xyz - coords(l), axis=1))), l[17:26])
                for l in excluded
            ),
            default=(None, None),
        )
        out = work / case
        out.mkdir(exist_ok=True, parents=True)
        details = {
            "closest_excluded_nonwater_heavy": closest[0],
            "closest_excluded_group": closest[1],
        }
        if closest[0] is not None and closest[0] < 8:
            details["status"] = "excluded_nearby_unmodeled_group"
        else:
            details.update(receptor(lines, out, case, xyz))
            details["status"] = "prepared"
        (out / "repair.json").write_text(json.dumps(details, indent=2) + "\n")


def dock(data, work, prior_source, mixed, soft):
    from replay import load_prior

    prior = load_prior(prior_source)
    prior.ARMS["soft"] = ("--soft_receptor_search",)
    rows = []
    for case in CASES:
        source = data / "posebusters_benchmark_set" / case
        out = work / case
        details = json.loads((out / "repair.json").read_text())
        if details["status"] != "prepared":
            continue
        ref = prior.read_sdf(source / (case + "_ligand.sdf"))
        pdb = source / (case + "_protein.pdb")
        inp = out / "ligand.pdbqt"
        prior.prepare_ligand(Chem.MolToSmiles(Chem.RemoveHs(ref)), inp)
        for seed in prior.SEEDS:
            for arm, binary in (("native", mixed), ("mixed", mixed), ("soft", soft)):
                pose, scores, vals = prior.docking(
                    binary,
                    out / "receptor.pdbqt",
                    inp,
                    prior.coords(ref).mean(axis=0),
                    out,
                    seed,
                    arm,
                    4,
                    3000,
                )
                stats = prior.pose_metrics(pose, scores, ref, pdb, out)
                rows.append(dict(case=case, seed=seed, arm=arm, **vals, **stats))
    with (work / "subset.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["prepare", "dock"])
    for name in ("data", "work"):
        p.add_argument(name, type=Path)
    for name in ("prior-source", "mixed-vina", "soft-vina"):
        p.add_argument("--" + name, type=Path)
    a = p.parse_args()
    if a.mode == "prepare":
        prepare(a.data, a.work)
    else:
        dock(
            a.data,
            a.work,
            a.prior_source,
            a.mixed_vina.resolve(),
            a.soft_vina.resolve(),
        )
