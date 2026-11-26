"""Same-assay size and charge controls for the native torsion divisor ablation."""

import argparse
import itertools
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from fxa_holdout import reported_scores
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "strain_transfer"))
from evaluate import SHARED, load, panels, write


def discovery(evidence):
    rows = []
    for (target, lig), c in panels(evidence):
        if c["receptor"] not in SHARED[target]:
            continue
        lines = c["path"].read_text().split("ENDMDL", 1)[0].splitlines()
        (smi,) = [
            l.split()[2]
            for l in lines
            if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
        ]
        mol = Chem.MolFromSmiles(smi)

        rows.append(
            dict(
                target=target,
                ligand=lig,
                receptor=c["receptor"],
                seed=int(c["seed"]),
                measurement=float(c.get("measurement", c.get("kd_uM"))),
                heavy=mol.GetNumHeavyAtoms(),
                charge=Chem.GetFormalCharge(mol),
                **reported_scores(lines),
            )
        )
    return rows


def summarize(rows):
    ligands = defaultdict(dict)
    for r in rows:
        ligands[r["target"]][r["ligand"]] = r
    cells = {(r["target"], r["receptor"], int(r["seed"]), r["ligand"]): r for r in rows}
    edges = []

    def credit(x):
        return 1 if x > 0 else 0.5 if x == 0 else 0

    for target, ligs in ligands.items():
        blocks = {(rec, seed) for t, rec, seed, lig in cells if t == target}
        for a, b in itertools.combinations(
            sorted(ligs, key=lambda l: (float(ligs[l]["measurement"]), l)), 2
        ):
            ratio = float(ligs[b]["measurement"]) / float(ligs[a]["measurement"])
            if ratio < 3:
                continue
            matched = int(ligs[a]["heavy"]) == int(ligs[b]["heavy"]) and int(
                ligs[a]["charge"]
            ) == int(ligs[b]["charge"])
            gaps = {}
            for arm in ("native", "no_divisor", "rotor"):
                # Integer milli-kcal differences preserve the precision of reported energies, including ties.
                raw = [
                    round(1000 * float(cells[target, rec, seed, b][arm]))
                    - round(1000 * float(cells[target, rec, seed, a][arm]))
                    for rec, seed in sorted(blocks)
                ]
                gaps[arm] = statistics.median(raw) / 1000
            gaps["size"] = int(ligs[a]["heavy"]) - int(ligs[b]["heavy"])
            edges.append(
                dict(
                    target=target,
                    strong=a,
                    weak=b,
                    ratio=ratio,
                    heavy_strong=ligs[a]["heavy"],
                    heavy_weak=ligs[b]["heavy"],
                    charge_strong=ligs[a]["charge"],
                    charge_weak=ligs[b]["charge"],
                    matched=int(matched),
                    **{arm + "_gap": v for arm, v in gaps.items()},
                    **{arm + "_credit": credit(v) for arm, v in gaps.items()},
                )
            )
    write(HERE / "results/size_pairs.csv", edges)
    summary = []
    for target in sorted(ligands):
        for subset in ("all", "equal_size_charge"):
            group = [
                r
                for r in edges
                if r["target"] == target and (subset == "all" or r["matched"])
            ]
            if not group:
                continue
            changes = [r["no_divisor_credit"] - r["native_credit"] for r in group]
            loo = []
            for lig in ligands[target]:
                left = [
                    d
                    for r, d in zip(group, changes)
                    if lig not in (r["strong"], r["weak"])
                ]
                if left:
                    loo.append(statistics.mean(left))
            summary.append(
                dict(
                    target=target,
                    subset=subset,
                    ligands=len(ligands[target]),
                    pairs=len(group),
                    **{
                        arm: sum(r[arm + "_credit"] for r in group)
                        for arm in ("native", "no_divisor", "size", "rotor")
                    },
                    improved=sum(d > 0 for d in changes),
                    worsened=sum(d < 0 for d in changes),
                    loo_delta_min=min(loo) if loo else None,
                    loo_delta_max=max(loo) if loo else None,
                )
            )
    write(HERE / "results/size_controls.csv", summary)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--evidence",
        type=Path,
        help="Regenerate development cells from archived raw, checksum-verified poses",
    )
    a = p.parse_args()
    orig = HERE / "results/discovery_cells.csv"
    if a.evidence:
        write(orig, discovery(a.evidence))
    first = load(HERE / "results/fxa_2P16_cells.csv")
    second = load(HERE / "results/fxa_1FJS_cells.csv")
    rotor_tables = load(HERE / "results/rotor_poses.csv") + load(
        HERE / "results/fxa_rotor_cells.csv"
    )

    def key(row):
        return tuple(row[k] for k in ("target", "ligand", "receptor", "seed"))

    rotor = {key(r): r["rotor"] for r in rotor_tables if r.get("model", "1") == "1"}
    merged = [{**r, "rotor": rotor[key(r)]} for r in load(orig) + first + second]
    summarize(
        merged
        + [
            {**r, "target": "fxa_" + r["receptor"]}
            for r in merged
            if r["target"] == "fxa"
        ]
    )


if __name__ == "__main__":
    main()
