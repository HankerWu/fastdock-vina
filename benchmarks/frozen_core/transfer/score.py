"""Apply the pre-existing frozen ranker, then measure receptor-frame pose RMSD."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "contacts"))
from evaluate import atoms, native_components, profile, shape, tanimoto

sys.path.insert(0, str(HERE.parents[1] / "geometry_comparison"))
from pose_run import load_meta, align_with_numbering


def rank(pool, rec, donor):
    # These arguments contain only the generated query poses and receiver-side inputs.
    known = profile(rec, donor)
    rows = []
    for i, block in enumerate(pool.split("ENDMDL"), 1):
        if not block.strip():
            continue
        x = atoms(block)
        score, divisor = native_components(block)
        a, b = [tanimoto(obs, ref) for obs, ref in zip(profile(rec, x), known)]
        terms = [t for t in (a, b) if t is not None]
        overlap = sum(terms) / len(terms) if terms else 0.0
        element_overlap = shape(x, donor)
        rows.append(
            dict(
                model=i,
                divisor=divisor,
                native=score,
                element=score - 6 * element_overlap / divisor,
                contact=score - 6 * overlap / divisor,
                general=score - 6 * (a or 0.0) / divisor,
                contact_only=-overlap,
                element_overlap=element_overlap,
                general_overlap=a,
                polar_overlap=b,
            )
        )
    return rows


def experiment(work, out):
    sites = {
        m["pdb"]: load_meta(work, m["pdb"])
        for m in json.loads((work / "manifest.json").read_text())
        if m["status"] == "ready"
    }
    modes, cells = [], []
    for path in sorted((work / "pools").glob("*.pdbqt")):
        q, r, seed = path.stem.split("-")
        qm, rm = sites[q], sites[r]
        assert q != r and qm["canonical_smiles"] != rm["canonical_smiles"]
        ligand = atoms((work / r / "ligand.pdbqt").read_text())
        ids = rm["serial"]
        donor = (
            np.asarray(rm["mapped"][0]),
            np.array([ligand[1][ligand[2].index(i)] for i in ids]),
            ids,
        )
        rec = atoms((work / r / "receptor.pdbqt").read_text())
        pool = path.read_text()
        rows = rank(pool, rec, donor)
        arms = ("native", "element", "contact", "general", "contact_only")
        selected = {
            arm: min(rows, key=lambda x: (x[arm], x["model"]))["model"] for arm in arms
        }
        rot, tran, alignment, count, offset = align_with_numbering(qm, rm)
        ref = np.asarray(qm["mapped"]) @ rot + tran
        case = dict(target=qm["target"], query=q, receiver=r, seed=int(seed))
        blocks = [b for b in pool.split("ENDMDL") if b.strip()]
        for block, row in zip(blocks, rows):
            xyz, _, serial = atoms(block)
            xyz = xyz[[serial.index(i) for i in qm["serial"]]]
            row["rmsd_A"] = float(
                np.sqrt(np.mean(((ref - xyz) ** 2).sum(axis=-1), axis=-1)).min()
            )
            modes.append(dict(**case, **row))
        cell = dict(
            **case,
            ca_count=count,
            ca_alignment_rmsd_A=alignment,
            residue_offset=offset,
            modes=len(rows),
            pool_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            best_A=min(row["rmsd_A"] for row in rows),
        )
        for arm, model in selected.items():
            cell[arm + "_rank"] = model
            cell[arm + "_A"] = rows[model - 1]["rmsd_A"]
        cells.append(cell)
    assert len(cells) == 8
    for name, rows in [("modes", modes), ("cells", cells)]:
        with (out / (name + ".csv")).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0])
            w.writeheader()
            w.writerows(rows)
    print(
        len(cells),
        "frozen pools",
        len(modes),
        "models",
        {
            a: sum(c[a + "_A"] <= 2 for c in cells)
            for a in ("native", "element", "contact", "general", "contact_only", "best")
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--out", type=Path, default=HERE)
    a = p.parse_args()
    experiment(a.work, a.out)
