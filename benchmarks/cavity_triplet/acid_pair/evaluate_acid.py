"""Read-only geometry analysis using original PDB, CCD and prepared microstate chemistry."""

from pathlib import Path
import argparse
import csv
import json
import sys
import tarfile

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from run_acid import BENCH, CASES, digest, table

sys.path.insert(0, str(BENCH / "official_audit"))
from verify_official import CODES, ligand, neutral_mol, transform

sys.path.insert(0, str(BENCH / "human_panel"))
from verify import reference


def xyz(block, ids):
    mapped = {
        int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
        for l in block.splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }
    return np.array([mapped[i] for i in ids])


def summary(p, out):
    out.mkdir(parents=True, exist_ok=True)
    dep = BENCH / "human_panel/depositions.tar.xz"
    with tarfile.open(dep) as archive:
        for rec, component in CODES.items():
            for n in (rec + ".pdb", component + ".cif"):
                archive.extract("raw/" + n, p, filter="data")
    refs = {r: reference(p, r, c) for r, c in CODES.items()}
    details = []
    cells = []
    stages = []
    chemical = []
    charge = []
    read = lambda fn: list(csv.DictReader(fn.open()))
    costs = {r["stem"]: r for r in read(p / "cost.csv")}
    for q, rec in CASES.items():
        rt = transform(q, rec, refs)
        ref = refs[q][2]
        site = rt(ref.GetConformer().GetPositions())
        for state in ("HA", "A"):
            tag = q + "-" + state
            start = (p / "inputs" / (tag + ".pdbqt")).read_text()
            m, _ = ligand(start)
            mapping = {}
            for l in start.splitlines():
                if l.startswith("REMARK SMILES IDX "):
                    v = list(map(int, l.split()[3:]))
                    mapping.update(zip(v[::2], v[1::2], strict=True))
            ids = [mapping[i + 1] for i in range(m.GetNumAtoms())]
            src = xyz(start, ids)
            neutral = neutral_mol(m)
            perms = ref.GetSubstructMatches(
                neutral, uniquify=False, useChirality=True, maxMatches=10001
            )
            assert 0 < len(perms) < 10001
            fixed = site[np.array(perms)]
            ((_, a, b),) = neutral.GetSubstructMatches(
                Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
            )
            swapped = fixed.copy()
            swapped[:, [a, b]] = fixed[:, [b, a]]
            equivalent = np.concatenate([fixed, swapped])
            chosen = fixed if state == "HA" else equivalent

            def rms(x, y):
                return np.sqrt(
                    np.min(
                        np.mean(
                            np.sum((x[..., None, :, :] - y) ** 2, axis=-1), axis=-1
                        ),
                        axis=-1,
                    )
                )

            edges = np.array(
                [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
            )
            lengths = lambda v: np.linalg.norm(v[edges[:, 0]] - v[edges[:, 1]], axis=1)
            stereo = lambda x: [
                np.sign(
                    np.linalg.det(
                        x[[j.GetIdx() for j in m.GetAtomWithIdx(i).GetNeighbors()][:3]]
                        - x[i]
                    )
                )
                for i, c in Chem.FindMolChiralCenters(m, includeUnassigned=False)
            ]
            initial_stereo = stereo(src)
            receiver = (p / "inputs" / rec / "receptor.pdbqt").read_text()
            heavy = [
                int(l[6:11])
                for l in receiver.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
                and l.split()[-1] not in ("H", "HD", "HS")
            ]
            tree = cKDTree(xyz(receiver, heavy))
            meta = dict(query=q, receiver=rec, state=state, seed=91813)
            qsum = sum(
                float(l[70:76])
                for l in start.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            )
            charge.append(
                dict(
                    **meta,
                    net_formal_charge=Chem.GetFormalCharge(m),
                    pdbqt_charge_e=qsum,
                    heavy_atoms=m.GetNumAtoms(),
                    graph_mappings=len(chosen),
                    torsdof=int(
                        next(
                            l.split()[1]
                            for l in start.splitlines()
                            if l.startswith("TORSDOF")
                        )
                    ),
                )
            )
            for arm in (
                "base",
                "high",
                "matched",
                "sham",
                "official",
                "official-double",
            ):
                stem = tag + "_" + arm
                pool = p / "raw" / (stem + ".pdbqt")
                assert digest(pool) == costs[stem]["pool_sha256"]
                errors = []
                for rank, block in enumerate(pool.read_text().split("ENDMDL")[:-1], 1):
                    x = xyz(block, ids)
                    r = float(rms(x, chosen))
                    errors.append(r)
                    v = float(
                        next(
                            l.split()[3]
                            for l in block.splitlines()
                            if l.startswith("REMARK VINA RESULT:")
                        )
                    )
                    details.append(
                        dict(
                            **meta,
                            arm=arm,
                            rank=rank,
                            rmsd_A=r,
                            fixed_O_A=float(rms(x, fixed)),
                            exchange_O_A=float(rms(x, equivalent)),
                            score=v,
                        )
                    )
                    chemical.append(
                        dict(
                            **meta,
                            arm=arm,
                            rank=rank,
                            heavy_bond_drift_A=max(abs(lengths(x) - lengths(src))),
                            centres_checked=len(initial_stereo),
                            inversions=sum(
                                u != v
                                for u, v in zip(initial_stereo, stereo(x), strict=True)
                            ),
                            closest_receptor_A=tree.query(x)[0].min(),
                        )
                    )
                cells.append(
                    dict(
                        **meta,
                        arm=arm,
                        modes=len(errors),
                        top1_A=errors[0],
                        best5_A=min(errors[:5]),
                        oracle_A=min(errors),
                        regret_A=errors[0] - min(errors),
                    )
                )
            name = tag + "_91813_t1_s1"
            with (
                np.load(p / "sources" / (name + ".npz")) as s,
                np.load(p / "candidates" / (name + ".npz")) as c,
            ):
                order = [s["serials"].tolist().index(i) for i in ids]
                for arm in ("matched", "sham"):
                    start_xyz = s["xyz"]
                    rot = c[arm + "_rotation"]
                    posed = (start_xyz - start_xyz.mean(0)) @ rot.transpose(
                        0, 2, 1
                    ) + c["center"][:, None, :]
                    vals = rms(posed[:, order], chosen)
                    index = int(c[arm + "_best"])
                    log = (
                        (p / "raw" / (tag + "_" + arm + ".stages"))
                        .read_text()
                        .splitlines()
                    )
                    assert list(map(int, log[0][1:].split())) == list(s["serials"])
                    events = {
                        l.split()[0]: np.array(l.split()[3:], float).reshape(-1, 3)[
                            order
                        ]
                        for l in log[1:]
                        if l.split()[0] in ("selected", "hunt", "second")
                    }
                    np.testing.assert_allclose(
                        posed[index, order], events["selected"], atol=1e-9
                    )
                    hist = (
                        (p / "raw" / (tag + "_" + arm + ".history"))
                        .read_text()
                        .splitlines()
                    )
                    stages.append(
                        dict(
                            **meta,
                            arm=arm,
                            source_A=float(rms(start_xyz[order], chosen)),
                            candidate_count=len(vals),
                            generator_best_A=min(vals),
                            generated_le2=sum(vals <= 2),
                            generated_le3=sum(vals <= 3),
                            selected_A=vals[index],
                            hunt_A=float(rms(events["hunt"], chosen)),
                            second_A=float(rms(events["second"], chosen)),
                            first_eviction_step=next(
                                (
                                    l.split()[0]
                                    for l in hist
                                    if "1" not in l.split()[4:]
                                ),
                                "",
                            ),
                        )
                    )
    for name, group in (
        ("modes", details),
        ("cells", cells),
        ("stages", stages),
        ("chemistry", chemical),
        ("microstates", charge),
    ):
        table(out / (name + ".csv"), group)
    sets = []
    for arm in ("base", "matched", "sham", "high", "official", "official-double"):
        for state in ("HA", "A", "all"):
            g = [
                r
                for r in cells
                if r["arm"] == arm and (r["state"] == state or state == "all")
            ]
            top = np.array([r["top1_A"] for r in g])
            sets.append(
                dict(
                    arm=arm,
                    state=state,
                    n=len(g),
                    mean=top.mean(),
                    **dict(
                        zip(
                            ("q25", "median", "q75", "p90", "p95"),
                            np.percentile(top, [25, 50, 75, 90, 95]),
                        )
                    ),
                    best5=np.mean([r["best5_A"] for r in g]),
                    oracle=np.mean([r["oracle_A"] for r in g]),
                    regret=np.mean([r["regret_A"] for r in g]),
                    severe_gt5=sum(top > 5),
                    severe_ge8=sum(top >= 8),
                    severe_gt10=sum(top > 10),
                    empty=sum(r["modes"] == 0 for r in g),
                )
            )
    table(out / "distributions.csv", sets)
    (out / "reference.json").write_text(
        json.dumps(
            dict(
                depositions_archive_sha256=digest(dep),
                donor_ids=CODES,
                source="CCD atom, charge and bond records; deposited chain A ligand coordinates; 14-A common-pocket C-alpha Kabsch fit",
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    summary(args.work.resolve(), args.out.resolve())
