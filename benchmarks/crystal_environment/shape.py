"""Compare the existing typed bidirectional overlap, a Gaussian volume integral and disorder sensitivity."""

import argparse
import itertools
import json
from collections import defaultdict

import numpy as np
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

from density import transform
from evaluate import CHAIN, ROOT, models, read, save
from packing import atoms

# Bondi radii, Angstroms; neutral elemental shape only.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "F": 1.47, "CL": 1.75}
I3 = np.eye(3)


def transfer(path, ref, rec, donor, target):
    if target != "t4l":
        return transform(path, ref, rec, donor)

    def backbone(p):
        return {
            int(a.number): a.xyz
            for a in atoms(p)
            if a.chain == "A"
            and a.polymer
            and a.name == "CA"
            and a.alt in ("", "A")
            and a.number.isdigit()
            and 80 <= int(a.number) <= 160
            and int(a.number) not in range(107, 116)
        }

    left, right = backbone(path), backbone(ref)
    ids = sorted(left.keys() & right.keys())
    a = np.array([left[k] for k in ids])
    b = np.array([right[k] for k in ids])
    rotation, _ = Rotation.align_vectors(b - b.mean(0), a - a.mean(0))
    return a.mean(0), rotation.as_matrix().T, b.mean(0)


def masses_variances(elements):
    m = np.array([4 * np.pi * VDW[e] ** 3 / 3 for e in elements])
    return m, m ** (2 / 3) / (2 * np.pi)


def kernel(m1, v1, m2, v2):
    covariance = v1[:, None] + v2[None, :]
    norm = m1[:, None] * m2[None, :] / np.sqrt(np.linalg.det(covariance))
    return np.linalg.inv(covariance), norm


def integral(x, y, k):
    delta = x[:, None] - y
    inv, norm = k
    exponent = np.einsum("...i,...ij,...j->...", delta, inv, delta)
    return np.sum(norm * np.exp(-0.5 * exponent))


class Template:
    def __init__(self, elems, atoms_ref, pdb):
        self.qelems = np.asarray(elems)
        self.x = np.array([a.xyz for a in atoms_ref])
        self.elements = np.array([a.element for a in atoms_ref])
        anis = {}
        for l in pdb.read_text().splitlines():
            if l.startswith("ANISOU"):
                u = (
                    np.array([int(l[k : k + 7]) for k in (28, 35, 42, 49, 56, 63)])
                    * 1e-4
                )
                u11, u22, u33, u12, u13, u23 = u
                anis[l[21], l[22:27].strip(), l[12:16].strip(), l[16].strip()] = (
                    np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
                )
        observed = [
            anis.get((a.chain, a.number, a.name, a.alt), a.b * I3 / (8 * np.pi**2))
            for a in atoms_ref
        ]
        if np.min(np.linalg.eigvalsh(observed)) < -0.0001:
            raise ValueError("Nonpositive deposited displacement tensor")
        qm, qv = masses_variances(elems)
        dm, dv = masses_variances(self.elements)
        qv = qv[:, None, None] * I3
        dv = dv[:, None, None] * I3
        self.qkernel = kernel(qm, qv, qm, qv)
        self.cross = []
        self.self_ref = []
        for width in (dv, dv + observed):
            self.cross.append(kernel(qm, qv, dm, width))
            self.self_ref.append(integral(self.x, self.x, kernel(dm, width, dm, width)))

    def scores(self, xyz):
        distances = cdist(xyz, self.x, "sqeuclidean")
        typed = np.where(
            self.qelems[:, None] == self.elements[None, :], distances, np.inf
        )
        one = float(np.exp(-0.5 * distances.min(1)).mean())
        loose = float(np.exp(-distances.min(1) / (2 * 1.5**2)).mean())
        both = float(
            (
                np.exp(-distances.min(1) / (2 * 1.5**2)).mean()
                + np.exp(-distances.min(0) / (2 * 1.5**2)).mean()
            )
            / 2
        )
        old = float(
            (
                np.exp(-typed.min(1) / (2 * 1.5**2)).mean()
                + np.exp(-typed.min(0) / (2 * 1.5**2)).mean()
            )
            / 2
        )
        qself = integral(xyz, xyz, self.qkernel)
        ts = []
        for k, dself in zip(self.cross, self.self_ref, strict=True):
            xy = integral(xyz, self.x, k)
            ts.append(float(xy / (qself + dself - xy)))
        return dict(
            one1=one,
            one15=loose,
            both15=both,
            existing=old,
            volume=ts[0],
            volume_u=ts[1],
        )


def evaluate(ca, brd, t4, ca_runs, brd_runs, t4_runs, out):
    first = read(out / "candidates.csv")
    pools = {
        "brd4": [r for r in first if r["target"] == "brd4"],
        "hcaii": [r for r in first if r["target"] == "hcaii"],
    }
    pools["t4l"] = [
        dict(
            **r,
            pool=r["arm"],
            slot=r["arm"],
            sign=r["sign_ok"],
            intra=0,
            factor=0,
            energy=r["native_affinity"],
            native=r["raw"],
        )
        for r in read(t4_runs / "modes.csv")
    ]
    donorlist = []
    cells = []
    dm = []
    for target, source, work in [
        ("hcaii", ca, ca_runs),
        ("brd4", brd, brd_runs),
        ("t4l", t4, t4_runs),
    ]:
        features = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        metadata = json.loads(
            (
                source
                / "prepared"
                / ("metadata.json" if target == "hcaii" else "manifest.json")
            ).read_text()
        )
        refs = json.loads(
            (
                ROOT
                / "benchmarks/mapping_ensemble/results"
                / ("hcaii_declared" if target == "hcaii" else target)
                / "graphs.json"
            ).read_text()
        )
        fingerprints = {
            k: features.GetFingerprint(Chem.MolFromSmiles(m["smiles"]))
            for k, m in metadata.items()
        }
        donor_map = {}
        for lig, rec in refs.items():
            sim = lambda k: DataStructs.TanimotoSimilarity(
                fingerprints[lig], fingerprints[k]
            )
            ref = rec["donor"]
            control = min(
                (k for k in refs if k not in (lig, ref)), key=lambda k: (sim(k), k)
            )
            donor_map[lig] = dict(scaffold=ref, distant=control)
            for kind, d in donor_map[lig].items():
                donorlist.append(
                    dict(target=target, ligand=lig, kind=kind, donor=d, tanimoto=sim(d))
                )
        byfile = defaultdict(list)
        for r in pools[target]:
            byfile[r["receptor"], r["ligand"], r["seed"], r["slot"]].append(r)
        comparisons = defaultdict(list)
        cache = {}
        label_key = "Ki_nM" if target == "hcaii" else "kd_uM"
        for (rec, lig, seed, slot), group in byfile.items():
            filename = work / f"{rec}_{lig}_{seed}_{slot}.pdbqt"
            poses = models(
                filename, None if target == "hcaii" else metadata[lig]["heavy_serials"]
            )
            coords, elements = poses[1]
            for kind, donor in donor_map[lig].items():
                key = (rec, lig, donor)
                pdb = source / "raw" / (donor + ".pdb")
                if key not in cache:
                    component = (
                        metadata[donor].get("ccd") or metadata[donor]["component"]
                    )
                    group_ref = [
                        a
                        for a in atoms(pdb)
                        if a.residue == component
                        and a.chain == CHAIN.get(donor, "A")
                        and (target != "hcaii" or a.number == "303")
                    ]
                    if len({a.number for a in group_ref}) != 1 or any(
                        a.alt for a in group_ref
                    ):
                        raise ValueError((donor, group_ref))
                    cache[key] = (
                        Template(elements, group_ref, pdb),
                        transfer(
                            source / "raw" / (rec + ".pdb"), pdb, rec, donor, target
                        ),
                    )
                template, (src, rot, dst) = cache[key]
                for r in group:
                    n = int(r["model"])
                    scores = template.scores((poses[n][0] - src) @ rot + dst)
                    comparisons[rec, lig, seed, r["pool"], kind].append((r, scores))
                    dm.append(
                        dict(
                            target=target,
                            receptor=rec,
                            ligand=lig,
                            seed=seed,
                            slot=slot,
                            model=n,
                            kind=kind,
                            donor=donor,
                            **scores,
                        )
                    )
        for (rec, lig, seed, pool, kind), group in comparisons.items():
            for metric in group[0][1]:

                def raw(r, s):
                    return float(r["native"]) - 6 * s[metric]

                for selection in ("geometric", "energy"):
                    r, score = min(
                        group,
                        key=lambda item: (
                            (
                                -item[1][metric]
                                if selection == "geometric"
                                else raw(*item)
                            ),
                            float(item[0]["native"]),
                            item[0]["slot"] != "common",
                            int(item[0]["model"]),
                        ),
                    )
                    energy = float(r["energy"])
                    if target == "hcaii":
                        unbound = float(
                            min(group, key=lambda it: float(it[0]["native"]))[0][
                                "intra"
                            ]
                        )
                        energy = (float(r["native"]) - unbound) * float(r["factor"])
                    cells.append(
                        dict(
                            target=target,
                            receptor=rec,
                            ligand=lig,
                            seed=seed,
                            pool=pool,
                            kind=kind,
                            metric=metric,
                            selection=selection,
                            model=r["model"],
                            slot=r["slot"],
                            rmsd=float(r["rmsd"]),
                            joint=int(float(r["rmsd"]) <= 2 and int(r["sign"])),
                            assay_value=metadata[lig][label_key],
                            energy=energy,
                        )
                    )
    save(out / "shape_donors.csv", donorlist)
    save(out / "shape_terms.csv", dm)
    save(out / "shape_cells.csv", cells)
    summary = []
    keys = sorted(
        {
            (r["target"], r["pool"], r["kind"], r["metric"], r["selection"])
            for r in cells
        }
    )
    for k in keys:
        subset = [
            r
            for r in cells
            if tuple(
                r[col] for col in ("target", "pool", "kind", "metric", "selection")
            )
            == k
        ]
        pairs = defaultdict(list)
        for block in sorted({(r["receptor"], r["seed"]) for r in subset}):
            local = sorted(
                (r for r in subset if (r["receptor"], r["seed"]) == block),
                key=lambda r: r["assay_value"],
            )
            for a, b in itertools.combinations(local, 2):
                if b["assay_value"] >= 3 * a["assay_value"]:
                    pairs[a["ligand"], b["ligand"]].append(
                        round((b["energy"] - a["energy"]) * 1e9)
                    )
        credit = lambda g: 1 if g > 0 else 0.5 if g == 0 else 0
        summary.append(
            dict(
                target=k[0],
                pool=k[1],
                kind=k[2],
                metric=k[3],
                selection=k[4],
                n=len(subset),
                joint=sum(r["joint"] for r in subset),
                median=np.median([r["rmsd"] for r in subset]),
                pair_credit=sum(credit(g) for gs in pairs.values() for g in gs),
                pair_n=sum(len(gs) for gs in pairs.values()),
                unique_credit=sum(credit(np.median(gs)) for gs in pairs.values()),
                unique_n=len(pairs),
            )
        )
    save(out / "shape_summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("ca", "brd", "t4", "ca_runs", "brd_runs", "t4_runs", "out"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    evaluate(a.ca, a.brd, a.t4, a.ca_runs, a.brd_runs, a.t4_runs, a.out)
