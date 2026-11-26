"""Independent fixed-coordinate geometry and affinity label/credit audit."""

import argparse
import json
import lzma
import xlrd
from itertools import combinations
from pathlib import Path
import re
import tarfile

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from rdkit import Chem

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "geometry_comparison" / "affinity"


def atom_records(text):
    out = []
    for line in text.splitlines():
        if line[:6].strip() in ("ATOM", "HETATM") and line.split()[-1] not in (
            "H",
            "HD",
            "HS",
        ):
            out.append(
                (int(line[6:11]), *map(float, (line[30:38], line[38:46], line[46:54])))
            )
    return np.array(out)


def read_archive(path):
    with tarfile.open(path) as archive:
        return {
            m.name: archive.extractfile(m).read().decode()
            for m in archive
            if m.isfile() and m.name.endswith(".pdbqt")
        }


def topology(prepared):
    mol = Chem.MolFromSmiles(
        re.search(r"^REMARK SMILES (?!IDX)(.+)$", prepared, re.M)[1]
    )
    idx = []
    for line in re.findall(r"^REMARK SMILES IDX (.+)$", prepared, re.M):
        idx.extend(map(int, line.split()))
    to_serial = dict(zip((i - 1 for i in idx[::2]), idx[1::2]))
    to_rdkit = {v: k for k, v in to_serial.items()}
    assert len(to_serial) == mol.GetNumAtoms()
    cuts = []
    for a, b in re.findall(r"^BRANCH\s+(\d+)\s+(\d+)\s*$", prepared, re.M):
        a, b = int(a), int(b)
        changed = Chem.RWMol(mol)
        changed.RemoveBond(to_rdkit[a], to_rdkit[b])
        frags = Chem.GetMolFrags(changed)
        assert len(frags) == 2
        if min(map(len, frags)) > 1:
            cuts.append((a, b, *[[to_serial[i] for i in frag] for frag in frags]))
    return cuts


def geometric_count(model, receptor, cuts):
    atoms = atom_records(model)
    serial_to_i = {int(k): i for i, k in enumerate(atoms[:, 0])}
    xyz = atoms[:, 1:]
    # Explicit dense Euclidean pair distances and projection identity; no production KD-tree or cut-tree.
    lo, hi = xyz.min(axis=0) - 4.01, xyz.max(axis=0) + 4.01
    near = receptor[np.all((receptor >= lo) & (receptor <= hi), axis=1)]
    if not len(near):
        return 0
    witnesses = set()
    for idx in np.flatnonzero(cdist(xyz, near).min(axis=1) <= 4.0):
        witnesses.add(int(atoms[idx, 0]))
    count = 0
    for a, b, left, right in cuts:
        origin = xyz[serial_to_i[a]]
        axis = xyz[serial_to_i[b]] - origin
        unit = axis / np.linalg.norm(axis)

        def anchored(side):
            for s in witnesses.intersection(side):
                d = xyz[serial_to_i[s]] - origin
                perpendicular_sq = np.dot(d, d) - np.dot(d, unit) ** 2
                if perpendicular_sq >= 1.0:
                    return True
            return False

        count += anchored(left) and anchored(right)
    return count


def verify_geometry(data, pooldict, substrates, target=None):
    n = 0
    data = data.sort_index()
    for path, text in pooldict.items():
        if target:
            state = path.split("/")[0]
            stem = Path(path).stem
            if "_" not in stem:
                continue
            compound, seed = stem.split("_")
            prepared = substrates[state + "/" + compound + ".pdbqt"]
            key = (state, compound, int(seed))
        else:
            state, compound, seed, opt = Path(path).stem.split("_")
            prepared = substrates["prepared/" + compound + ".pdbqt"]
            key = (opt, state, compound, int(seed))
        recpath = (state + "/" if target else "prepared/") + state + ".pdbqt"
        receptor = atom_records(substrates[recpath])[:, 1:]
        cuts = topology(prepared)
        group = data.loc[key].sort_index()
        for rank, mode in enumerate(
            re.findall(r"MODEL \d+\n(.+?)ENDMDL", text, re.S), 1
        ):
            expected = group.loc[rank]
            assert len(cuts) == expected.eligible_rotors
            assert geometric_count(mode, receptor, cuts) == expected.constrained, (
                key,
                rank,
            )
            n += 1
    assert n == len(data)
    return n


def audit_labels(source, out):
    panel = pd.read_csv(source / "panel.csv", dtype={"compound": str})
    with tarfile.open(source / "inputs.tar.xz") as archived:
        table = xlrd.open_workbook(
            file_contents=archived.extractfile("raw/table4.xls").read()
        ).sheet_by_index(0)
        kinase = {
            table.cell_value(i, 2): table.row_values(i) for i in range(1, table.nrows)
        }
        bace = {}
        for m in archived:
            if m.name.endswith("_Ki.json"):
                for activity in json.load(archived.extractfile(m))["activities"]:
                    if (
                        activity["standard_type"],
                        activity["standard_relation"],
                        activity["standard_units"],
                    ) == ("Ki", "=", "nM"):
                        bace[
                            activity["assay_chembl_id"], activity["molecule_chembl_id"]
                        ] = float(activity["standard_value"])
        for r in panel.itertuples():
            if r.endpoint == "Ki":
                assert bace[r.assay, r.compound] == r.value_nM
            else:
                assert (
                    float(
                        kinase[r.assay.removeprefix("Davis-")][
                            table.row_values(0).index(r.name)
                        ]
                    )
                    == r.value_nM
                )
    actual = pd.read_csv(out / "pairs.csv.xz", dtype={"stronger": str, "weaker": str})
    poses = pd.read_csv(out / "poses.csv.xz", dtype={"compound": str})
    for cap in sorted(poses.budget.unique()):
        subset = poses.query("budget==@cap").set_index(
            ["optimizer", "receptor", "compound", "seed", "rank"]
        )
        verify_geometry(
            subset,
            read_archive(source / f"cap{cap}" / "outputs.tar.xz"),
            read_archive(source / "inputs.tar.xz"),
        )
    dims = ["budget", "optimizer", "policy", "score", "receptor", "compound", "seed"]
    selected = pd.read_csv(out / "selected.csv.xz", dtype={"compound": str}).set_index(
        dims
    )
    for score in ("native", "raw", "uniform", "contact"):
        for policy in ("all", "first"):
            rows = poses if policy == "all" else poses.query("rank==1")
            target = rows.loc[
                rows.groupby(["budget", "optimizer", "receptor", "compound", "seed"])[
                    score
                ].idxmin()
            ]
            sel = selected.reset_index().query("score==@score and policy==@policy")
            compare = target.merge(
                sel, on=["budget", "optimizer", "receptor", "compound", "seed"]
            )
            assert np.array_equal(compare["rank"], compare["selected"])
            assert np.allclose(compare[score], compare["value"])
    med = selected.reset_index().groupby(dims[:-1]).value.median().reset_index()
    keys = [
        "budget",
        "optimizer",
        "policy",
        "score",
        "target",
        "assay",
        "receptor",
        "stronger",
        "weaker",
    ]
    expected = []
    for (cap, opt, policy, score), grp in med.groupby(dims[:4]):
        for (target, assay), labs in panel.groupby(["target", "assay"]):
            states = actual.loc[
                (actual.assay == assay) & (actual.receptor != "mean_receptors"),
                "receptor",
            ].unique()
            sc = (
                grp.loc[grp.receptor.isin(states)]
                .pivot(index="compound", columns="receptor", values="value")
                .loc[labs.compound]
            )
            sc["mean_receptors"] = sc.mean(axis=1)
            labs = labs.set_index("compound").value_nM.to_dict()
            for a, b in combinations(sorted(labs, key=labs.get), 2):
                if labs[b] < 3 * labs[a]:
                    continue
                for rec, gap in (sc.loc[b] - sc.loc[a]).items():
                    gap = round(gap, 10)
                    expected.append(
                        (
                            cap,
                            opt,
                            policy,
                            score,
                            target,
                            assay,
                            rec,
                            a,
                            b,
                            (np.sign(gap) + 1) / 2,
                        )
                    )
    check = pd.DataFrame(expected, columns=[*keys, "credit_ref"]).merge(
        actual, on=keys, validate="one_to_one"
    )
    assert (
        len(check) == len(expected) == len(actual)
        and (check.credit_ref == check.credit).all()
    )
    print(
        len(poses),
        "affinity poses validated with independent graph cuts and distances;",
        len(check),
        "main panel pair credits agree",
    )


def audit_fxa(out):
    with lzma.open(out / "fxa_activity.json.xz", "rt") as s:
        records = json.load(s)["activities"]
    labels = {}
    for r in records:
        assert (
            r["assay_chembl_id"],
            r["standard_relation"],
            r["standard_type"],
            r["standard_units"],
        ) == ("CHEMBL661804", "=", "Ki", "nM")
        mol = Chem.MolFromSmiles(r["canonical_smiles"])
        if any(
            tag == "?"
            for _, tag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        ):
            continue
        labels[r["molecule_chembl_id"]] = float(r["standard_value"])
    assert len(records) == 34 and len(labels) == 28
    poses = pd.read_csv(out / "fxa_poses.csv.xz")
    expected = (
        pd.read_csv(out / "fxa_pairs.csv.xz")
        .set_index(["policy", "score", "receptor", "stronger", "weaker"])
        .sort_index()
    )
    n = 0
    for policy in ("first", "all"):
        modes = poses if policy == "all" else poses.query("rank==1")
        for score in ("native", "raw", "uniform", "contact"):
            wide = (
                modes.groupby(["compound", "receptor", "seed"])[score]
                .min()
                .unstack("compound")
            )
            for state in ("1FJS", "2P16", "both"):
                block = wide if state == "both" else wide.loc[state]
                for a, b in combinations(sorted(labels, key=labels.get), 2):
                    if labels[b] < 3 * labels[a]:
                        continue
                    ref = expected.loc[policy, score, state, a, b]
                    gap = round(float(np.median(block[b] - block[a])), 10)
                    assert abs(ref.fold - labels[b] / labels[a]) < 1e-8
                    assert (
                        abs(ref.gap - gap) < 1e-8
                        and (np.sign(gap) + 1) / 2 == ref.credit
                    )
                    n += 1
    assert n == len(expected)
    print(
        n, "additional external FXa pair credits verified against original assay export"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    audit_labels(args.source, args.out)
    audit_fxa(args.out)
    blob = read_archive(args.out / "fxa_frozen.tar.xz")
    data = pd.read_csv(args.out / "fxa_poses.csv.xz").set_index(
        ["receptor", "compound", "seed", "rank"]
    )
    print(
        verify_geometry(data, blob, blob, target="fxa"),
        "additional external frozen FXa geometries independently agree",
    )
