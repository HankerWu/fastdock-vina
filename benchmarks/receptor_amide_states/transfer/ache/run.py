"""Matched physical atom identities for a two-ligand AChE amine protonation counterfactual."""

import argparse, csv, json, re, subprocess, sys, time
from pathlib import Path
import numpy as np
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import paired

CODES = ("4EY6", "4EY7")
SEEDS = (71993, 72019)


def atoms(text):
    return {
        int(l[6:11]): (
            np.array([float(l[i : i + 8]) for i in (30, 38, 46)]),
            l.split()[-1],
        )
        for l in text.splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }


def extract_map(text, prefix):
    m = {}
    for l in text.splitlines():
        if l.startswith(prefix):
            a = list(map(int, l[len(prefix) :].split()))
            for i, j in zip(a[::2], a[1::2], strict=True):
                if i in m:
                    raise ValueError(f"Repeated map key {i}")
                m[i] = j
    return m


def branch_layout(text, serials):
    index = {s: i for i, s in enumerate(serials)}
    stack = []
    branches = []
    for line in text.splitlines():
        if line.startswith("BRANCH"):
            i, j = map(int, line.split()[1:])
            stack.append((index[i], index[j], []))
        elif line.startswith(("ATOM  ", "HETATM")):
            atom = int(line[6:11])
            if atom in index:
                for _, _, children in stack:
                    children.append(index[atom])
        elif line.startswith("ENDBRANCH"):
            i, j, children = stack.pop()
            branches.append((i, j, tuple(sorted(children))))
    if stack:
        raise ValueError("Unterminated branches")
    return sorted(branches)


def charged(base, work, q):
    d = work / q
    d.mkdir(exist_ok=True)
    native = (base / q / "ligand.pdbqt").read_text()
    m = json.loads((base / q / "preparation.json").read_text())
    mol = Chem.MolFromSmiles(m["canonical_smiles"])
    original = atoms(native)
    source = m["serial"]
    xyz = np.asarray([original[i][0] for i in source])
    conf = Chem.Conformer(len(source))
    for i, pos in enumerate(xyz):
        conf.SetAtomPosition(i, pos)
    mol.AddConformer(conf)
    amines = [
        a
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == 7
        and a.GetDegree() == 3
        and not a.GetIsAromatic()
        and not a.GetFormalCharge()
    ]
    if len(amines) != 1:
        raise ValueError("Expected one tertiary amine")
    n = amines[0].GetIdx()
    amines[0].SetFormalCharge(1)
    amines[0].SetNumExplicitHs(1)
    Chem.SanitizeMol(mol)
    withh = Chem.AddHs(mol, addCoords=True)

    # Preserve existing bound polar hydrogen vectors in the blind prepared conformer.
    def hydrogen_owner(text):
        sm = extract_map(text, "REMARK SMILES IDX ")
        return {sm[k]: v for k, v in extract_map(text, "REMARK H PARENT ").items()}

    orig_h = hydrogen_owner(native)
    for i, serial in enumerate(source):
        if serial not in orig_h:
            continue
        hi = [
            h.GetIdx()
            for h in withh.GetAtomWithIdx(i).GetNeighbors()
            if h.GetAtomicNum() == 1
        ]
        if len(hi) != 1:
            raise ValueError("Ambiguous donor hydrogen")
        withh.GetConformer().SetAtomPosition(hi[0], original[orig_h[serial]][0])
    (setup,) = MoleculePreparation().prepare(withh)
    out, ok, e = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(e)
    dest = atoms(out)
    ser = extract_map(out, "REMARK INDEX MAP ")
    ids = [ser[i + 1] for i in range(len(source))]
    new_xyz = np.asarray([dest[s][0] for s in ids])
    shift = float(np.linalg.norm(new_xyz - xyz, axis=1).max())
    assert shift < 1e-8 and original[source[n]][1] == "NA" and dest[ids[n]][1] == "N"
    assert len(hydrogen_owner(out)) == len(orig_h) + 1
    oldcg = [i for i, k in enumerate(source) if original[k][1].startswith("CG")]
    newcg = [i for i, k in enumerate(ids) if dest[k][1].startswith("CG")]
    if oldcg != newcg:
        raise ValueError("Ring closure changed")
    if branch_layout(native, source) != branch_layout(out, ids):
        raise ValueError("Heavy atom rotor layout changed")
    changed_types = [
        i
        for i, (si, sj) in enumerate(zip(source, ids, strict=True))
        if original[si][1] != dest[sj][1]
    ]
    if changed_types != [n]:
        raise ValueError("Unexpected changes in physical heavy atom types")
    (d / "ligand.pdbqt").write_text(out)
    (d / "receptor.pdbqt").write_bytes((base / q / "receptor.pdbqt").read_bytes())
    m["serial"] = ids
    (d / "preparation.json").write_text(json.dumps(m) + "\n")
    return m, dict(
        query=q,
        canonical_smiles=Chem.MolToSmiles(mol),
        n=n + 1,
        neutral_type=original[source[n]][1],
        charged_type=dest[ids[n]][1],
        physical_heavy_atoms=len(source),
        initial_max_shift_A=shift,
        heavy_rotor_layout_equal=True,
        heavy_ad_type_changes=len(changed_types),
        closure_atom_indices=[i + 1 for i in oldcg],
        neutral_torsdof=native.split("TORSDOF ")[1].split()[0],
        charged_torsdof=out.split("TORSDOF ")[1].split()[0],
    )


def smiles_order(text):
    m = extract_map(text, "REMARK SMILES IDX ")
    owners = extract_map(text, "REMARK H PARENT ")
    return {j: m[k] for k, j in owners.items()}


def write_xyz(l, xyz):
    return l[:30] + "".join(f"{v:8.3f}" for v in xyz) + l[54:]


def transplant(block, donor, acceptor, dmeta, ameta, same):
    # A native same-state block already has all physical, closure and H coordinates.
    if same:
        return block
    dat = atoms(block)
    dtemp = atoms(donor)
    atemp = atoms(acceptor)
    ds, as_ = dmeta["serial"], ameta["serial"]
    common = dict(zip(as_, ds, strict=True))
    # G pseudo-atoms follow the break carbon whose input coordinate they coincide with.
    for k, (p, typ) in atemp.items():
        if re.fullmatch(r"G[0-3]", typ):
            (a,) = [
                a
                for a in as_
                if atemp[a][1] == "C" + typ and np.linalg.norm(p - atemp[a][0]) < 0.002
            ]
            (g,) = [
                j
                for j, (v, t) in dtemp.items()
                if t == typ and np.linalg.norm(v - dtemp[common[a]][0]) < 0.002
            ]
            common[k] = g
    do, ao = smiles_order(donor), smiles_order(acceptor)
    ref = {common_amine: h for h, common_amine in do.items()}
    newpos = {}
    mol = Chem.MolFromSmiles(dmeta["canonical_smiles"])
    for hi, pa in ao.items():
        heavy = common[pa]
        if heavy in ref:
            common[hi] = ref[heavy]
            continue
        atom = mol.GetAtomWithIdx(ds.index(heavy))
        if atom.GetAtomicNum() != 7 or atom.GetDegree() != 3:
            raise ValueError("Unmatched hydrogen is not on the tertiary amine")
        bondvec = np.asarray(
            [dat[ds[k.GetIdx()]][0] - dat[heavy][0] for k in atom.GetNeighbors()]
        )
        v = -(bondvec / np.linalg.norm(bondvec, axis=1)[:, None]).sum(axis=0)
        newpos[hi] = dat[heavy][0] + v / np.linalg.norm(v) * 1.02
    assert set(atemp) == set(common) | set(newpos)
    rows = []
    for l in acceptor.splitlines():
        if l.startswith(("ATOM  ", "HETATM")):
            k = int(l[6:11])
            pos = newpos[k] if k in newpos else dat[common[k]][0]
            l = write_xyz(l, pos)
        rows.append(l)
    return "\n".join(rows) + "\n"


def write(path, items):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=items[0], lineterminator="\n")
        w.writeheader()
        w.writerows(items)


def geometry(meta, block):
    parsed = atoms(block)
    serials = meta["serial"]
    mol = Chem.MolFromSmiles(meta["canonical_smiles"])
    conf = Chem.Conformer(len(serials))
    for i, serial in enumerate(serials):
        conf.SetAtomPosition(i, parsed[serial][0])
    observed = Chem.Mol(mol)
    observed.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(observed)
    assigned = dict(Chem.FindMolChiralCenters(observed))
    expected = Chem.FindMolChiralCenters(mol)
    cg = [s for s in serials if parsed[s][1] == "CG0"]
    closure = ""
    if cg:
        first, second = cg
        assert (
            mol.GetBondBetweenAtoms(serials.index(first), serials.index(second))
            is not None
        )
        closure = np.linalg.norm(parsed[first][0] - parsed[second][0])
    return dict(
        closure_A=closure,
        specified_stereocenters=len(expected),
        stereo_mismatches=sum(assigned.get(k) != v for k, v in expected),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("base", type=Path, help="restored novel input/output directory")
    p.add_argument("work", type=Path)
    p.add_argument("--binary", required=True, type=Path)
    p.add_argument("--output", type=Path, default=HERE / "results")
    p.add_argument("--replication", action="store_true")
    a = p.parse_args()
    base = a.base.resolve()
    work = a.work.resolve()
    binary = a.binary.resolve()
    a.output.mkdir(exist_ok=True, parents=True)
    (work / "runs").mkdir(parents=True, exist_ok=True)
    charged_meta = {}
    diagnostics = []
    for q in CODES:
        m, diag = charged(base, work, q)
        diagnostics.append(diag)
        charged_meta[q] = m
    write(a.output / "preparation.csv", diagnostics)
    neutral = base
    if a.replication:
        neutral = work / "neutral"
        (neutral / "runs").mkdir(parents=True, exist_ok=True)
        for q in CODES:
            (neutral / q).mkdir(exist_ok=True)
            for name in ["receptor.pdbqt", "ligand.pdbqt", "preparation.json"]:
                (neutral / q / name).write_bytes((base / q / name).read_bytes())
    states = {"neutral": neutral, "cation": work}
    metadata = {
        "cation": charged_meta,
        "neutral": {
            q: json.loads((base / q / "preparation.json").read_text()) for q in CODES
        },
    }
    cells = []
    results = []
    scores = []
    for q in CODES:
        for r in CODES:
            for s in (43103, 43117) if a.replication else SEEDS:
                pools = {}
                for arm, d in states.items():
                    row, modes = paired.execute(
                        d,
                        metadata[arm],
                        q,
                        r,
                        s,
                        12000,
                        "baseline" if arm == "neutral" else "cation",
                        binary,
                        arm == "neutral" and not a.replication,
                    )
                    pools[arm] = [
                        (b, z)
                        for b, z in zip(
                            (d / "runs" / (row["key"] + ".pdbqt"))
                            .read_text()
                            .split("ENDMDL")[:-1],
                            modes,
                            strict=True,
                        )
                    ]
                    cells.append(dict(state=arm, **row))
                for source, pool in pools.items():
                    for b, z in pool:
                        record = dict(
                            query=q,
                            receiver=r,
                            seed=s,
                            source=source,
                            rank=z["mode"],
                            rmsd_A=z["rmsd_A"],
                            saved_score=z["score"],
                            clashes=z["clash_pairs"],
                            **geometry(metadata[source][q], b),
                        )
                        for state in states:
                            path = work / "projected.pdbqt"
                            text = transplant(
                                b,
                                (states[source] / q / "ligand.pdbqt").read_text(),
                                (states[state] / q / "ligand.pdbqt").read_text(),
                                metadata[source][q],
                                metadata[state][q],
                                state == source,
                            )
                            if state == source:
                                text = "\n".join(
                                    l
                                    for l in text.splitlines()
                                    if not l.startswith("MODEL")
                                )
                            path.write_text(text)
                            start = time.perf_counter()
                            cmd = [
                                str(binary),
                                "--cpu",
                                "1",
                                "--seed",
                                "71993",
                                "--score_only",
                                "--receptor",
                                str(base / r / "receptor.pdbqt"),
                                "--ligand",
                                str(states[state] / q / "ligand.pdbqt"),
                                "--pose_coordinates",
                                str(path),
                            ]
                            out = subprocess.run(cmd, capture_output=True, text=True)
                            if out.returncode:
                                raise ValueError(out.stdout + out.stderr)
                            en = re.findall(
                                r":\s+(-?\d+\.\d+)\s+\(kcal/mol\)",
                                out.stdout.split("Estimated Free Energy of Binding")[
                                    -1
                                ],
                            )
                            if len(en) != 10:
                                raise ValueError(out.stdout)
                            record[state + "_inter"] = float(en[1])
                            record[state + "_intra"] = float(en[4])
                            record[state + "_energy"] = float(en[1]) + float(en[4])
                            record[state + "_seconds"] = time.perf_counter() - start
                        scores.append(record)
                current = [
                    z
                    for z in scores
                    if (z["query"], z["receiver"], z["seed"]) == (q, r, s)
                ]
                for source in ["neutral", "cation", "union"]:
                    allowed = [
                        x for x in current if source == "union" or x["source"] == source
                    ]
                    for state in states:
                        ranked = sorted(
                            allowed,
                            key=lambda x: (
                                x[state + "_energy"],
                                x["source"] != "neutral",
                                x["rank"],
                            ),
                        )
                        first = ranked[0]
                        results.append(
                            dict(
                                query=q,
                                receiver=r,
                                seed=s,
                                pool=source,
                                scoring=state,
                                selected_source=first["source"],
                                rank=first["rank"],
                                top1_A=first["rmsd_A"],
                                best5_A=min(x["rmsd_A"] for x in ranked[:5]),
                                oracle_A=min(x["rmsd_A"] for x in ranked),
                                n=len(allowed),
                            )
                        )
    write(a.output / "cells.csv", cells)
    write(a.output / "fixed_modes.csv", scores)
    write(a.output / "fixed_pools.csv", results)


if __name__ == "__main__":
    main()
