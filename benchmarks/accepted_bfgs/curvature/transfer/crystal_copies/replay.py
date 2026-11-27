"""Measure VHL fixed-frame RMSD to every complete asymmetric-unit ligand copy."""

import argparse
import csv
import itertools
import json
import tarfile
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation

PARENT = Path(__file__).resolve().parents[1]
PDBS = ("4W9H", "5LLI")
CHAINS = "CFIL"
ARMS = ("official", "parent", "fixed", "parent_double")


def read_text(archive, name):
    return archive.extractfile(name).read().decode()


def xyz(line):
    return [float(line[k : k + 8]) for k in (30, 38, 46)]


def write(path, rows):
    with path.open("w", newline="") as f:
        out = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def read_structure(archive, pdb):
    lines = read_text(archive, f"raw/{pdb}.pdb").splitlines()
    meta = json.loads(read_text(archive, f"{pdb}/metadata.json"))
    template = Chem.MolFromSmiles(meta["smiles"])
    crystals, ca = {}, {}
    for ch in CHAINS:
        lig = [
            line
            for line in lines
            if line.startswith("HETATM")
            and line[21:27] == f"{ch} 301 "
            and line[17:20] == meta["component"]
        ]
        assert len(lig) == template.GetNumAtoms()
        assert all(l[16] == " " and float(l[54:60]) == 1.0 for l in lig)
        atom_ids = {int(l[6:11]) for l in lig}
        bonds = [
            l for l in lines if l.startswith("CONECT") and int(l[6:11]) in atom_ids
        ]
        raw = Chem.MolFromPDBBlock(
            "\n".join(lig + bonds + ["END", ""]),
            sanitize=False,
            removeHs=False,
            proximityBonding=False,
        )
        ref = AllChem.AssignBondOrdersFromTemplate(template, raw)
        Chem.AssignStereochemistryFrom3D(ref, replaceExistingTags=True)
        assert Chem.MolToInchiKey(ref) == meta["inchikey"]
        matches = ref.GetSubstructMatches(
            template, useChirality=True, uniquify=False, maxMatches=10001
        )
        assert 0 < len(matches) < 10001 and len(matches) == len(meta["crystal"])
        crystals[ch] = ref.GetConformer().GetPositions()[np.array(matches)]
        atoms = {}
        for l in lines:
            if (
                l.startswith(("ATOM  ", "HETATM"))
                and l[21] == ch
                and l[12:16].strip() == "CA"
                and l[16] in " A"
                and (l.startswith("ATOM") or l[17:20] == "CAS")
            ):
                # The historical registration labels the deposited modified Cys as CYS.
                name = "CYS" if l[17:20] == "CAS" else l[17:20]
                key = (l[22:27], name)
                if key in atoms:
                    raise ValueError(f"ambiguous C-alpha {pdb} {ch} {key}")
                atoms[key] = xyz(l)
        ca[ch] = atoms
    pocket = {(key[1:6], key[6:]) for key in meta["ca"] if key[0] == "C"}
    return meta["serial"], crystals, ca, pocket


def transform(a, b, x):
    origin, dest = a.mean(0), b.mean(0)
    rot, _ = Rotation.align_vectors(b - dest, a - origin)
    return rot.apply((x - origin).reshape(-1, 3)).reshape(x.shape) + dest


def register(data, q, r, kind):
    _, crystals, qca, qpocket = data[q]
    _, _, rca, rpocket = data[r]
    # Use identical labelled positions for all four query copies within the method.
    common = set(rca["C"]).intersection(*(set(qca[c]) for c in CHAINS))
    keys = sorted(common if kind == "chain" else common & qpocket & rpocket)
    if kind == "pocket":
        assert len(keys) == len(qpocket & rpocket) == 55
    refs, fits = {}, []
    b = np.array([rca["C"][k] for k in keys])
    for ch in CHAINS:
        a = np.array([qca[ch][k] for k in keys])
        refs[ch] = transform(a, b, crystals[ch])
        fits.append(
            {
                "query": q,
                "receiver": r,
                "frame": kind,
                "copy": ch,
                "ca_count": len(keys),
                "ca_rmsd_A": rmsd(transform(a, b, a), b),
                "heavy_atoms": crystals[ch].shape[1],
                "stereo_maps": len(crystals[ch]),
            }
        )
    return refs, fits


def rmsd(x, y):
    # Fixed receptor frame; an optional leading dimension enumerates graph maps.
    return float(np.sqrt(np.square(x - y).sum(axis=-1).mean(axis=-1)).min())


def poses(text, serial):
    for block in text.split("ENDMDL"):
        if block.strip():
            atoms = {
                int(l[6:11]): xyz(l)
                for l in block.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            yield np.array([atoms[k] for k in serial])


def replay(out):
    fits, separations, values, checks, deltas = [], [], [], [], []
    with (PARENT / "results/modes.csv").open() as f:
        archived = {
            (v["query"], v["arm"], int(v["seed"]), int(v["rank"])): float(v["rmsd"])
            for v in csv.DictReader(f)
        }
    with (
        tarfile.open(PARENT / "evidence.tar.xz") as inputs,
        tarfile.open(PARENT / "witness.tar.xz") as wit,
    ):
        data = {p: read_structure(inputs, p) for p in PDBS}
        for q, r in itertools.permutations(PDBS):
            for kind in ("pocket", "chain"):
                refs, rows = register(data, q, r, kind)
                fits.extend(rows)
                for c, d in itertools.combinations(CHAINS, 2):
                    separations.append(
                        {
                            "query": q,
                            "receiver": r,
                            "frame": kind,
                            "copy1": c,
                            "copy2": d,
                            "rmsd_A": rmsd(refs[c][0], refs[d]),
                        }
                    )
                first = {}
                for arm, seed in itertools.product(ARMS, (27183, 91439)):
                    name = f"runs/{q}_{r}_{seed}_{arm}.pdbqt"
                    for rank, pt in enumerate(
                        poses(read_text(inputs, name), data[q][0]), 1
                    ):
                        for ch in CHAINS:
                            err = rmsd(pt, refs[ch])
                            values.append(
                                {
                                    "query": q,
                                    "receiver": r,
                                    "frame": kind,
                                    "copy": ch,
                                    "arm": arm,
                                    "seed": seed,
                                    "rank": rank,
                                    "rmsd_A": err,
                                }
                            )
                            if kind == "pocket" and ch == "C":
                                checks.append(abs(err - archived[q, arm, seed, rank]))
                            if rank == 1:
                                first[arm, seed, ch] = err
                (local,) = poses(
                    read_text(wit, f"witness/{q}_fixed_local.pdbqt"), data[q][0]
                )
                for ch in CHAINS:
                    values.append(
                        {
                            "query": q,
                            "receiver": r,
                            "frame": kind,
                            "copy": ch,
                            "arm": "guided_local",
                            "seed": "",
                            "rank": 1,
                            "rmsd_A": rmsd(local, refs[ch]),
                        }
                    )
                    for seed in (27183, 91439):
                        old, new = first["parent", seed, ch], first["fixed", seed, ch]
                        deltas.append(
                            {
                                "query": q,
                                "receiver": r,
                                "frame": kind,
                                "copy": ch,
                                "seed": seed,
                                "parent_A": old,
                                "fixed_A": new,
                                "fixed_minus_parent_A": new - old,
                            }
                        )
    assert len(checks) == len(archived) == 91 and max(checks) < 1e-8
    for name, rows in (
        ("registrations", fits),
        ("crystal_pairs", separations),
        ("poses", values),
        ("paired_top1", deltas),
    ):
        for row in rows:
            for key in row:
                if isinstance(row[key], float):
                    row[key] = f"{row[key]:.8f}"
        write(out / f"{name}.csv", rows)
    print(
        f"{len(checks)} archived RMSDs reproduced; max |difference| {max(checks):.3g} Å"
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    RDLogger.DisableLog(
        "rdApp.warning"
    )  # Template symmetry can have equivalent mappings.
    replay(args.out)
