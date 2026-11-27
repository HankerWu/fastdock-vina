"""Rescore six frozen FABP4 poses after changing only the ligand acid microstate."""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import tarfile
from pathlib import Path

import meeko
import rdkit
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

HERE = Path(__file__).resolve().parent
REFERENCE = "684173_rank1"
FIELDS = (
    "Estimated Free Energy of Binding",
    "Final Intermolecular Energy",
    "Final Total Internal Energy",
    "Torsional Free Energy",
    "Unbound System's Energy",
)


def atom_rows(pdbqt):
    return {
        int(s[6:11]): s for s in pdbqt.splitlines() if s.startswith(("ATOM", "HETATM"))
    }


def xyz(s):
    return tuple(float(s[k : k + 8]) for k in (30, 38, 46))


def index_map(pdbqt, prefix):
    mapping = {}
    for line in pdbqt.splitlines():
        if line.startswith(prefix):
            values = list(map(int, line[len(prefix) :].split()))
            for index, serial in zip(values[::2], values[1::2], strict=True):
                if index - 1 in mapping:
                    raise ValueError("Duplicate atom index")
                mapping[index - 1] = serial
    return mapping


def prepare(source, state, meta):
    smiles = next(
        x.split()[2] for x in source.splitlines() if x.startswith("REMARK SMILES ")
    )
    assert smiles == meta["canonical_smiles"]
    mol = Chem.MolFromSmiles(smiles)
    old_serial = index_map(source, "REMARK SMILES IDX ")
    assert [old_serial[i] for i in range(mol.GetNumAtoms())] == meta["serial"]
    old = atom_rows(source)
    acid = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"))
    assert len(acid) == 1
    hydroxyl = mol.GetAtomWithIdx(acid[0][2])
    if state == "carboxylate":
        hydroxyl.SetFormalCharge(-1)
        hydroxyl.SetNumExplicitHs(0)
        hydroxyl.SetNoImplicit(True)
        Chem.SanitizeMol(mol)
    n_heavy = mol.GetNumAtoms()
    coords = Chem.Conformer(n_heavy)
    for i in range(n_heavy):
        atom = mol.GetAtomWithIdx(i)
        record = old[old_serial[i]]
        assert (record.split()[-1] in ("C", "A") and atom.GetSymbol() == "C") or (
            record.split()[-1] == "OA" and atom.GetSymbol() == "O"
        )
        coords.SetAtomPosition(i, xyz(record))
    mol.AddConformer(coords)
    mol = Chem.AddHs(mol, addCoords=True)
    if state == "neutral":
        from_h = index_map(source, "REMARK H PARENT ")
        assert list(from_h) == [acid[0][2]]
        (h,) = [
            a.GetIdx()
            for a in mol.GetAtomWithIdx(acid[0][2]).GetNeighbors()
            if a.GetAtomicNum() == 1
        ]
        mol.GetConformer().SetAtomPosition(h, xyz(old[from_h[acid[0][2]]]))
    (setup,) = MoleculePreparation().prepare(mol)
    text, ok, error = PDBQTWriterLegacy.write_string(setup, add_index_map=True)
    if not ok:
        raise ValueError(error)
    new = atom_rows(text)
    serial = index_map(text, "REMARK INDEX MAP ")
    assert len({serial[i] for i in range(n_heavy)}) == n_heavy == 13
    for i in range(n_heavy):
        assert xyz(new[serial[i]]) == xyz(old[old_serial[i]])
        assert old[old_serial[i]].split()[-1] == new[serial[i]].split()[-1]
    charge = sum(float(r[70:76]) for r in new.values())
    formal_charge = Chem.GetFormalCharge(mol)
    assert abs(charge - formal_charge) < 0.005 and formal_charge == (
        0 if state == "neutral" else -1
    )
    assert sum(r.split()[-1] == "HD" for r in new.values()) == int(state == "neutral")
    return text, {
        "state": state,
        "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
        "formal_charge": formal_charge,
        "pdbqt_charge_sum": round(charge, 3),
        "heavy": n_heavy,
        "max_heavy_displacement_A": 0,
        "torsdof": int(
            next(l.split()[1] for l in text.splitlines() if l.startswith("TORSDOF"))
        ),
    }


def score(binary, receptor, center, pose, unbound, log):
    args = [
        binary,
        "--receptor",
        receptor,
        "--ligand",
        pose,
        "--score_only",
        f"--unbound_energy={unbound}",
        "--seed",
        "1",
        "--cpu",
        "1",
        "--scoring",
        "vina",
    ]
    for axis, value in zip("xyz", center, strict=True):
        args.extend(["--center_" + axis, str(value), "--size_" + axis, "22.5"])
    result = subprocess.run(args, check=True, capture_output=True, text=True)
    log.write_text(result.stdout + result.stderr)
    return result.stdout


def values(log):
    return [
        float(re.search(re.escape(field) + r"\s*:\s*([\d.\-]+)", log)[1])
        for field in FIELDS
    ]


def run(binary, out):
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(HERE.parents[1] / "inputs.tar.xz") as inputs:
        receptor = out / "receptor.pdbqt"
        receptor.write_bytes(inputs.extractfile("3P6E/receptor.pdbqt").read())
        meta = json.load(inputs.extractfile("3P6E/preparation.json"))
    with tarfile.open(HERE.parent / "native_scores.tar.xz") as archive:
        originals = {
            Path(m.name).stem: archive.extractfile(m).read().decode()
            for m in archive
            if m.name.endswith(".pdbqt")
        }
        original_logs = {
            key: archive.extractfile(key + ".log").read().decode() for key in originals
        }
    assert len(originals) == 6
    rows = []
    for state in ("neutral", "carboxylate"):
        preparations = {}
        for key, source in originals.items():
            pdbqt, info = prepare(source, state, meta)
            path = out / (key + "-" + state + ".pdbqt")
            path.write_text(pdbqt)
            preparations[key] = (path, info)
        # The reference is the same physical ligand geometry in both states.
        # Use just its intra field, shared among every pose in this state.
        reference = values(
            score(
                binary,
                receptor,
                meta["center"],
                preparations[REFERENCE][0],
                0,
                out / (state + "-reference.log"),
            )
        )[2]
        for key, (path, info) in preparations.items():
            raw = score(
                binary,
                receptor,
                meta["center"],
                path,
                reference,
                path.with_suffix(".log"),
            )
            energy = values(raw)
            zero = out / (path.stem + "-zero.pdbqt")
            zero.write_text(
                "\n".join(
                    line[:70] + " 0.000" + line[76:]
                    if line.startswith(("ATOM", "HETATM"))
                    else line
                    for line in path.read_text().splitlines()
                )
                + "\n"
            )
            qcontrol = score(
                binary,
                receptor,
                meta["center"],
                zero,
                reference,
                zero.with_suffix(".log"),
            )
            terms = lambda s: re.findall(r":\s*([-+.\deE]+)\s+\(kcal/mol\)", s)
            assert len(terms(raw)) == 10 and terms(raw) == terms(qcontrol)
            if state == "neutral":
                assert terms(raw) == terms(original_logs[key]), key
            rows.append(
                {
                    "pose": key,
                    **info,
                    **dict(
                        zip(
                            ("score", "inter", "intra", "torsion", "unbound"),
                            energy,
                            strict=True,
                        )
                    ),
                    "zero_charge_terms_equal": len(terms(raw)),
                }
            )
    with (out / "scores.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (out / "identity.json").write_text(
        json.dumps(
            {
                "rdkit": rdkit.__version__,
                "meeko": meeko.__version__,
                "upstream_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    print((out / "scores.csv").read_text())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("upstream", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    run(args.upstream.resolve(), args.output)
