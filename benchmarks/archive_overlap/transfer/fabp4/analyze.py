"""Fixed-coordinate audit of the two 3P6E archive observations and both deposited conformers."""

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_transfer import independent_reference, xyz
from run_panel import paired, sha


def local_frame(positions, acid_ids):
    carbon, carbonyl, hydroxyl = (positions[j] for j in acid_ids)
    ex = carbon - hydroxyl
    ex = ex / np.linalg.norm(ex)
    ey = carbonyl - hydroxyl
    ey = ey - np.dot(ey, ex) * ex
    ey = ey / np.linalg.norm(ey)
    return np.array([ex, ey, np.cross(ex, ey)])


def main(args):
    work, result, out = args.inputs, args.results, args.out
    out.mkdir(parents=True, exist_ok=True)
    q = work / "3P6E"
    template = (q / "ligand.pdbqt").read_text().splitlines()
    smiles = next(s.split()[2] for s in template if s.startswith("REMARK SMILES "))
    mol = Chem.MolFromSmiles(smiles)
    refs, _, center, serial, alt, _info = independent_reference(work, "3P6E")
    acid = mol.GetSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"))
    methoxy = mol.GetSubstructMatch(Chem.MolFromSmarts("[CH3]-[O]-[c]"))
    assert len(acid) == 3 and len(methoxy) == 3
    acid_ids = [serial[j] for j in acid]
    ether = serial[methoxy[1]]
    hydrogen = int(
        next(s.split()[-1] for s in template if s.startswith("REMARK H PARENT"))
    )
    source = {
        int(l[6:11]): xyz(l) for l in template if l.startswith(("ATOM", "HETATM"))
    }
    old_oh = source[hydrogen] - source[acid_ids[2]]
    internal_h = local_frame(source, acid_ids) @ old_oh
    structures = {}
    for label in ("A", "B"):
        coordinates = {k: v for k, v in zip(serial, refs[alt.index(label)])}
        coordinates[hydrogen] = coordinates[acid_ids[2]] + internal_h @ local_frame(
            coordinates, acid_ids
        )
        structures["crystal_" + label] = [
            line[:30]
            + "".join(f"{v:8.3f}" for v in coordinates[int(line[6:11])])
            + line[54:]
            if line.startswith(("ATOM", "HETATM"))
            else line
            for line in template
        ]
    for seed, correct in [(684173, 9), (961751, 5)]:
        modes = (
            (result / "runs" / f"3P6E_3P6E_{seed}_candidate.pdbqt")
            .read_text()
            .split("ENDMDL")
        )
        for rank in (1, correct):
            structures[f"{seed}_rank{rank}"] = [
                l
                for l in modes[rank - 1].splitlines()
                if l and not l.startswith("MODEL")
            ]
    raw = (work / "raw/3P6E.pdb").read_text().splitlines()
    protein = [
        l
        for l in raw
        if l.startswith("ATOM  ")
        and l[21] == "A"
        and l[16] in " A"
        and l[76:78].strip() not in ("H", "D")
    ]
    waters = [
        l
        for l in raw
        if l.startswith("HETATM")
        and l[17:20] == "HOH"
        and l[21] == "A"
        and l[76:78].strip() == "O"
    ]
    arg = np.array(
        [
            xyz(l)
            for l in protein
            if l[17:20] == "ARG"
            and int(l[22:26]) == 126
            and l[12:16].strip() in ("NE", "NH1", "NH2")
        ]
    )
    tyr = np.array(
        [
            xyz(l)
            for l in protein
            if l[17:20] == "TYR" and int(l[22:26]) == 128 and l[12:16].strip() == "OH"
        ]
    )
    c_atoms = np.array([xyz(l) for l in protein if l[76:78].strip() == "C"])
    wxyz = np.array([xyz(l) for l in waters])
    rows = []
    fixed_unbound = -0.230
    for key, lines in structures.items():
        path = out / (key + ".pdbqt")
        path.write_text("\n".join(lines) + "\n")
        atoms = {
            int(l[6:11]): xyz(l) for l in lines if l.startswith(("ATOM", "HETATM"))
        }
        oxy = np.array([atoms[j] for j in acid_ids[1:]])
        coff = np.array([atoms[ether]])
        lc = np.array(
            [atoms[serial[a.GetIdx()]] for a in mol.GetAtoms() if a.GetSymbol() == "C"]
        )
        heavy = np.array([atoms[j] for j in serial])
        water_d = cdist(oxy, wxyz).min(axis=0)
        cmd = [
            str(args.upstream.resolve()),
            "--receptor",
            str(q / "receptor.pdbqt"),
            "--ligand",
            str(path),
            "--score_only",
            f"--unbound_energy={fixed_unbound}",
            "--cpu",
            "1",
            "--seed",
            "1",
        ]
        for axis, value in zip("xyz", center):
            cmd.extend(["--center_" + axis, str(value), "--size_" + axis, "22.5"])
        log = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
        (out / (key + ".log")).write_text(log)

        def field(label):
            return float(
                re.search(re.escape(label) + r"\s*:\s*([\d.\-]+)", log).group(1)
            )

        rows.append(
            dict(
                pose=key,
                rmsd_A=round(
                    min(np.linalg.norm(heavy - p) / np.sqrt(len(heavy)) for p in refs),
                    4,
                ),
                score=field("Estimated Free Energy of Binding"),
                inter=field("Final Intermolecular Energy"),
                intra=field("Final Total Internal Energy"),
                torsion=field("Torsional Free Energy"),
                unbound=fixed_unbound,
                acid_to_Arg126_N_A=round(cdist(oxy, arg).min(), 3),
                acid_to_Tyr128_OH_A=round(cdist(oxy, tyr).min(), 3),
                ether_to_Arg126_N_A=round(cdist(coff, arg).min(), 3),
                ether_to_Tyr128_OH_A=round(cdist(coff, tyr).min(), 3),
                acid_to_water_A=round(water_d.min(), 3),
                nearest_acid_water=waters[water_d.argmin()][22:26].strip(),
                water_occupancy=float(waters[water_d.argmin()][54:60]),
                ligand_to_any_water_A=round(cdist(heavy, wxyz).min(), 3),
                ligand_protein_CC_pairs_within_4p5A=int(
                    (cdist(lc, c_atoms) <= 4.5).sum()
                ),
            )
        )
    paired.write(out / "fixed.csv", rows)
    (out / "native_sha256").write_text(sha(args.upstream) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("inputs", "results", "upstream", "out"):
        p.add_argument("--" + k, type=Path, required=True)
    main(p.parse_args())
