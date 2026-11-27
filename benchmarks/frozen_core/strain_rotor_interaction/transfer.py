"""Recompute the ligand reference and first-pose components on Sörme Gal3 and Davis S6K1."""
import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from tempfile import TemporaryDirectory

from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path[:0] = [str(BENCH / "terminal_rotor_score"), str(HERE.parent / "affinity"),
                str(HERE.parent / "strain")]
import assay
import test_entropy as entropy
from pose import physical, score_pose


def score_compound(task):
    group, history = task
    rows = []
    trees = {}
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "first.pdbqt"
        prepared = Path(tmp) / "ligand.pdbqt"
        for target, ligand, measured, rec, seed, receptor, ligand_bytes, pool in group:
            raw = pool.decode().split("ENDMDL", 1)[0] + "ENDMDL\n"
            path.write_text(raw)
            prepared.write_bytes(ligand_bytes)
            pose = physical.read_poses(path)[0]
            physical.fixed_coordinate_inputs([pose], prepared, None)
            molecule, smiles = physical.decode_ligand(pose)
            if not rows:
                identity = smiles
                free, count = physical.gas_reference(molecule)
            if rec not in trees:
                atoms = map(entropy.atom, entropy.atom_lines(receptor.decode()))
                trees[rec] = cKDTree([a[2] for a in atoms if entropy.heavy(a)])
            saved = history[target, ligand, rec, seed]
            tors = float(saved["old_rotors"])
            scale = 1 / (1 + assay.WEIGHT * tors)
            energies = score_pose(pose, identity, free, count, scale)
            # Production contact-rotor scoring has no model for ring-closure or hydration pseudo-atoms.
            types = [a[1] for a in map(entropy.atom, entropy.atom_lines(ligand_bytes.decode()))]
            pseudo = any(t == "W" or t.startswith(("G", "CG")) for t in types)
            if pseudo:
                rotors = dict(native=pose.native_score, raw=float(saved["saved_inter"]), contact="", constrained="")
            else:
                rotors = entropy.pose(raw, trees[rec], entropy.ligand_graph(ligand_bytes.decode()), tors)
            if rotors["native"] != float(saved["saved_score"]):
                raise ValueError(f"Historical native score differs for {target} {ligand}")
            rows.append(dict(target=target, compound=ligand, receptor=rec, seed=seed,
                             measurement=measured, pool_sha256=hashlib.sha256(pool).hexdigest(),
                             native=rotors["native"], raw=rotors["raw"], unit=energies["full_unit"],
                             scale=scale, contact=rotors["contact"], tors=tors,
                             active_heavy_tors=float(saved["fixed_rotors"]), constrained=rotors["constrained"],
                             heavy=molecule.GetNumHeavyAtoms(), charge=physical.Chem.GetFormalCharge(molecule),
                             fixed_mmff=energies["fixed_mmff"], free_mmff=free, nfree=count,
                             rotor_status="unsupported_pseudo_atoms" if pseudo else "ok"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    history = {(r["target"], r["ligand"], r["receptor"], int(r["seed"])): r
               for r in assay.read(BENCH / "terminal_rotor_score/results/poses.csv")}
    groups = defaultdict(list)
    for entry in assay.assay_inputs():
        if entry[0] in {"Gal3", "S6K1"}:
            groups[entry[:2]].append(entry)
    with ProcessPoolExecutor(args.jobs) as pool:
        rows = [r for batch in pool.map(score_compound, [(v, history) for v in groups.values()]) for r in batch]
    if len(rows) != 122 or len(groups) != 61:
        raise ValueError("Transfer panel incomplete")
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
