"""Localize the observed 2PRG rank reversal to moving individual receptor hydroxyl hydrogens."""
from copy import copy
import argparse
import csv
import json
from pathlib import Path
import sys

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "frozen_core/hydroxyl_partition"))
from scorer import Receptor, score_pose


def main(work):
    rec = Receptor((work / "2PRG/receptor.pdbqt").read_text())
    meta = json.loads((work / "2PRG/preparation.json").read_text())
    ligand = (work / "2PRG/ligand.pdbqt").read_text()
    atoms = (work / "runs/2PRG_12000_43117_pro.pdbqt").read_text().split("ENDMDL")
    rows = []
    for rank, pose in enumerate(atoms[:3], 1):
        full = score_pose(rec, ligand, pose)
        changes = []
        for serial, values in rec.rotors.items():
            active = copy(rec)
            active.rotors = {serial: values}
            one = score_pose(active, ligand, pose)
            if one["rotors"] == 0:
                continue
            atom = rec.atoms[serial]
            shift = one["raw_ensemble"] - one["raw_fixed"]
            changes.append(shift)
            rows.append(dict(native_rank=rank, atom=atom.name, resname=atom.residue[:3],
                             deposition_chain=atom.residue[4], deposition_number=meta["receptor_ordinal_to_deposition"][int(atom.residue[5:])-1].strip(),
                             full_raw_fixed=full["raw_fixed"], full_raw_ensemble=full["raw_ensemble"],
                             raw_fixed_to_ensemble_shift=shift, rotor_contacts=one["mobile_contacts"],
                             receptor_contact_prior_range=float(values[1].max()-values[1].min())))
        assert abs(sum(changes) - (full["raw_ensemble"]-full["raw_fixed"])) < 1.e-9
    with (work / "polar_switch.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    main(p.parse_args().work)
