"""Conditional O-H torsion tests in two published neutron complexes."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys

import numpy as np
from scipy.special import softmax

from bonded import proper_energy, unit
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "hydroxyl_partition"))
from scorer import RT, quality

# The two named sites and partners come from the experimental articles. Distances
# are recomputed. Residue 106 also has a second carboxylate oxygen, outside the ramp.
CASES = {
    "4G0C": dict(res="THR", site=199, carbon="CB", oxygen="OG1", h="DG1", frame="CA",
                 ligand=(302, "N1"), ligand_h=(302, "DN11"), internal=[(106,"OE1"),(106,"OE2")],
                 checksum="6fed0b6c7df1228e32a2b14f01074c71bd8181720e3d91c3096ed462bc422711"),
    "2R24": dict(res="TYR", site=48, carbon="CZ", oxygen="OH", h="DH", frame="CE1",
                 ligand=(317, "O33"), ligand_h=None, internal=[],
                 checksum="7168bfc1f52af846d663e70b9f8554d1bf77a1dec82181dd86320da0eaecf0cb"),
}


def angle(vertex, a, b):
    return float(np.degrees(np.arccos(np.clip(unit(a-vertex) @ unit(b-vertex), -1., 1.))))


def evaluate(pdb, src, spec):
    block = (src / (pdb + ".pdb")).read_bytes()
    if hashlib.sha256(block).hexdigest() != spec["checksum"]:
        raise ValueError("Full wwPDB file differs from audited revision: " + pdb)
    atoms = {}
    for line in block.decode().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[21] == "A" and line[16] in " A":
            atoms[int(line[22:26]), line[12:16].strip()] = (
                np.array([float(line[q:q+8]) for q in (30, 38, 46)]),
                line[76:78].strip(), float(line[54:60]))
    coords = {name: value[0] for (idx, name), value in atoms.items() if idx == spec["site"]}
    c, o = coords[spec["carbon"]], coords[spec["oxygen"]]
    axis = unit(o - c)
    anchor = coords[spec["frame"]] - c
    u = unit(anchor - np.dot(anchor, axis) * axis)
    v = np.cross(axis, u)
    phi = np.radians(np.arange(0, 360, 5))
    nodes = o + .97 * (axis / 3 + np.sqrt(8/9) * (np.cos(phi)[:, None]*u + np.sin(phi)[:, None]*v))
    # Only the reporting steps below use the observed side-chain D.
    observed = coords[spec["h"]]
    azimuth = np.arctan2((observed-o) @ v, (observed-o) @ u)
    err = np.abs(np.degrees(np.angle(np.exp(1j * (phi-azimuth)))))
    ligand = atoms[spec["ligand"]][0]
    donor = spec["ligand_h"]
    reverse = 0. if donor is None else float(quality(ligand, atoms[donor][0], o))
    contacts = []
    internal = np.zeros(len(nodes))
    for partner in [spec["ligand"], *spec["internal"]]:
        a, element, _ = atoms[partner]
        r = np.linalg.norm(a-o)
        magnitude = .587439 * np.clip((1.7 + {"N":1.8, "O":1.7}[element] - r) / .7, 0., 1.)
        gate = quality(o, nodes, a)
        contacts.append(dict(pdb=pdb, residue=spec["res"]+str(spec["site"]), partner=str(partner[0])+":"+partner[1],
                             occupancy=atoms[spec["site"],spec["h"]][2], heavy_A=r, D_to_partner_A=np.linalg.norm(a-observed),
                             O_D_A_deg=angle(observed,o,a), observed_forward_gate=quality(o,observed,a),
                             observed_reverse_gate=reverse if partner==spec["ligand"] else 0., radial_magnitude=magnitude))
        if partner == spec["ligand"]:
            ligand_only = -magnitude * gate
            dual_ligand = -magnitude * np.maximum(gate, reverse)
        else:
            internal -= magnitude * gate
    torsion = proper_energy(spec["res"], coords, nodes)
    models = dict(flat=np.zeros(len(nodes)), protein_only=internal, force_receptor_donation=ligand_only,
                  ligand_with_reverse=dual_ligand, joint=internal+dual_ligand, bonded_only=torsion,
                  joint_bonded=internal+dual_ligand+torsion)
    rows, scans = [], []
    for name, energy in models.items():
        prob = softmax(-energy/RT)
        minima = np.isclose(energy, energy.min(), atol=1.e-9, rtol=0)
        rows.append(dict(pdb=pdb, model=name, nodes=len(nodes), observed_azimuth_deg=np.degrees(azimuth),
                         span_kcal_mol=np.ptp(energy), min_error_deg=(min(err[minima]) if not minima.all() else ""),
                         max_minimizer_error_deg=(max(err[minima]) if not minima.all() else ""),
                         posterior_within_30deg=sum(prob[err<=30]), flat_fraction=np.mean(err<=30)))
        scans.extend(dict(pdb=pdb, model=name, phi_deg=round(np.degrees(p)), energy=val, weight=pr)
                     for p,val,pr in zip(phi,energy,prob))
    return contacts, rows, scans


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("pdb_dir", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    outs = [[] for _ in range(3)]
    for pdb, case in CASES.items():
        for dest, rows in zip(outs, evaluate(pdb, args.pdb_dir, case)):
            dest.extend(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, rows in zip(("neutron_contacts", "neutron_modes", "neutron_scan"), outs):
        with (args.output / (name+".csv")).open("w", newline="") as stream:
            w = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
