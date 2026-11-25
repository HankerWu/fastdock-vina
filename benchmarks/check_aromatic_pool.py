#!/usr/bin/env python3
"""Exercise fixed-pool aromatic ordering, the score boundary and ring validation."""

from math import cos, pi, sin
from pathlib import Path
from tempfile import TemporaryDirectory
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from score_physical_pool import PHENYL_RING, receptor_ring_planes, score_aromatic_pool


def ring_atoms(names, residue, z, *, ruffle=0):
    rows = []
    for index, name in enumerate(names):
        theta = pi * index / 3
        xyz = (1.4 * cos(theta), 1.4 * sin(theta), z + (ruffle if index == 2 else 0))
        rows.append(f"ATOM  {index + 1:5d} {name:>4} {residue:>3} A   1    "
                    f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00     0.000 A\n")
    return "".join(rows)


CARBONS = [f"C{i}" for i in range(1, 7)]
HYDROGENS = [f"H{i}" for i in range(1, 7)]


def component():
    return ("data_BEN\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n"
            + "".join(f"{atom} {atom[0]} 0\n" for atom in CARBONS + HYDROGENS)
            + "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n"
            + "".join(f"C{i} C{i % 6 + 1} {'DOUB' if i % 2 else 'SING'}\nC{i} H{i} SING\n"
                      for i in range(1, 7)))


def pose(model, energy, z, *, ruffle=0):
    serials = " ".join(f"{i} {i}" for i in range(1, 7))
    return (f"MODEL {model}\nREMARK VINA RESULT: {energy} 0 0\n"
            f"REMARK SMILES c1ccccc1\nREMARK SMILES IDX {serials}\nROOT\n"
            + ring_atoms(CARBONS, "BEN", z, ruffle=ruffle) + "ENDROOT\nTORSDOF 0\nENDMDL\n")


def main():
    with TemporaryDirectory(prefix="aromatic-pool-check-") as work:
        directory = Path(work)
        pool, receptor, ccd = [directory / name for name in ("pool.pdbqt", "receptor.pdbqt", "BEN.cif")]
        ccd.write_text(component())
        receptor.write_text(ring_atoms(PHENYL_RING, "PHE", 0))
        before = pose(1, "-8.300", 8)
        edge = pose(2, "-7.800", 4)
        outside = pose(3, "-7.799", 4)
        pool.write_text(before + edge + outside)
        clean = score_aromatic_pool(pool, receptor, ccd)
        assert clean["aromatic_geometry"]["eligible_models"] == [1, 2], clean
        assert clean["orders"]["native"] == [1, 2, 3] and clean["orders"]["aromatic"] == [2, 1, 3]
        assert clean["poses"][1]["aromatic_ring_pairs"][0]["kind"] == "parallel"

        for displacement, reason in ((.85, "Nonplanar aromatic ring"), (3, "Implausible aromatic bond")):
            # An out-of-plane displacement can preserve plausible individual
            # bond lengths, so test planarity and malformed edges separately.
            pool.write_text(before + pose(2, "-7.800", 4, ruffle=displacement) + outside)
            invalid = score_aromatic_pool(pool, receptor, ccd)
            assert invalid["orders"]["aromatic"] == invalid["orders"]["native"]
            assert reason in invalid["abstentions"]["aromatic"], invalid["abstentions"]
            assert invalid["poses"][1]["aromatic_contacts"] is None

        receptor.write_text(ring_atoms(PHENYL_RING, "PHE", 0, ruffle=3))
        valid, omitted = receptor_ring_planes(receptor)
        assert not valid and "Implausible aromatic bond" in omitted[0]["reason"]
    print("Inclusive saved-score threshold, cyclic bond length and ring planarity passed.")


if __name__ == "__main__":
    main()
