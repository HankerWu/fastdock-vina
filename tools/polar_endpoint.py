"""Frozen-coordinate Coulomb and OBC-II polar transfer using native PDBQT charges."""

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import unit
from scipy.spatial import cKDTree, distance
from select_zinc_pose import read_models

# OpenMM OBC2 standard mbondi2 radii (Å) and HCT scale parameters.
PARAMETERS = {
    "H": (1.2, 0.85),
    "C": (1.7, 0.72),
    "N": (1.55, 0.79),
    "O": (1.5, 0.85),
    "F": (1.5, 0.88),
    "P": (1.85, 0.86),
    "S": (1.8, 0.96),
    "Cl": (1.7, 0.8),
}
ELEMENT = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H", "HS": "H"}
COULOMB = 138.935456 * 10 / 4.184  # kcal Å / mol e²


@dataclass
class Atoms:
    records: tuple
    xyz: np.ndarray
    charges: np.ndarray
    radii: np.ndarray
    scales: np.ndarray
    omitted_zinc: int


def atoms(lines, *, receptor=False):
    records, xyz, charges, elements = [], [], [], []
    omitted = 0
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial, atom_type = int(line[6:11]), line[77:].strip()
        if receptor and atom_type == "Zn":
            omitted += 1
            continue
        element = ELEMENT.get(atom_type, atom_type)
        if element not in PARAMETERS:
            raise ValueError(f"Unsupported PDBQT atom type {atom_type}")
        records.append((serial, atom_type))
        xyz.append([float(line[i : i + 8]) for i in (30, 38, 46)])
        charges.append(float(line[70:76]))
        elements.append(element)
    if not xyz or len({i for i, _ in records}) != len(xyz):
        raise ValueError("Empty or duplicate-serial PDBQT atom set")
    xyz, charges = np.asarray(xyz), np.asarray(charges)
    if not (np.isfinite(xyz).all() and np.isfinite(charges).all()):
        raise ValueError("Nonfinite coordinates or charges")
    elements = np.asarray(elements)
    radii, scales = np.transpose([PARAMETERS[e] for e in elements])
    hs, heavy = np.flatnonzero(elements == "H"), np.flatnonzero(elements != "H")
    if len(hs):
        d, nearest = cKDTree(xyz[heavy]).query(xyz[hs])
        if np.any(d > 1.55):
            raise ValueError("Cannot assign a covalent parent to a PDBQT hydrogen")
        radii[hs[elements[heavy[nearest]] == "N"]] = 1.3
    return Atoms(tuple(records), xyz, charges, radii, scales, omitted)


class PolarContext:
    def __init__(self, molecules, platform="CPU", *, solute_dielectric=4):
        system = mm.System()
        force = mm.GBSAOBCForce()
        force.setSoluteDielectric(solute_dielectric)
        force.setSolventDielectric(78.5)
        force.setSurfaceAreaEnergy(0)
        force.setNonbondedMethod(mm.GBSAOBCForce.NoCutoff)
        for molecule in molecules:
            for q, r, s in zip(
                molecule.charges, molecule.radii, molecule.scales, strict=True
            ):
                system.addParticle(1)
                force.addParticle(float(q), float(r / 10), float(s))
        system.addForce(force)
        self.context = mm.Context(
            system,
            mm.VerletIntegrator(0.001),
            mm.Platform.getPlatformByName(platform),
            {"Threads": "1"} if platform == "CPU" else {},
        )

    def energy(self, xyz):
        self.context.setPositions(np.asarray(xyz) / 10)
        return (
            self.context.getState(getEnergy=True)
            .getPotentialEnergy()
            .value_in_unit(unit.kilocalories_per_mole)
        )


class PolarEndpoint:
    """One fixed receptor and one prepared ligand with invariant atom identity."""

    def __init__(self, receptor, ligand, platform="CPU"):
        self.receptor, self.ligand = receptor, ligand
        apo = PolarContext([receptor], platform)
        self.apo_energy = apo.energy(receptor.xyz)
        self.bound = PolarContext([receptor, ligand], platform)
        self.free = PolarContext([ligand], platform)

    def evaluate(self, pose):
        if (
            pose.records != self.ligand.records
            or not np.array_equal(pose.charges, self.ligand.charges)
            or not np.array_equal(pose.radii, self.ligand.radii)
        ):
            raise ValueError(
                "PDBQT pose changes prepared atom identities, charges or radii"
            )
        r = self.receptor
        d = distance.cdist(r.xyz, pose.xyz)
        if not np.all(d > 0):
            raise ValueError("Coincident receptor and ligand charges")
        vacuum = COULOMB * np.sum((r.charges[:, None] * pose.charges) / d)
        complex_g = self.bound.energy(np.concatenate((r.xyz, pose.xyz)))
        ligand_g = self.free.energy(pose.xyz)
        delta_g = complex_g - self.apo_energy - ligand_g
        return {
            "coulomb": float(vacuum) / 4,
            "complex_gb": complex_g,
            "apo_gb": self.apo_energy,
            "ligand_gb": ligand_g,
            "delta_gb": delta_g,
            "polar_transfer": float(vacuum) / 4 + delta_g,
        }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("receptor", type=Path)
    p.add_argument("ligand", type=Path)
    p.add_argument("poses", type=Path)
    p.add_argument("--platform", choices=("CPU", "Reference"), default="CPU")
    a = p.parse_args()
    receptor = atoms(a.receptor.read_text().splitlines(), receptor=True)
    calc = PolarEndpoint(receptor, atoms(a.ligand.read_text().splitlines()), a.platform)
    writer = None
    for model, lines in read_models(a.poses.read_text().splitlines()).items():
        row = {
            "model": model,
            "omitted_receptor_zinc": receptor.omitted_zinc,
            **calc.evaluate(atoms(lines)),
        }
        if writer is None:
            writer = csv.DictWriter(sys.stdout, fieldnames=row, lineterminator="\n")
            writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
