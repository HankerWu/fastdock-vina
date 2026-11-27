"""Independently check the raw Aurora A sequence and the ff14SB gap bond energy/force."""

import argparse
import hashlib
import itertools
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit


def deposition(path):
    lines = path.read_text().splitlines()
    missing = {}
    observed = {}
    seqres = []
    for line in lines:
        record = re.fullmatch(r"REMARK 465\s+([A-Z]{3}) A\s+(\d+)\s*", line)
        if record:
            name, number = record.groups()
            missing[int(number)] = name
        if line.startswith("ATOM  ") and line[21] == "A":
            assert line[26] == " "
            observed[int(line[22:26])] = line[17:20]
        if line.startswith("SEQRES") and line[11] == "A":
            seqres.extend(line[19:70].split())
    assert not (missing.keys() & observed.keys())
    all_residues = observed | missing
    sequence = [
        all_residues[i] for i in range(min(all_residues), max(all_residues) + 1)
    ]
    assert sequence == seqres
    return observed, missing, len(seqres)


def atom_type(xml, residue, atom):
    names = ("HID", "HIE", "HIP") if residue == "HIS" else (residue,)
    selected = {
        xml.find(f"./Residues/Residue[@name='{name}']/Atom[@name='{atom}']").get("type")
        for name in names
    }
    assert len(selected) == 1
    return selected.pop()


def native_term(types, xml):
    for node in xml.findall("./HarmonicBondForce/Bond"):
        if {node.get("type1"), node.get("type2")} == set(types):
            return float(node.attrib["length"]), float(node.attrib["k"])
    raise KeyError(types)


def check_term(a, b, r0, k):
    # Isolate the single term; this is never the total force on a protein atom.
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, r0, k)
    system = mm.System()
    system.addParticle(12)
    system.addParticle(14)
    system.addForce(force)
    step = mm.VerletIntegrator(0.001)
    context = mm.Context(system, step, mm.Platform.getPlatformByName("Reference"))
    context.setPositions([a, b])
    data = context.getState(energy=True, forces=True)
    energy = data.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = data.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    separation = np.linalg.norm(b - a)
    expected_force = k * (separation - r0)
    expected_energy = 0.5 * k * (separation - r0) ** 2
    assert np.allclose(
        forces[0], expected_force * (b - a) / separation, rtol=1e-11, atol=1e-9
    )
    assert np.allclose(forces[0], -forces[1]) and math.isclose(
        energy, expected_energy, rel_tol=1e-11
    )
    return {
        "distance_A": 10 * separation,
        "equilibrium_A": 10 * r0,
        "k_kJ_mol_nm2": k,
        "energy_kJ_mol": energy,
        "force_magnitude_kJ_mol_nm": float(np.linalg.norm(forces[0])),
    }


def verify(work):
    xmlpath = Path(app.__file__).parent / "data/amber14/protein.ff14SB.xml"
    xml = ET.parse(xmlpath)
    report = {
        "amber14_protein_ff14sb_sha256": hashlib.sha256(
            xmlpath.read_bytes()
        ).hexdigest(),
        "openmm_version": mm.version.version,
        "structures": {},
    }
    for pdb in ("2J4Z", "2X6E"):
        observed, missing, length = deposition(work / "raw" / f"{pdb}.pdb")
        model = app.PDBFile(str(work / "prepared" / f"{pdb}-protein.pdb"))
        residues = list(model.topology.residues())
        assert len({r.chain for r in residues}) == 1 and residues[0].chain.id == "A"
        ids = [int(r.id) for r in residues]
        assert ids == sorted(observed)
        seq_gaps = {(a, b) for a, b in itertools.pairwise(ids) if b > a + 1}
        topology_gaps = set()
        xyz = model.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        entries = []
        for a, b in model.topology.bonds():
            if a.name == "N" and b.name == "C":
                a, b = b, a
            if (a.name, b.name) != ("C", "N"):
                continue
            i, j = int(a.residue.id), int(b.residue.id)
            if j <= i + 1:
                continue
            assert (i, j) in seq_gaps
            topology_gaps.add((i, j))
            names = [missing[n] for n in range(i + 1, j)]
            types = (
                atom_type(xml, a.residue.name, "C"),
                atom_type(xml, b.residue.name, "N"),
            )
            r0, k = native_term(types, xml)
            entries.append(
                {
                    "from": i,
                    "to": j,
                    "missing": names,
                    "types": types,
                    **check_term(xyz[a.index], xyz[b.index], r0, k),
                }
            )
        assert topology_gaps == seq_gaps
        report["structures"][pdb] = {
            "SEQRES_count": length,
            "observed_residues": len(ids),
            "artificial_links": entries,
        }
    assert [len(v["artificial_links"]) for v in report["structures"].values()] == [0, 2]
    return report


def compare_exports(old, new):
    """Compare chemical graphs and atom coordinates after full preparation runs."""
    results = {}
    for pdb in ("2X6E", "2J4Z"):

        def read(root, pdb):
            base = root / "prepared"
            m = app.PDBFile(str(base / f"{pdb}-protein.pdb"))
            atoms = [
                (a.residue.chain.id, a.residue.id, a.residue.insertionCode, a.name)
                for a in m.topology.atoms()
            ]
            bonds = {
                tuple(sorted((atoms[a.index], atoms[b.index])))
                for a, b in m.topology.bonds()
            }
            meta = json.loads((base / f"{pdb}.json").read_text())
            return (
                m.topology.getNumChains(),
                atoms,
                bonds,
                m.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
                hashlib.sha256((base / f"{pdb}.pdbqt").read_bytes()).hexdigest(),
                meta,
            )

        a, b = read(old, pdb), read(new, pdb)
        removed = a[2] - b[2]
        assert (
            a[1] == b[1]
            and np.array_equal(a[3], b[3])
            and a[4] == b[4]
            and not b[2] - a[2]
        )
        expected = {("280", "292"), ("302", "306")} if pdb == "2X6E" else set()
        assert {
            tuple(sorted((left[1], right[1]), key=int)) for left, right in removed
        } == expected
        assert {
            (g["after"], g["before"]) for g in b[5]["unresolved_internal_segments"]
        } == expected
        results[pdb] = {
            "segments_before_after": [a[0], b[0]],
            "identical_heavy_coordinates": len(a[1]),
            "removed_edges": sorted(removed),
            "rigid_pdbqt_sha256_both": a[4],
            "deposited_gaps": b[5]["reported_gaps"],
        }
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--exports", nargs=2, type=Path, metavar=("PARENT", "FIXED"))
    args = parser.parse_args()
    result = compare_exports(*args.exports) if args.exports else verify(args.work)
    print(json.dumps(result, indent=2))
