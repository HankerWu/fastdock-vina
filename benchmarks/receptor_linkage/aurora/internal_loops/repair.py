"""Keep the measured scaffold fixed and restore L-peptide stereochemistry to a loop proposal."""

import argparse
import json
import math
import os
import random
import time
from pathlib import Path

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
import numpy as np
import openmm as mm
from openmm import app, unit

LOOP = set(range(281, 292)) | {303, 304, 305}
BACKBONE = ("CA", "N", "C", "CB")
THREONINE = ("CB", "CA", "OG1", "CG2")


def signed_volume(xyz):
    return float(np.linalg.det(xyz[1:] - xyz[0]))


def model(work, seed, rebuild_h=False):
    start = time.perf_counter()
    out = work / f"loop{seed}"
    pdb = app.PDBFile(str(out / "completed.pdb"))
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    cpu = mm.Platform.getPlatformByName("CPU")
    random.seed(seed)
    m = app.Modeller(pdb.topology, pdb.positions)
    m.addHydrogens(ff, pH=7, platform=cpu)
    system = ff.createSystem(
        m.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometer,
        constraints=None,
        removeCMMotion=False,
    )
    atoms = list(m.topology.atoms())
    atom_id = {(int(a.residue.id), a.name): a.index for a in atoms}
    xyz = np.asarray(m.positions.value_in_unit(unit.nanometer))
    for a in atoms:
        if int(a.residue.id) not in LOOP and a.element != app.element.hydrogen:
            system.setParticleMass(a.index, 0)
    # V is the signed determinant in nm^3; V/v0 is dimensionless.
    v = (
        "(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))"
        "-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))"
        "+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))"
    )
    hand = mm.CustomCompoundBondForce(4, "0.5*Khand*(V/v0-1)^2;V=" + v)
    hand.addGlobalParameter("Khand", 4000.0)
    hand.addPerBondParameter("v0")
    n = 0
    centers = []
    for residue in m.topology.residues():
        i = int(residue.id)
        if i not in LOOP or residue.name == "GLY":
            continue
        for quartet in [BACKBONE, THREONINE] if residue.name == "THR" else [BACKBONE]:
            templates = []
            for other in m.topology.residues():
                j = int(other.id)
                if (
                    j not in LOOP
                    and other.name == residue.name
                    and all((j, a) in atom_id for a in quartet)
                ):
                    templates.append(
                        signed_volume(xyz[[atom_id[j, a] for a in quartet]])
                    )
            target = float(np.median(templates))
            if abs(target) < 0.001 or not templates:
                raise ValueError(("no stereochemical reference", i, quartet))
            indices = [atom_id[i, a] for a in quartet]
            hand.addBond(indices, [target])
            centers.append((i, quartet[0], indices, target))
            n += 1
    system.addForce(hand)
    # Use a trans initial peptide model, including the one X–Pro junction.
    omega = mm.CustomTorsionForce("Komega*(1+cos(theta))")
    omega.addGlobalParameter("Komega", 500.0)
    peptides = []
    for i in sorted(LOOP | {280, 302}):
        ii = [
            atom_id[i, "CA"],
            atom_id[i, "C"],
            atom_id[i + 1, "N"],
            atom_id[i + 1, "CA"],
        ]
        omega.addTorsion(*ii)
        peptides.append(ii)
    system.addForce(omega)
    itg = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, itg, cpu)
    ctx.setPositions(m.positions)
    stages = []
    ch_moves = []
    neighbors = {a.index: [] for a in atoms}
    for a, b in m.topology.bonds():
        neighbors[a.index].append(b.index)
        neighbors[b.index].append(a.index)

    for kh, ko, budget in [(4000, 500, 750), (4000, 500, 1500), (0, 0, 1200)]:
        if kh == 0 and rebuild_h:
            # A heavy signed volume can be correct with the old H on the wrong face.
            coords = np.array(
                ctx.getState(getPositions=True)
                .getPositions(asNumpy=True)
                .value_in_unit(unit.nanometer)
            )
            for i, name, indices, target in centers:
                central, *heavy = indices
                (hydrogen,) = [
                    j
                    for j in neighbors[central]
                    if atoms[j].element == app.element.hydrogen
                ]
                axis = coords[heavy] - coords[central]
                direction = -np.sum(
                    axis / np.linalg.norm(axis, axis=1)[:, None], axis=0
                )
                old = coords[hydrogen] - coords[central]
                projected = direction / np.linalg.norm(direction) * np.linalg.norm(old)
                ch_moves.append(
                    {
                        "residue": i,
                        "atom": name,
                        "displacement_A": 10 * np.linalg.norm(projected - old),
                        "cosine": float(np.dot(old, projected) / np.dot(old, old)),
                    }
                )
                coords[hydrogen] = coords[central] + projected
            ctx.setPositions(coords * unit.nanometer)
        ctx.setParameter("Khand", kh)
        ctx.setParameter("Komega", ko)
        mm.LocalEnergyMinimizer.minimize(ctx, 10.0, budget)
        state = ctx.getState(getPositions=True, getEnergy=True, getForces=True)
        coords = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        residual = state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
        mobile = [
            a.index
            for a in atoms
            if system.getParticleMass(a.index).value_in_unit(unit.dalton) > 0
        ]
        stages.append(
            {
                "hand_k": kh,
                "omega_k": ko,
                "iteration_cap": budget,
                "potential_kJ_mol": state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                ),
                "rms_per_mobile_component_kJ_mol_nm": float(
                    np.sqrt(np.mean(residual[mobile] ** 2))
                ),
                "volume_A3": {
                    f"{i}-{name}": 1000 * signed_volume(coords[idx])
                    for i, name, idx, t in centers
                },
                "min_relative_hand": min(
                    signed_volume(coords[idx]) / t for i, name, idx, t in centers
                ),
            }
        )
    result = app.Modeller(m.topology, state.getPositions())
    result.delete(
        [a for a in result.topology.atoms() if a.element == app.element.hydrogen]
    )
    prefix = "repaired" if rebuild_h else "first_release"
    with (out / (prefix + ".pdb")).open("w") as f:
        app.PDBFile.writeFile(result.topology, result.positions, f, keepIds=True)

    def torsion(idx):
        u, v, w = np.diff(coords[idx], axis=0)
        a, b = np.cross(u, v), np.cross(v, w)
        return math.degrees(
            math.atan2(np.dot(np.cross(a, b), v) / np.linalg.norm(v), np.dot(a, b))
        )

    report = {
        "seed": seed,
        "stereocenters": n,
        "hydrogen_reset": ch_moves,
        "final_omegas_degrees": list(map(torsion, peptides)),
        "stages": stages,
        "elapsed_s": time.perf_counter() - start,
    }
    (out / (prefix + ".json")).write_text(json.dumps(report, indent=2) + "\n")
    print(
        report["elapsed_s"],
        stages[-1]["min_relative_hand"],
        report["final_omegas_degrees"],
        flush=True,
    )
    return stages[-1]["min_relative_hand"] > 0


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("seed", type=int)
    p.add_argument("--reset-h", action="store_true")
    a = p.parse_args()
    model(a.work.resolve(), a.seed, a.reset_h)
