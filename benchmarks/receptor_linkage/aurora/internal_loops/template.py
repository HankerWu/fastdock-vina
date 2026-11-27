"""Close an observed Aurora A backbone against fixed flanks without a ligand."""

import argparse
import json
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENMM_CPU_THREADS", "1")
import check_predictions as checks
import numpy as np
import openmm as mm
from build import atom_key, lines
from donor_check import alignment, atoms, record
from openmm import app, unit
from scipy.optimize import least_squares, minimize
from scipy.spatial.transform import Rotation
from scipy.stats import qmc

SEGMENTS = ((280, 292), (302, 306))
LOOP = checks.NEW
SUPPORT = {"1OL5", "1OL7"}
BACKBONE = {"N", "CA", "C", "O"}


def prepare_donor(work, name):
    model = atoms(record(work, name)[0])
    complete = atoms(record(work, "1OL5")[0])
    groups = {}
    for i in sorted(LOOP | {v for span in SEGMENTS for v in span}):
        src = {
            n: p
            for (j, r, n), p in model.items()
            if i == j and n not in {"P", "O1P", "O2P", "O3P"}
        }
        ref = {
            n: p
            for (j, r, n), p in complete.items()
            if i == j and n not in {"P", "O1P", "O2P", "O3P"}
        }
        if src.keys() != ref.keys():
            # Missing side chains all have a measured C-beta. Transfer the complete distal branch.
            anchors = ["N", "CA", "CB"]
            move = alignment(
                np.array([ref[n] for n in anchors]), np.array([src[n] for n in anchors])
            )
            src.update(
                {n: move(p) for n, p in ref.items() if n not in BACKBONE | {"CB"}}
            )
        groups[i] = src
    return groups


def close(left, right, groups, host, branch_search):
    anchors_left = ("CA", "C", "O")
    anchors_right = ("N", "CA", "C")
    move = alignment(
        np.array([groups[left][n] for n in anchors_left]),
        np.array([host[left, n] for n in anchors_left]),
    )
    labels = [(i, n) for i in range(left, right + 1) for n in groups[i]]
    x0 = move([groups[i][n] for i, n in labels])
    iids = np.array([i for i, n in labels])
    moving = np.array([left < i < right for i, n in labels])
    end = np.array([labels.index((right, n)) for n in anchors_right])
    target = np.array([host[right, n] for n in anchors_right])
    axes = []
    for i in range(left + 1, right + 1):
        for a, b in (("N", "CA"), ("CA", "C")):
            if (a == "N" and i == 282) or (a == "CA" and i == right):
                continue
            mask = np.array(
                [
                    j > i
                    or (j == i and (n not in {"N", "CA"} if a == "N" else n == "O"))
                    for j, n in labels
                ]
            )
            axes.append((labels.index((i, a)), labels.index((i, b)), mask))
    # Include nearby host and all nonadjacent intrasegment atoms in the soft clash objective.
    env = [
        (i, p)
        for (i, n), p in host.items()
        if i not in range(left + 1, right)
        and np.min(np.linalg.norm(x0 - p, axis=1)) < 14
    ]
    envi, envxyz = zip(*env)
    allxyz = np.concatenate([x0, np.array(envxyz)])
    h, k = np.where(
        (np.abs(iids[:, None] - np.array(envi)[None, :]) > 1) & moving[:, None]
    )
    a, b = np.where(
        np.triu(
            (np.abs(iids[:, None] - iids[None, :]) > 1)
            & moving[:, None]
            & moving[None, :],
            k=1,
        )
    )
    u = np.concatenate((h, a))
    v = np.concatenate((len(labels) + k, b))

    def coordinates(theta):
        x = x0.copy()
        for angle, (a, b, mask) in zip(theta, axes, strict=True):
            d = x[b] - x[a]
            x[mask] = (
                Rotation.from_rotvec(angle * d / np.linalg.norm(d)).apply(
                    x[mask] - x[a]
                )
                + x[a]
            )
        return x

    def objective(theta):
        x = coordinates(theta)
        allxyz[: len(labels)] = x
        vec = allxyz[u] - allxyz[v]
        sq = np.einsum("ij,ij->i", vec, vec)
        active = sq < 2.9**2
        d = np.sqrt(np.maximum(sq[active], 1e-16))
        overlap = d - 2.9
        energy = 0.5 * (
            400 * np.square(x[end] - target).sum()
            + 0.25 * theta @ theta
            + 4 * overlap @ overlap
        )
        grad = np.zeros_like(allxyz)
        grad[end] += 400 * (x[end] - target)
        f = (4 * overlap / d)[:, None] * vec[active]
        np.add.at(grad, u[active], f)
        np.add.at(grad, v[active], -f)
        slope = []
        for angle, (a, b, mask) in zip(theta, axes, strict=True):
            axis = x[b] - x[a]
            axis /= np.linalg.norm(axis)
            slope.append(
                np.sum(grad[: len(x)][mask] * np.cross(axis, x[mask] - x[a]))
                + 0.25 * angle
            )
        return energy, np.asarray(slope)

    starts = [np.zeros(len(axes))]
    candidates = []
    branches = []
    if branch_search and right - left == 4:
        starts.extend(
            np.pi
            * (2 * qmc.Sobol(len(axes), scramble=True, seed=280306).random_base2(4) - 1)
        )
    for index, initial in enumerate(starts):
        start = initial
        if len(starts) > 1:

            def residual(theta):
                return (coordinates(theta)[end] - target).ravel()

            def jac(theta):
                pts = coordinates(theta)
                derivatives = []
                for a, b, mask in axes:
                    unit_axis = pts[b] - pts[a]
                    unit_axis /= np.linalg.norm(unit_axis)
                    derivatives.append(
                        (
                            np.cross(unit_axis, pts[end] - pts[a]) * mask[end, None]
                        ).ravel()
                    )
                return np.array(derivatives).T

            root = least_squares(
                residual, start, jac=jac, max_nfev=150, ftol=1e-9, gtol=1e-8, xtol=1e-9
            )
            start = (root.x + np.pi) % (2 * np.pi) - np.pi
        candidate = minimize(
            objective,
            start,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": 1800, "ftol": 1e-12},
        )
        candidates.append(candidate)
        branches.append(
            {
                "index": index,
                "start_rad": initial.tolist(),
                "score": float(candidate.fun),
                "final_endpoint_rmsd_A": float(
                    np.linalg.norm(coordinates(candidate.x)[end] - target) / np.sqrt(3)
                ),
            }
        )
    opt = min(candidates, key=lambda opt: opt.fun)
    objective(opt.x)  # places the optimized coordinates in allxyz
    closure = float(np.sqrt(np.mean(np.sum((allxyz[end] - target) ** 2, axis=1))))
    return {
        l: p.copy()
        for l, p, m in zip(labels, allxyz[: len(labels)], moving, strict=True)
        if m
    }, {
        "branches": branches,
        "iterations": int(opt.nit),
        "optimizer_status": opt.message,
        "endpoint_rmsd_A": closure,
        "rotation_radians": opt.x.tolist(),
    }


def build(work, donor, out, branch_search, trans=False):
    start = time.perf_counter()
    scaffold = work / "prepared/2X6E-protein.pdb"
    raw = checks.atoms(scaffold)
    host = {(i, n): p for (i, r, n), p in raw.items()}
    ref = prepare_donor(work, donor)
    additions, log = {}, {}
    for lo, hi in SEGMENTS:
        positions, report = close(lo, hi, ref, host, branch_search)
        additions.update(positions)
        host.update(positions)
        log[f"{lo}-{hi}"] = report
    # Take heavy atom identity and PDB formatting from the complete observed template.
    # Substitute THR for TPO after deleting phosphate. Preserve all original host coordinates.
    keep = {}
    for l in lines(scaffold):
        keep.setdefault(int(l[22:26]), []).append(l)
    for l in record(work, "1OL5")[0]:
        if l.startswith(("ATOM  ", "HETATM")) and l[21] == "A" and l[16] in " A":
            i, n = int(l[22:26]), l[12:16].strip()
            if (
                (i, n) in additions
                and l[17:20] not in {"HOH"}
                and (l.startswith("ATOM") or l[17:20] == "TPO")
            ):
                x, y, z = additions.pop((i, n))
                name = "THR" if l[17:20] == "TPO" else l[17:20]
                keep.setdefault(i, []).append(
                    f"ATOM  {l[6:17]}{name}{l[20:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}"
                )
    assert not additions
    pdb = out / "closed.pdb"
    rows = [l for i in sorted(keep) for l in keep[i]]
    pdb.write_text(
        "\n".join(l[:6] + f"{j:5d}" + l[11:] for j, l in enumerate(rows, 1)) + "\nEND\n"
    )
    report = {
        "kinematics": log,
        "geometry_before_relaxation": checks.geometry(pdb, raw),
    }
    (out / "closure.json").write_text(json.dumps(report, indent=2) + "\n")
    initial = app.PDBFile(str(pdb))
    cpu = mm.Platform.getPlatformByName("CPU")
    ff = app.ForceField("amber14/protein.ff14SB.xml")
    random.seed(280306)
    mod = app.Modeller(initial.topology, initial.positions)
    mod.addHydrogens(ff, pH=7.0, platform=cpu)
    system = ff.createSystem(
        mod.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometer,
        constraints=None,
        removeCMMotion=False,
    )
    seq = list(mod.topology.atoms())
    for a in seq:
        if int(a.residue.id) not in LOOP and a.element != app.element.hydrogen:
            system.setParticleMass(a.index, 0)
    if trans:
        # Same stereochemistry and peptide priors as the registered Cartesian controls.
        ids = {(int(a.residue.id), a.name): a.index for a in seq}
        x = np.array(mod.positions.value_in_unit(unit.nanometer))
        v = (
            "(x2-x1)*((y3-y1)*(z4-z1)-(z3-z1)*(y4-y1))"
            "-(y2-y1)*((x3-x1)*(z4-z1)-(z3-z1)*(x4-x1))"
            "+(z2-z1)*((x3-x1)*(y4-y1)-(y3-y1)*(x4-x1))"
        )
        h = mm.CustomCompoundBondForce(4, "0.5*Khand*(V/v0-1)^2;V=" + v)
        h.addGlobalParameter("Khand", 4000.0)
        h.addPerBondParameter("v0")
        for residue in mod.topology.residues():
            i = int(residue.id)
            if i not in LOOP or residue.name == "GLY":
                continue
            quartets = [("CA", "N", "C", "CB")]
            if residue.name == "THR":
                quartets.append(("CB", "CA", "OG1", "CG2"))
            for names in quartets:
                samples = []
                for j in mod.topology.residues():
                    n = int(j.id)
                    if (
                        n not in LOOP
                        and j.name == residue.name
                        and all((n, a) in ids for a in names)
                    ):
                        pos = x[[ids[n, a] for a in names]]
                        samples.append(np.linalg.det(pos[1:] - pos[0]))
                h.addBond([ids[i, a] for a in names], [np.median(samples)])
        system.addForce(h)
        omega = mm.CustomTorsionForce("Kt*(1+cos(theta))")
        omega.addGlobalParameter("Kt", 500.0)
        for i in checks.LINKS:
            omega.addTorsion(
                ids[i, "CA"], ids[i, "C"], ids[i + 1, "N"], ids[i + 1, "CA"]
            )
        system.addForce(omega)
    itg = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, itg, cpu)
    ctx.setPositions(mod.positions)
    for stage in ("restrained", "released") if trans else ("direct",):
        if stage == "released":
            ctx.setParameter("Kt", 0)
            ctx.setParameter("Khand", 0)
        mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 3000)
        target = out / stage if trans else out
        target.mkdir(exist_ok=True, parents=True)
        assess(
            work,
            donor,
            target,
            ctx,
            system,
            seq,
            raw,
            mod,
            dict(report, stage=stage),
            start,
        )


def assess(work, donor, out, ctx, system, seq, raw, mod, report, start):
    state = ctx.getState(getPositions=True, getEnergy=True, getForces=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    movable = [
        a.index
        for a in seq
        if system.getParticleMass(a.index).value_in_unit(unit.dalton) > 0
    ]
    force = state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    report["mobile_component_force_rms_kJ_mol_nm"] = float(
        np.sqrt(np.square(force[movable]).mean())
    )
    err = []
    for f in system.getForces():
        if not isinstance(f, mm.HarmonicBondForce):
            continue
        for k in range(f.getNumBonds()):
            i, j, r, _ = f.getBondParameters(k)
            if any(seq[z].element == app.element.hydrogen for z in (i, j)) or all(
                int(seq[z].residue.id) not in LOOP for z in (i, j)
            ):
                continue
            err.append(
                (
                    10
                    * abs(
                        np.linalg.norm(pos[i] - pos[j])
                        - r.value_in_unit(unit.nanometer)
                    ),
                    str(seq[i]),
                    str(seq[j]),
                )
            )
    report["worst_heavy_bond_error_A"] = max(err)
    output = app.Modeller(mod.topology, state.getPositions())
    output.delete(
        a for a in output.topology.atoms() if a.element == app.element.hydrogen
    )
    result = out / "repaired.pdb"
    with result.open("w") as f:
        app.PDBFile.writeFile(output.topology, output.positions, f, keepIds=True)
    g = checks.geometry(result, raw)
    report["after"] = g
    now = checks.atoms(result)
    displacement = max(np.linalg.norm(now[k] - p) for k, p in raw.items())
    omega = max(180 - abs(v) for v in g["peptide_omega_deg"].values())
    valid = (
        displacement < 0.002
        and len(g["wrong_centres_in_chain"]) == 0
        and g["centres_checked_new"] == 14
        and max(err)[0] <= 0.15
        and min(g["peptide_CN_range_A"]) >= 1.20
        and max(g["peptide_CN_range_A"]) <= 1.47
        and omega <= 35
        and not g["heavy_pairs_lt_2A_graph_separation_gt_3_involving_new_residues"]
    )
    report.update(
        valid_geometry=bool(valid),
        max_host_displacement_A=float(displacement),
        max_trans_error_deg=omega,
    )
    published = out / "receptor.pdbqt"
    published.unlink(missing_ok=True)
    if valid:
        typed = out / "typed.pdbqt"
        with (out / "typing.log").open("w") as f:
            typing = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "meeko.cli.mk_prepare_receptor",
                    "--read_pdb",
                    str(result),
                    "-p",
                    str(typed),
                ],
                check=False,
                stdout=f,
                stderr=subprocess.STDOUT,
            )
        report["typing_exit_code"] = typing.returncode
        if typing.returncode == 0:
            qt = lines(work / "prepared/2X6E.pdbqt")
            extra = [l for l in lines(typed) if atom_key(l)[1] in LOOP]
            assert {
                (i, n)
                for _, i, n in (
                    atom_key(l) for l in extra if l.split()[-1] not in ("H", "HD")
                )
            } == {(i, n) for i, r, n in now if i in LOOP}
            highest = max(int(l[6:11]) for l in qt)
            qt.extend(
                l[:6] + f"{j:5d}" + l[11:] for j, l in enumerate(extra, highest + 1)
            )
            (out / "receptor.pdbqt").write_text("\n".join(qt) + "\n")
    if report["stage"] == "released":
        canonical = out.parents[1] / "receptor.pdbqt"
        if published.exists():
            shutil.copyfile(published, canonical)
        else:
            canonical.unlink(missing_ok=True)
    report["wall_s"] = time.perf_counter() - start
    (out / "template.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        donor,
        report["valid_geometry"],
        g["peptide_CN_range_A"],
        max(err),
        omega,
        flush=True,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", type=Path)
    ap.add_argument("donor", choices=sorted(SUPPORT))
    ap.add_argument(
        "--initial",
        action="store_true",
        help="Reproduce the initial single-start attempt",
    )
    ap.add_argument(
        "--trans-model",
        action="store_true",
        help="Evaluate restrained geometry and subsequent release",
    )
    args = ap.parse_args()
    dest = args.work.resolve() / f"template{args.donor}"
    if args.initial:
        dest = dest / "initial"
    if args.trans_model:
        dest = dest / "trans_model"
    dest.mkdir(exist_ok=True, parents=True)
    build(args.work.resolve(), args.donor, dest, not args.initial, args.trans_model)
