"""Labelled diagnostic of the gradient lookahead on two earlier development failures.

Uses SciPy's trilinear interpolation and an SVD of the full Cartesian rigid-body
Jacobian, independently of the prediction kernel. No evaluation labels influence
the scan; the two deliberately oracle-selected poses are marked in the output.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
SCAN = HERE.parents[1] / "rigid_translation_scan"
sys.path.insert(0, str(SCAN))
from run import paired, parse_pose, rms, save
from scan import read_ligand, receptor_maps, translate_search, write_start

sys.path.insert(0, str(SCAN.parents[1] / "tools"))


def atom_field(field, xyz, origin, spacing, gradient=False):
    """Take slopes as face-value differences; trilinear fields are affine along each axis."""
    idx = ((xyz - origin) / spacing).T

    def sample(at):
        return map_coordinates(
            field, at, order=1, prefilter=False, mode="constant", cval=np.inf
        )

    energy = sample(idx)
    if not gradient:
        return energy
    slopes = []
    for axis, size in enumerate(field.shape):
        lo = idx.copy()
        lo[axis] = np.minimum(np.floor(idx[axis]), size - 2)
        hi = lo.copy()
        hi[axis] += 1
        slopes.append((sample(hi) - sample(lo)) / spacing)
    return energy, np.stack(slopes, axis=-1)


def jacobian(offsets):
    """Map (translation vector, rotation vector) to lab-frame atom displacements."""
    axes = np.eye(3)
    turn = np.cross(axes[:, None, :], offsets[None, :, :]).transpose(1, 2, 0)
    return np.concatenate((np.broadcast_to(axes, turn.shape), turn), axis=2).reshape(
        -1, 6
    )


def tangent_step(offsets, gradients, length):
    jac = jacobian(offsets)
    flat = gradients.reshape(len(gradients), -1)
    step = -flat @ np.linalg.pinv(jac).T
    scale = np.linalg.norm(step @ jac.T, axis=1) / np.sqrt(len(offsets))
    return np.divide(
        length * step,
        scale[:, None],
        out=np.zeros_like(step),
        where=scale[:, None] > 1e-12,
    )


def lookahead(fields, origin, spacing, offsets, types, positions):
    centroid = offsets.mean(axis=0)
    local = offsets - centroid
    gradients = []
    initial = np.zeros(len(positions))
    for off, typ in zip(offsets, types):
        energy, grad = atom_field(fields[typ], positions + off, origin, spacing, True)
        initial += energy
        gradients.append(grad)
    step = tangent_step(local, np.stack(gradients, axis=1), spacing)
    turn = Rotation.from_rotvec(step[:, 3:]).as_matrix()
    moved = (positions + centroid + step[:, :3])[:, None, :] + np.einsum(
        "bki,ai->bak", turn, local
    )
    trial = np.zeros(len(positions))
    for j, typ in enumerate(types):
        trial += atom_field(fields[typ], moved[:, j], origin, spacing)
    accepted = trial < initial
    step[~accepted] = 0
    return initial, np.where(accepted, trial, initial), step


def transform(rotated, positions, steps, heavy):
    centroid = rotated[heavy].mean(0)
    rotation = Rotation.from_rotvec(steps[:, 3:]).as_matrix()
    return (positions + centroid + steps[:, :3])[:, None, :] + np.einsum(
        "bki,ai->bak", rotation, rotated - centroid
    )


def rmsds(xyz, refs):
    return np.sqrt(
        np.minimum.reduce(
            [
                np.square(xyz - reference).sum(axis=-1).mean(axis=-1)
                for reference in refs
            ]
        )
    )


def run_case(args, meta, query, rec, seed, ligand, fields, origin, spacing):
    center = np.array(meta[rec]["center"])
    rigid, shift, *_ = paired.alignment(meta[query], meta[rec])
    refs = np.asarray(meta[query]["mapped"]) @ rigid + shift
    indices = {serial: i for i, (_, serial, _) in enumerate(ligand[1])}
    idx = [indices[s] for s in meta[query]["serial"]]
    heavy = ligand[3]
    blocks = []
    all_xyz = []
    errors = []

    def observe(orient, rotated, positions, parent_energy):
        before, after, step = lookahead(
            fields, origin, spacing, rotated[heavy], ligand[4], positions
        )
        np.testing.assert_allclose(parent_energy, before, rtol=1e-13, atol=1e-12)
        if args.check_kernel:
            from rigid_lookahead import rigid_trial

            prod = rigid_trial(
                fields,
                origin,
                spacing,
                rotated[heavy],
                ligand[4],
                positions,
                "gradient",
                np.random.default_rng(0),
            )
            prod_step = np.column_stack((prod["translation"], prod["rotvec"]))
            for field, expected, actual in (
                ("parent", before, prod["parent_energy"]),
                ("accepted_energy", after, prod["energy"]),
                ("rigid_step", step, prod_step),
            ):
                np.testing.assert_allclose(expected, actual, rtol=2e-11, atol=2e-10)
                errors.append((field, np.max(np.abs(expected - actual))))
        base = rotated + positions[:, None, :]
        end = transform(rotated, positions, step, heavy)
        blocks.append(
            (
                np.full(len(before), orient),
                before,
                after,
                rmsds(base[:, idx], refs),
                rmsds(end[:, idx], refs),
            )
        )
        all_xyz.append((base, end))

    translate_search(fields, origin, spacing, center, ligand, seed, "lattice", observe)
    orient, before, after, d0, d1 = (np.concatenate(values) for values in zip(*blocks))
    stages = np.array([np.concatenate(values) for values in zip(*all_xyz)])
    low = {
        i: np.flatnonzero(orient == i)[np.argmin(after[orient == i])]
        for i in np.unique(orient)
    }
    released = []
    for k in sorted(low.values(), key=lambda j: (after[j], orient[j])):
        if all(
            np.square(stages[1, k, heavy] - stages[1, j, heavy]).sum(axis=1).mean() >= 1
            for j in released
        ):
            released.append(k)
        if len(released) == 9:
            break
    pre, post = d0 <= 2, d1 <= 2
    keys = [("closest_before", int(d0.argmin()))]
    if pre.any():
        old_near = int(np.where(pre, before, np.inf).argmin())
        keys.append(("lowest_before_within2", old_near))
        keys.append(("new_winner_same_orientation", int(low[orient[old_near]])))
    if post.any():
        keys.append(
            ("lowest_after_within2", int(np.where(post, after, np.inf).argmin()))
        )
    keys.append(("blind_lookahead_minimum", int(np.argmin(after))))
    out = args.work / f"{query}_{rec}_{seed}"
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    refine_ids = {(j, state) for _, j in keys for state in (0, 1)}
    refine_ids.update((j, 1) for j in released)
    for j, state in sorted(refine_ids):
        path = out / f"{j}_{state}.pdbqt"
        write_start(path, *ligand[:2], stages[state, j])
        paths[j, state] = path
    log = subprocess.check_output(
        [
            str(args.probe),
            "refine",
            str(args.dev / rec / "receptor.pdbqt"),
            str(args.dev / query / "ligand.pdbqt"),
            *map(str, center),
            *[path.name for path in paths.values()],
        ],
        text=True,
        cwd=out,
    )
    (out / "oracle.log").write_text(log)
    refined = {
        line.split()[1]: float(line.split()[3])
        for line in log.splitlines()
        if line.startswith("result ")
    }
    final_coords = {
        j: parse_pose(Path(str(paths[j, 1]) + ".refined"), meta[query]["serial"])
        for j in released
    }
    final_order = []
    for j in sorted(released, key=lambda i: refined[paths[i, 1].name]):
        if all(
            np.square(final_coords[j] - final_coords[k]).sum(axis=1).mean() >= 1
            for k in final_order
        ):
            final_order.append(j)
    final_rank = {j: k + 1 for k, j in enumerate(final_order)}
    info = dict(query=query, receiver=rec, seed=seed)
    cell = dict(
        **info,
        probes=len(before),
        initial_near=int(pre.sum()),
        retained_near=int((pre & post).sum()),
        lost_near=int((pre & ~post).sum()),
        gained_near=int((~pre & post).sum()),
        final_near=int(post.sum()),
        near_orientation_winners=int(post[list(low.values())].sum()),
        near_refine_starts=int(post[released].sum()),
        reported=len(final_rank),
        top1_rmsd=rms(final_coords[final_order[0]], refs),
        top5_rmsd=min(rms(final_coords[j], refs) for j in final_order[:5]),
        oracle_rmsd=min(rms(final_coords[j], refs) for j in final_order),
    )
    rows = []
    atoms = []
    for label, j in keys:
        scope = orient == orient[j]
        row = dict(
            **info,
            diagnostic=label,
            orientation=int(orient[j]),
            orientation_index=int(scope[:j].sum()),
            rank_before=1 + int((before < before[j]).sum()),
            rank_after=1 + int((after < after[j]).sum()),
            in_orientation_before=1 + int((before[scope] < before[j]).sum()),
            in_orientation_after=1 + int((after[scope] < after[j]).sum()),
            in_refine_starts=int(j in released),
            final_rank=final_rank.get(j, ""),
        )
        for state, name, e, d in [(0, "before", before, d0), (1, "after", after, d1)]:
            path = paths[j, state]
            xyz = parse_pose(Path(str(path) + ".refined"), meta[query]["serial"])
            row.update(
                {
                    name + "_coarse": e[j],
                    name + "_rmsd": d[j],
                    name + "_refined_energy": refined[path.name],
                    name + "_refined_rmsd": rms(xyz, refs),
                }
            )
        rows.append(row)
        if label in ("lowest_before_within2", "new_winner_same_orientation"):
            heavy_serials = [ligand[1][a][1] for a in heavy]
            for k, (typ, serial) in enumerate(zip(ligand[4], heavy_serials)):
                coords = stages[:, j, heavy[k], :]
                prior, trial = atom_field(fields[typ], coords, origin, spacing)
                atoms.append(
                    dict(
                        **info,
                        diagnostic=label,
                        serial=serial,
                        xs=typ,
                        prior=prior,
                        trial=trial,
                        change=trial - prior,
                    )
                )
    maxima = [
        dict(
            **info,
            observable=field,
            max_abs_difference=max(x for name, x in errors if name == field),
        )
        for field in sorted({name for name, x in errors})
    ]
    winners = [
        dict(
            **info,
            slot=k + 1,
            orientation=int(orient[j]),
            coarse=after[j],
            rmsd=d1[j],
            final_rank=final_rank.get(j, ""),
            final_rmsd=rms(final_coords[j], refs),
        )
        for k, j in enumerate(released)
    ]
    return cell, rows, winners, atoms, maxima


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("dev", "work", "probe", "vina"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument(
        "--check-kernel",
        action="store_true",
        help="Compare the production kernel with the independent SVD reference on every parent",
    )
    args = p.parse_args()
    args.dev, args.probe, args.vina = (
        args.dev.resolve(),
        args.probe.resolve(),
        args.vina.resolve(),
    )
    args.work.mkdir(parents=True, exist_ok=True)
    meta = json.loads((args.dev / "metadata.json").read_text())
    cells, rows, winners, atom_changes, maxima = [], [], [], [], []
    for query, rec in [("5EPS", "5C7N"), ("3WKE", "4OD0")]:
        fields, origin, spacing = receptor_maps(
            args.vina,
            args.dev / rec / "receptor.pdbqt",
            meta[rec]["center"],
            args.work / "maps" / rec,
        )
        ligand = read_ligand(args.probe, args.dev / query / "ligand.pdbqt")
        for seed in (43103, 43117):
            cell, data, chosen, atoms, errors = run_case(
                args, meta, query, rec, seed, ligand, fields, origin, spacing
            )
            cells.append(cell)
            rows.extend(data)
            winners.extend(chosen)
            atom_changes.extend(atoms)
            maxima.extend(errors)
    for name, data in [
        ("coverage", cells),
        ("counterfactual", rows),
        ("blind_starts", winners),
        ("atom_changes", atom_changes),
        ("kernel_errors", maxima),
    ]:
        if data:
            save(args.work / (name + ".csv"), data)


if __name__ == "__main__":
    main()
