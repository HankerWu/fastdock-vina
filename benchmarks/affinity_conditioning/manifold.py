"""Project archived ligand-6 poses onto ligand 5's actual PDBQT rigid/torsional manifold."""
import argparse
import csv
import itertools
import json
import time
from pathlib import Path

import numpy as np

from verify import BENCH, REC, Fields, Ligand, active_pools, members, modes


def tree(pdbqt, ligand):
    """Read direct rigid groups and the complete downstream set of each rotatable bond."""
    graph_index = {serial: i for i, serial in enumerate(ligand.order)}
    parents, axes, groups = [-1], [None], [[]]
    stack = [0]
    for line in pdbqt.splitlines():
        tokens = line.split()
        if not tokens:
            continue
        if tokens[0] == "BRANCH":
            axes.append(tuple(graph_index[int(s)] for s in tokens[1:3]))
            parents.append(stack[-1])
            stack.append(len(groups))
            groups.append([])
        elif tokens[0] == "ENDBRANCH":
            assert axes[stack.pop()] == tuple(graph_index[int(s)] for s in tokens[1:3])
        elif tokens[0] == "ATOM":
            groups[stack[-1]].append(graph_index[int(line[6:11])])
    assert stack == [0] and len(axes) == ligand.tors + 1
    assert sorted(itertools.chain.from_iterable(groups)) == list(range(len(ligand.order)))
    descendants = [list(g) for g in groups]
    local = [list(g) for g in groups]
    for i in range(len(axes)-1, 0, -1):
        a, b = axes[i]
        assert b in groups[i] and a in groups[parents[i]]
        descendants[parents[i]].extend(descendants[i])
        local[parents[i]].append(b)
        local[i].append(a)
    return axes, groups, local, descendants


def align(coords, target, indices):
    """Independent proper least-squares rigid motion for each batch entry."""
    source, dest = coords[:, indices], target[:, indices]
    source_mean, dest_mean = source.mean(axis=1, keepdims=True), dest.mean(axis=1, keepdims=True)
    covariance = (source-source_mean).transpose(0, 2, 1) @ (dest-dest_mean)
    u, _, vt = np.linalg.svd(covariance)
    u[:, :, -1] *= np.linalg.det(u @ vt)[:, None]
    return (coords-source_mean) @ (u @ vt) + dest_mean


def rotate_branch(coords, target, axis, fit, move):
    """Rotate a complete subtree, fitting selected points about its existing physical axis."""
    a, b = axis
    origin = coords[:, a:a+1]
    unit = coords[:, b:b+1] - origin
    unit /= np.linalg.norm(unit, axis=2, keepdims=True)
    def project(points):
        rel = points - origin
        axial = (rel*unit).sum(axis=2, keepdims=True)*unit
        return axial, rel-axial
    _, radial = project(coords[:, fit])
    dest = target[:, fit] - origin
    c = (radial * dest).sum(axis=(1, 2))
    s = (np.cross(unit, radial) * dest).sum(axis=(1, 2))
    radius = np.hypot(c, s)
    # Degenerate axial-only groups give no rotation preference.
    cosine = np.divide(c, radius, out=np.ones_like(c), where=radius>1e-13)[:, None, None]
    sine = np.divide(s, radius, out=np.zeros_like(s), where=radius>1e-13)[:, None, None]
    axial, radial = project(coords[:, move])
    coords[:, move] = origin + axial + cosine*radial + sine*np.cross(unit, radial)


def solve(reference, targets, architecture, start):
    axes, _, local, descendants = architecture
    all_atoms = np.arange(len(reference))
    coords = np.broadcast_to(reference, targets.shape).copy()
    if start == "tree":
        coords = align(coords, targets, local[0])
        for i in range(1, len(axes)):
            rotate_branch(coords, targets, axes[i], local[i], descendants[i])
    coords = align(coords, targets, all_atoms)
    def mse(x, y):
        return np.square(x-y).sum(axis=2).mean(axis=1)
    error = mse(coords, targets)
    initial = np.sqrt(error)
    active = np.arange(len(coords))
    previous_small = np.zeros(len(coords), bool)
    stopped = np.full(len(coords), 100)
    converged = np.zeros(len(coords), bool)
    for epoch in range(1, 101):
        if not len(active):
            break
        x = coords[active]
        y = targets[active]
        for i in range(1, len(axes)):
            rotate_branch(x, y, axes[i], descendants[i], descendants[i])
        x = align(x, y, all_atoms)
        curr = mse(x, y)
        gain = error[active] - curr
        if np.min(gain) < -1e-10:
            raise ValueError("Analytic rotation/rigid fit increased geometric objective")
        use = gain >= 0
        coords[active[use]] = x[use]
        error[active[use]] = curr[use]
        small = np.abs(gain) < 1e-10
        done = small & previous_small[active]
        previous_small[active] = small
        stopped[active[done]] = epoch
        converged[active[done]] = True
        active = active[~done]
    return coords, initial, np.sqrt(error), stopped, converged


def csv_out(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(out):
    init = time.perf_counter(), time.process_time()
    files = members(BENCH/"cubic_grid/human_panel/affinity-pools.tar.xz",
                    ["affinity/5.pdbqt", "affinity/6.pdbqt"])
    a, b = [Ligand(files[f"affinity/{i}.pdbqt"]) for i in (5,6)]
    architecture = tree(files["affinity/5.pdbqt"], a)
    matches = np.array(b.graph.GetSubstructMatches(a.graph, useChirality=True, uniquify=False))
    maps = {rec: Fields(rec) for rec in REC}
    target_xyz, keys, five_modes, six_modes = [], [], [], []
    for panel, arm, seed, pools in active_pools():
        for receiver in REC:
            for rank, (five, six) in enumerate(zip(modes(pools[receiver,5]), modes(pools[receiver,6]), strict=True), 1):
                sites = b.coords(six)
                target_xyz.extend(sites[matches])
                keys.append(dict(panel=panel, arm=arm, seed=seed, receiver=receiver, rank=rank))
                five_modes.append(a.coords(five))
                six_modes.append(sites)
    targets = np.array(target_xyz)
    starts = ("rigid", "tree")
    # Two self-geometry controls: the input itself and an actual native 5 pose.
    direct = np.stack([a.reference, five_modes[0]])
    from_input = solve(a.reference, direct, architecture, "tree")[2]
    if from_input[0] > 1e-10 or from_input[1] > 0.005:
        raise ArithmeticError(f"Failed self-geometry controls: {from_input}")

    candidates, initial, distance, steps, converged = zip(
        *(solve(a.reference, targets, architecture, start) for start in starts), strict=True)
    n, m = len(keys), len(matches)
    shape = (2,n,m)
    distance = np.array(distance).reshape(shape)
    candidates = np.array(candidates).reshape(*shape, len(a.order), 3)
    initial, steps, converged = [np.array(x).reshape(shape) for x in (initial, steps, converged)]
    trials, rows, chosen = [], [], []
    bonds = a.bonds
    # All distances within every rigid group plus flanking branch endpoints must be fixed.
    constraints = set()
    for grp in architecture[2]:
        constraints.update(tuple(sorted(p)) for p in itertools.combinations(grp,2))
    iedge, jedge = np.array(sorted(constraints)).T
    edge_ref = np.linalg.norm(a.reference[iedge]-a.reference[jedge], axis=1)
    rigid_error = 0.
    for i, key in enumerate(keys):
        for mi, si in itertools.product(range(m), range(2)):
            trials.append(dict(**key, mapping=mi, start=starts[si], initial_rmsd=initial[si,i,mi],
                               feasible_rmsd=distance[si,i,mi], sweeps=steps[si,i,mi],
                               converged=int(converged[si,i,mi])))
        mi, si = min(itertools.product(range(m),range(2)), key=lambda k: distance[k[1],i,k[0]])
        xyz = candidates[si,i,mi]
        chosen.append(xyz)
        rigid_error = max(rigid_error,np.max(np.abs(np.linalg.norm(xyz[iedge]-xyz[jedge], axis=1)-edge_ref)))
        a.check_center(xyz)
        official = key["arm"].startswith("official")
        score = maps[key["receiver"]].score
        s5 = score(a, five_modes[i], official)
        s6 = score(b, six_modes[i], official)
        if key["rank"] == 1:
            first5 = s5
        raw = targets[i*m+mi]
        raw5 = score(a,raw,official)
        try:
            legal5 = score(a,xyz,official)
        except ValueError as exc:
            if "outside map" not in str(exc):
                raise
            legal5 = float("nan")
        legal_sub, raw_sub = s6-legal5, s6-raw5
        legal_place, raw_place = legal5-first5, raw5-first5
        net = s6-first5
        if np.isfinite(legal5) and abs(legal_sub+legal_place-net) > 1e-9:
            raise ArithmeticError("Score identity failed")
        lower = []
        for match in matches:
            target=six_modes[i][match]
            def length(x):
                return np.linalg.norm(x[bonds[:,0]]-x[bonds[:,1]],axis=1)
            lower.append(max(abs(length(target)-length(a.reference)))/np.sqrt(2*len(a.order)))
        rows.append(dict(**key, mapping=mi, start=starts[si], sweeps=steps[si,i,mi],
                         converged=int(converged[si,i,mi]), bond_lower_bound=min(lower),
                         feasible_rmsd=distance[si,i,mi], s5_native=s5, s6_native=s6, s5_first=first5,
                         s5_deleted=raw5, s5_feasible=legal5, score_shift=legal5-raw5,
                         substitution_deleted=raw_sub, substitution_feasible=legal_sub,
                         placement_deleted=raw_place, placement_feasible=legal_place, net=net))
    if rigid_error>1e-8:
        raise ArithmeticError(f"Rigid fragment displacement {rigid_error}")
    out.mkdir(parents=True, exist_ok=True)
    csv_out(out/"manifold_poses.csv",rows)
    csv_out(out/"manifold_starts.csv",trials)
    np.savez_compressed(out/"manifold_xyz.npz", xyz=np.array(chosen),
                        pdbqt_serials=np.array(a.order), mappings=matches)
    meta = dict(poses=n, branches=a.tors, candidate_fits=len(trials), converged_starts=int(converged.sum()),
                native_self_control_rmsd=from_input.tolist(),
                max_rigid_group_distance_change=rigid_error,
                finite_scores=sum(np.isfinite(r["s5_feasible"]) for r in rows).item(),
                wall_s=time.perf_counter()-init[0],cpu_s=time.process_time()-init[1])
    (out/"manifold_verification.json").write_text(json.dumps(meta, indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    run(parser.parse_args().out)
