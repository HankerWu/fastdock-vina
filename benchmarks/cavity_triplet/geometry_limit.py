"""Post hoc exact bound on re-pairing the generated centers and orientations."""
import argparse
import csv
import json
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from pose_runs import OLD, table
sys.path.insert(0, str(OLD / "rigid_coupled"))
import study


def references(system, serials):
    mol, order, crystal = system[-3:]
    matches = crystal.GetSubstructMatches(mol, uniquify=False, maxMatches=100000)
    if not 0 < len(matches) < 100000:
        raise ValueError("Incomplete atom mapping")
    indices = np.array(matches)[:, [order.index(i) for i in serials]]
    return crystal.GetConformer().GetPositions()[indices]


def bound(x, centers, rotations, refs, shift):
    parent = x - x.mean(axis=0)
    ref_center = refs[0].mean(axis=0)
    target = refs - ref_center
    # MSE = squared centroid distance + centered MSE; permutations preserve centroid.
    translation = np.sum((centers - ref_center) ** 2, axis=1)
    rotated = parent @ rotations.transpose(0, 2, 1)
    orientation = np.minimum.reduce([
        ((rotated - ref) ** 2).sum(axis=2).mean(axis=1) for ref in target
    ])
    k = len(centers)
    index = np.arange(k)
    all_shifts = (index[None, :] - index[:, None]) % k
    err2 = translation[None, :] + orientation[all_shifts]
    min_by_shift = np.sqrt(err2.min(axis=1))
    # Also optimize over all orientations to measure the internal-shape floor.
    rigid = []
    for ref in target:
        u, _, vt = np.linalg.svd(parent.T @ ref)
        u[:, -1] *= np.linalg.det(u @ vt)
        rigid.append(np.linalg.norm(parent @ u @ vt - ref) / np.sqrt(len(x)))
    result = dict(rigid_shape_min_A=min(rigid), closest_centroid_A=np.sqrt(translation.min()),
                  centered_orientation_min_A=np.sqrt(orientation.min()),
                  best_any_pair_A=min_by_shift.min(), matched_best_A=min_by_shift[0],
                  shuffle_best_A=min_by_shift[shift], median_other_shift_A=np.median(min_by_shift[1:]),
                  better_or_equal_other_shifts=int(np.sum(min_by_shift[1:] <= min_by_shift[0]+1e-10)),
                  total_other_shifts=k-1,total_pairs=k*k)
    for threshold in (2, 3):
        near = err2 <= threshold**2
        result.update({f"pairs_within_{threshold}A":int(np.sum(near)),
                       f"shifts_within_{threshold}A":int(np.sum(near.any(axis=1))),
                       f"matched_within_{threshold}A":int(near[0].sum()),
                       f"shuffle_within_{threshold}A":int(near[shift].sum())})
    return result


def common_eviction(system):
    path = OLD/"rigid_capture/eviction.tar.xz"
    with tarfile.open(path) as archive:
        coords = archive.extractfile("stages.xyz").read().decode().splitlines()
        events = archive.extractfile("events.txt").read().decode().splitlines()
    refs = references(system, list(map(int, coords[0][1:].split())))
    offers = {int(line.split()[2]):np.array(line.split()[3:],float).reshape(-1,3)
              for line in coords[1:] if line.startswith("offered ")}
    errors = {i:np.sqrt(((x - refs)**2).sum(axis=2).mean(axis=1)).min() for i,x in offers.items()}
    energies, current = {}, []
    for line in events:
        step, energy, capacity, exclusion, *items = line.split()
        i, size, eps = int(step), int(capacity), float(exclusion)
        previous = current[:]
        energies[i] = float(energy)
        # Independent greedy replay on all residents and the current offer.
        current = []
        for candidate in sorted(previous + [i], key=energies.get):
            if all(np.linalg.norm(offers[candidate]-offers[j])/np.sqrt(len(offers[j])) >= eps
                   for j in current):
                current.append(candidate)
                if len(current) == size:
                    break
        assert len(current) == len(items)
        assert all(np.array_equal(offers[a],offers[int(b)]) for a,b in zip(current,items))
        if i == 137:
            lost = next(j for j in previous if errors[j] <= 2 and j not in current)
            witness = dict(step=i, lost_id=lost, lost_A=errors[lost], new_A=errors[i], size=size,
                           new_score=energies[i],lost_score=energies[lost],
                           direct_distance_A=np.linalg.norm(offers[i]-offers[lost])/np.sqrt(len(offers[i])),
                           exclusion_A=eps,remaining_best_A=min(errors[j] for j in current))
    return dict(offers=len(offers), final_best_A=min(errors[j] for j in current),witness=witness)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        systems=study.audit.inputs(args.work/"aurora",args.work/"brd",Path(tmp))
    rows=[]
    for event in csv.DictReader((args.work/"sources.csv").open()):
        key=event["source"]
        with np.load(args.work/"sources"/(key+".npz")) as source, np.load(args.work/"candidates"/(key+".npz")) as c:
            refs=references(systems[event["drug"]],source["serials"].tolist())
            rows.append(dict(source=key,**bound(source["xyz"],c["center"],c["matched_rotation"],refs,int(c["shift"]))))
    table(args.out/"geometry_limit.csv",rows)
    (args.out/"historical_eviction.json").write_text(json.dumps(common_eviction(systems["at9283"]),indent=2)+"\n")
