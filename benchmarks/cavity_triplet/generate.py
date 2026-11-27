"""Generate and freeze both blind placement arms before opening any reference pose."""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from pocket_triangle import Pocket


def write(path, rows):
    with path.open("w", newline="") as f:
        out = csv.DictWriter(f, rows[0].keys())
        out.writeheader()
        out.writerows(rows)


def candidate_coords(x, result, arm):
    return (x - x.mean(axis=0)) @ result[arm + "_rotation"].transpose(0, 2, 1) + result["center"][:, None, :]


def profiles(table, out):
    out.mkdir(exist_ok=True, parents=True)
    cache, sites, costs, selections = {}, [], [], []
    for row in csv.DictReader(table.open()):
        receptor = row["receptor"]
        center = tuple(float(row["center_" + a]) for a in "xyz")
        size = float(row["size"])
        key = receptor, center, size
        if key not in cache:
            t0, c0 = time.perf_counter(), time.process_time()
            pocket = Pocket.from_pdbqt(receptor, center, size)
            index = len(cache)
            cache[key] = index, pocket
            sites.append(dict(receptor_index=index, receptor=Path(receptor).name,
                              wall_s=time.perf_counter()-t0, cpu_s=time.process_time()-c0,
                              nearest_queries=pocket.n_grid, burial_queries=pocket.n_burial,
                              retained_sites=len(pocket.sites), triangles=len(pocket.triangles)))
            np.save(out / f"sites{index}.npy", pocket.sites)
        index, pocket = cache[key]
        source = row["source"]
        with np.load(row["npz"]) as data:
            x, q = data["xyz"], data["q"]
        timings = {}
        res = pocket.place(x, int(row["seed"]), timings)
        np.savez_compressed(out / (source + ".npz"), **res)
        for arm in ("matched", "sham"):
            i = res[arm + "_best"]
            R = Rotation.from_matrix(res[arm + "_rotation"][i])
            T = R.apply(q[:3] - x.mean(axis=0)) + res["center"][i]
            quat = R * Rotation.from_quat(q[[4, 5, 6, 3]])
            selections.append(dict(source=source, arm=arm, index=i,
                **dict(zip(("qx", "qy", "qz", "qw", "qi", "qj", "qk"),
                           np.r_[T, quat.as_quat()[[3, 0, 1, 2]]], strict=True))))
        costs.append(dict(source=source, receptor_index=index, heavy=len(x),
                          n=len(res["center"]), triangles_tested=len(pocket.triangles),
                          unique_unordered_triangles=len(np.unique(np.sort(res["site_triangles"],axis=1),axis=0)),
                          unique_centers=len(np.unique(np.round(res["center"],8),axis=0)),
                          unique_rotations=len(np.unique(np.round(res["matched_rotation"],8),axis=0)),
                          field_queries_per_arm=len(x)*len(res["center"]), **timings))
    write(out.parent / "selections.csv", selections)
    write(out / "generation_cost.csv", costs)
    write(out / "site_cost.csv", sites)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sources",type=Path)
    p.add_argument("out",type=Path)
    a = p.parse_args()
    profiles(a.sources,a.out)
