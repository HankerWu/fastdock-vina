"""Verify five new gate banks directly with full atom coordinates and Horn's quaternions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"benchmarks/fragment_anchor/engine"))
from verify_banks import Pocket, horn, independent_anchor
from coverage_gate import candidates


def filehash(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def direct_gate(poses, feasible, first, budget):
    """All 3N components contribute directly, with no covariance or rotation embedding."""
    a = np.ascontiguousarray(poses[feasible].reshape(len(feasible), -1))
    selected, current = [], int(first)
    near = np.full(len(a), np.inf)
    for j in range(budget):
        selected.append(current)
        d2 = cdist(a, a[current:current+1], metric="sqeuclidean")[:, 0] / poses.shape[1]
        near = np.minimum(near, d2)
        near[current] = -np.inf
        if j+1 < budget:
            current = int(np.argmax(np.round(near, 10)))
    return np.sort(feasible[selected]), max(0., float(np.max(near)))


def verify(work, out):
    begin, cpu = time.perf_counter(), time.process_time()
    meta = json.loads((work/"inputs.json").read_text())
    for path, checksum in meta["archives"].items():
        assert filehash(ROOT / path) == checksum
    results, boundary = [], []
    source_files = sorted((work/"dump").glob("*_coverage/source.xyz"))
    if not source_files:
        raise ValueError("missing prospective dumps")
    for source_file in source_files:
        stem = source_file.parent.name.removesuffix("_coverage")
        name, seed = stem.rsplit("_", 1)
        cfg = meta["cells"][name]
        for path, checksum in cfg["sha256"].items():
            assert filehash(work/path) == checksum
        headers, coords = source_file.read_text().splitlines()
        src_ids = np.array(headers.split()[1:], dtype=int)
        serials = np.sort(src_ids)
        x = np.array(coords.split()[3:], float).reshape(-1, 3)[np.argsort(src_ids)]
        p = Pocket.from_pdbqt(work/cfg["receptor"], cfg["center"], cfg["size"])
        rows = {s:i for i,s in enumerate(serials)}
        _, anchors = independent_anchor(x, rows, [list(serials)])
        a = x[[rows[i] for i in anchors]]
        corners = p.sites[p.triangles]
        rot = horn(a, corners)
        dest = corners.mean(axis=1)
        # Independent affine convention: subtract anchor centroid before rotation.
        placed = (x-a.mean(axis=0)) @ rot.transpose(0, 2, 1) + dest[:, None, :]
        mask = ((placed >= p.center - p.half_width) & (placed <= p.center + p.half_width)).all(axis=(1,2))
        alternative = mask.copy()
        frozen = candidates(p, x, serials)
        mask = frozen["eligible"]
        # Cross-solver exact box membership can differ at a nominal box face.
        # Retain that observation, then use the frozen SVD support to verify
        # selection, using independent Horn/Cartesian positions for distances.
        native_xyz = (x-x.mean(axis=0)) @ frozen["matched_rotation"].transpose(0,2,1) + frozen["center"][:,None,:]
        all_error = float(np.abs(native_xyz - placed).max())
        assert all_error < 1e-8
        def excess(arr):
            return (np.abs(arr - p.center) - p.half_width).max(axis=(1,2))
        changed = np.where(mask != alternative)[0]
        for i, da, db in zip(changed, excess(native_xyz[changed]), excess(placed[changed]), strict=True):
            assert max(abs(da), abs(db)) < 1e-9
            boundary.append(dict(cell=stem, triangle_id=i, svd=int(mask[i]), horn=int(alternative[i]),
                                 svd_signed_wall_A=da, horn_signed_wall_A=db))
        feasible = np.flatnonzero(mask)
        edge = sum((np.linalg.norm(a[i]-a[j]) - np.linalg.norm(corners[:,i]-corners[:,j], axis=1))**2
                   for i,j in ((0,1), (0,2), (1,2)))
        edge_order = np.argsort(np.round(edge, 10), kind="stable")
        eligible_order = edge_order[mask[edge_order]]
        cover, radius = direct_gate(placed, feasible, np.searchsorted(feasible,eligible_order[0]), 1024)
        uniform = np.sort(np.random.Generator(np.random.PCG64(int(seed))).choice(feasible, size=1024, replace=False))
        chosen = dict(edge=edge_order[:1024], feasible_edge=eligible_order[:1024], coverage=cover,
                      uniform=uniform, uniform_work=uniform)
        alt_ids = np.flatnonzero(alternative)
        alt_order = edge_order[alternative[edge_order]]
        if len(changed):
            alt_cover, _ = direct_gate(placed, alt_ids, np.searchsorted(alt_ids, alt_order[0]), 1024)
            alt_uniform = np.sort(np.random.Generator(np.random.PCG64(int(seed))).choice(alt_ids, size=1024, replace=False))
        else:
            alt_cover, alt_uniform = cover, uniform
        horn_selections = dict(edge=edge_order[:1024], feasible_edge=alt_order[:1024],
                               coverage=alt_cover, uniform=alt_uniform, uniform_work=alt_uniform)
        for arm, expect in chosen.items():
            dump = work/"dump"/(stem+"_"+arm)
            assert source_file.read_bytes() == (dump/"source.xyz").read_bytes()
            with np.load(dump/"bank.npz") as b:
                np.testing.assert_array_equal(expect, b["triangle_id"], err_msg=f"{stem} {arm} feasible={len(feasible)} native={b['feasible_count']}")
                np.testing.assert_array_equal(p.triangles[expect], b["site_triangles"])
                np.testing.assert_array_equal(anchors, b["anchor_serials"])
                if arm != "edge":
                    assert int(b["feasible_count"]) == len(feasible)
                error = float(np.abs(placed[expect] - ((x-x.mean(0)) @ b["matched_rotation"].transpose(0,2,1)
                                                        + b["center"][:,None])).max())
                assert error < 1e-9
                if arm in ("coverage", "uniform_work"):
                    np.testing.assert_allclose(radius, b["radius_A2"], atol=1e-9, rtol=0)
            scores = np.genfromtxt(dump/"scores.csv", names=True, delimiter=",")
            field = next(l.split() for l in (dump/"inject.log").read_text().splitlines() if l.startswith("field_rank"))
            assert int(field[3]) == int(np.argmin(scores["score"]))
            assert len(scores) == 1024 and np.array_equal(scores["index"], np.arange(1024))
            if arm == "uniform_work":
                assert (work/"raw"/(stem+"_uniform_work.pdbqt")).read_bytes() == (work/"raw"/(stem+"_uniform.pdbqt")).read_bytes()
            results.append(dict(cell=stem, arm=arm, source_sha256=filehash(source_file), all=len(placed),
                                feasible=len(feasible), horn_feasible=int(alternative.sum()), boundary_flips=len(changed),
                                all_max_coordinate_error_A=all_error, selected=len(expect),
                                both_supports_selected=len(np.intersect1d(expect, horn_selections[arm])), max_coordinate_error_A=error,
                                explicit_radius_A2=radius if arm=="coverage" else "", winner=int(field[3]),
                                winner_ties=int(np.sum(scores["score"] == scores["score"].min())),
                                winner_original_triangle=int(expect[int(field[3])])))
        print(stem, len(feasible), flush=True)
    with (out/"coordinate_check.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, results[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    with (out/"boundary_flips.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, ["cell","triangle_id","svd","horn","svd_signed_wall_A","horn_signed_wall_A"],
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(boundary)
    (out/"coordinate_check_cost.json").write_text(json.dumps(dict(cells=len(source_files), banks=len(results),
        wall_s=time.perf_counter()-begin, cpu_s=time.process_time()-cpu), indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    verify(args.work, args.output)
