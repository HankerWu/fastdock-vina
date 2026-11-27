"""Synthetic direct-coordinate checks and regression against all ancestral whole-ligand banks."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation

from coverage_gate import _farthest, embedding, place, select
from pocket_triangle import Pocket


def check(work, dest):
    start, cpu = time.perf_counter(), time.process_time()
    rng = np.random.default_rng(13199)
    rotations = Rotation.random(40, rng=rng).as_matrix()
    centers = rng.normal(size=(40,3)) * 9
    errors = {}
    gen = rng.normal(size=(27,3)) * [15, 3, 0.3]
    for label, xyz in (("full",gen), ("planar",gen*[1,1,0]), ("linear",gen*[1,0,0]), ("point",gen*0)):
        explicit = (xyz-xyz.mean(axis=0)) @ rotations.transpose(0,2,1) + centers[:,None,:]
        sq_direct = pdist(explicit.reshape(len(centers), -1), metric="sqeuclidean")/len(xyz)
        sq_embed = pdist(embedding(xyz, centers, rotations), metric="sqeuclidean")
        errors[label] = float(np.abs(sq_direct-sq_embed).max())
        np.testing.assert_allclose(sq_direct, sq_embed, rtol=1e-13, atol=1e-12)
    tiny = np.array([[0.,0], [-1,0], [1,0], [1,0]])
    chosen, rad = _farthest(tiny, 0, 4)
    # stable tie by lowest ID and zero-distance duplicate consumes a distinct slot
    np.testing.assert_array_equal(chosen, [0,1,2,3])
    assert rad == 0.
    try:
        select(dict(eligible=np.array([True, False])), gen, 0, "uniform", 2)
    except ValueError as err:
        assert str(err) == "uniform: 1 available correspondences, require 2"
    else:
        raise AssertionError("expected empty-domain failure")
    cells = json.loads((work/"inputs.json").read_text())["cells"]
    old_keys, tested = set(), []
    for file in sorted((work/"dump").glob("*_whole/bank.npz")):
        name, seed, _ = file.parent.name.split("_")
        cfg = cells[name]
        for path, sha in cfg["sha256"].items():
            assert hashlib.sha256((work/path).read_bytes()).hexdigest() == sha
        header, row = (file.parent/"source.xyz").read_text().splitlines()
        s = np.array(header.split()[1:], int)
        x = np.array(row.split()[3:], float).reshape(-1,3)[np.argsort(s)]
        pocket = Pocket.from_pdbqt(work/cfg["receptor"], cfg["center"], cfg["size"])
        bank = place(pocket, x, np.sort(s), work/cfg["ligand"], work/cfg["bonds"], int(seed), "edge")
        with np.load(file) as old:
            for key in old:
                np.testing.assert_array_equal(old[key], bank[key], err_msg=f"{file.parent.name}:{key}")
                old_keys.add(key)
        tested.append(file.parent.name)
    assert len(tested) == 16
    dest.write_text(json.dumps(dict(max_pair_distance_error_A2=errors, pairs_per_geometry=len(sq_embed),
        sources=tested, equal_arrays=sorted(old_keys), ties_and_short_pool="pass",
        wall_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu), indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("old_pose_work", type=Path)
    p.add_argument("json", type=Path)
    args = p.parse_args()
    check(args.old_pose_work, args.json)
