"""Reconstruct all retained pose candidates from rigid PDBQT bodies and Horn quaternions."""
import argparse
from collections import defaultdict
import csv
import io
import json
from itertools import combinations
from pathlib import Path
import tarfile
import time

import numpy as np
from scipy.spatial.transform import Rotation

from anchors import ROOT
from pocket_triangle import Pocket


def pdbqt_bodies(path):
    # The ordinary PDBQT nested torsion records give an independent check on cutting a bond graph.
    bodies, active = defaultdict(list), []
    for line in path.read_text().splitlines():
        words = line.split()
        if not words:
            continue
        if words[0] == "BRANCH":
            active.append(tuple(map(int, words[1:])))
        elif words[0] == "ENDBRANCH":
            assert active.pop() == tuple(map(int, words[1:]))
        elif words[0] in ("ATOM", "HETATM") and words[-1] not in ("H", "HD", "HS"):
            bodies[tuple(active)].append(int(line[6:11]))
    assert not active
    return sorted((sorted(v) for v in bodies.values()), key=lambda v: (-len(v), v))


def independent_anchor(x, rows, bodies):
    for component in bodies:
        if len(component) < 3:
            continue
        pairs = list(combinations(component, 2))
        d = [np.linalg.norm(x[rows[b]] - x[rows[a]]) for a, b in pairs]
        i = max(range(len(d)), key=lambda i: round(d[i], 10))
        if d[i] < 1e-8:
            continue
        a, b = pairs[i]
        unit = (x[rows[b]] - x[rows[a]]) / d[i]
        offsets = x[[rows[n] for n in component]] - x[rows[a]]
        h = np.linalg.norm(offsets - (offsets @ unit)[:, None] * unit, axis=1)
        k = max(range(len(h)), key=lambda j: round(h[j], 10))
        if h[k]*d[i] >= 1e-8:
            return component, [a, b, component[k]]
    return None, None


def horn(x, target):
    centered = x - x.mean(0)
    d = target - target.mean(axis=1)[:, None]
    s = np.einsum("ij,nik->njk", centered, d)
    trace = np.trace(s, axis1=1, axis2=2)
    k = np.zeros((len(d), 4, 4))
    k[:, 0, 0] = trace
    v = np.array([s[:,1,2] - s[:,2,1], s[:,2,0] - s[:,0,2], s[:,0,1] - s[:,1,0]]).T
    k[:, 0, 1:] = k[:, 1:, 0] = v
    k[:, 1:, 1:] = s + s.transpose(0,2,1) - trace[:,None,None] * np.eye(3)
    q = np.linalg.eigh(k)[1][:, :, -1]
    return Rotation.from_quat(q[:, [1,2,3,0]]).as_matrix()


def verify(work, target):
    start, cpu = time.perf_counter(), time.process_time()
    cases = json.loads((work/"inputs.json").read_text())["cells"]
    output = []
    archive = ROOT/"benchmarks/cavity_field_rank/pose/raw.tar.xz"
    with tarfile.open(archive) as tar:
        for path in sorted((work/"dump").glob("*/bank.npz")):
            name, seed, arm = path.parent.name.rsplit("_", 2)
            c = cases[name]
            srcrows = (path.parent/"source.xyz").read_text().splitlines()
            serials = np.asarray(srcrows[0].split()[1:], int)
            order = np.argsort(serials)
            serials = serials[order]
            x = np.asarray(srcrows[1].split()[3:], float).reshape(-1,3)[order]
            rows = {s:i for i,s in enumerate(serials)}
            bodies = pdbqt_bodies(work/c["ligand"])
            if arm == "whole":
                bodies = [list(serials)]
            component, anchor = independent_anchor(x, rows, bodies)
            fallback = component is None
            if fallback:
                component, anchor = independent_anchor(x, rows, [list(serials)])
            pocket = Pocket.from_pdbqt(work/c["receptor"], c["center"], c["size"])
            with np.load(path) as bank:
                np.testing.assert_array_equal(component, bank["component_serials"])
                np.testing.assert_array_equal(anchor, bank["anchor_serials"])
                assert bool(bank["anchor_fallback"]) == fallback
                ligand = x[[rows[s] for s in anchor]]
                corners = pocket.sites[pocket.triangles]
                residual = sum((np.linalg.norm(ligand[a] - ligand[b]) -
                                np.linalg.norm(corners[:,a] - corners[:,b],axis=1))**2
                               for a,b in ((0,1),(0,2),(1,2)))
                ids = np.argsort(np.round(residual,10),kind="stable")[:1024]
                np.testing.assert_array_equal(ids,bank["triangle_id"])
                np.testing.assert_array_equal(pocket.triangles[ids], bank["site_triangles"])
                target_xyz = corners[ids]
                fit = horn(ligand, target_xyz)
                center = (fit @ (x.mean(0) - ligand.mean(0))) + target_xyz.mean(axis=1)
                if arm == "sham":
                    shift = int(np.random.default_rng(int(seed)).integers(1,1024))
                    assert shift == int(bank["shift"])
                    fit = np.roll(fit, shift, axis=0)
                expected = (x - x.mean(0)) @ fit.transpose(0,2,1) + center[:,None]
                actual = (x-x.mean(0)) @ bank["matched_rotation"].transpose(0,2,1) + bank["center"][:,None]
                err = float(np.max(np.abs(expected - actual)))
                angle_error = float(np.max(np.abs(fit-bank["matched_rotation"])))
                assert err < 1e-8 and angle_error < 1e-9
                prior = f"banks/{name}_{seed}.npz"
                checked = 0
                if arm == "whole" and prior in tar.getnames():
                    with np.load(io.BytesIO(tar.extractfile(prior).read())) as archived:
                        for key in archived:
                            np.testing.assert_array_equal(archived[key], bank[key])
                            checked += 1
                scores = np.genfromtxt(path.with_name("scores.csv"),delimiter=",",names=True)
                assert np.array_equal(scores["index"],np.arange(1024))
                line = next(l.split() for l in path.with_name("inject.log").read_text().splitlines()
                            if l.startswith("field_rank"))
                assert int(line[3]) == scores["score"].argmin()
                output.append(dict(cell=path.parent.name,first_component=" ".join(map(str,component)),
                                   anchor=" ".join(map(str,anchor)),fallback=int(fallback),
                                   candidates=1024,max_matrix_abs=angle_error,max_coord_abs_A=err,
                                   old_whole_arrays_equal=checked,selected=int(line[3])))
    if len(output) != 48 or sum(row["old_whole_arrays_equal"] > 0 for row in output) != 10:
        raise ValueError("expected 16 cells, three arms each and ten frozen whole banks")
    with target.open("w",newline="") as f:
        w = csv.DictWriter(f,output[0].keys())
        w.writeheader()
        w.writerows(output)
    print("verified",len(output),"banks; wall",time.perf_counter()-start,"cpu",time.process_time()-cpu)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("csv",type=Path)
    args = p.parse_args()
    verify(args.work,args.csv)
