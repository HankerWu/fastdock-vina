"""Independent Gemmi/RDKit crystal-frame checks on published witnesses, or every raw coordinate."""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit.Chem import rdMolAlign

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_predictions as independent


def verify(work, directory, raw=None):
    rows = list(csv.DictReader((directory/"runs.csv").open()))
    stages = list(csv.DictReader((directory/"stages.csv").open()))
    by_drug = {}
    totals = Counter()
    max_raw = max_witness = 0.0
    for row in rows:
        drug, graph, seed = (row[s] for s in ("drug", "graph", "seed"))
        case = "_".join((drug, graph, seed))
        crystal, anchors = independent.refmol(work, drug)
        mol, serials = independent.frame.read_query(work, drug)
        perms = crystal.GetSubstructMatches(mol, uniquify=False, maxMatches=100000)
        if not perms or len(perms) == 100000:
            raise ValueError(("Graph correspondence", drug, len(perms)))
        by_drug[drug] = {"anchor_residues": anchors, "equivalent_mappings": len(perms)}
        with np.load(directory/(case+".npz"), allow_pickle=False) as z:
            a = dict(z.items())
        reorder = [a["serials"].tolist().index(s) for s in serials]
        r = a["rmsd_A"]
        for i, xyz in zip(a["witness_index"], a["witness_xyz"], strict=True):
            for j, coord in enumerate(xyz[reorder]):
                mol.GetConformer().SetAtomPosition(j, coord)
            actual = rdMolAlign.CalcRMS(mol, crystal, symmetrizeConjugatedTerminalGroups=False)
            max_witness = max(max_witness, abs(r[i] - actual))
        for stage in (v for v in stages if (v["drug"],v["graph"],v["seed"]) == (drug,graph,seed)):
            i = int(stage["min_row"])
            if i not in a["witness_index"] or abs(r[i] - float(stage["min_A"])) > 1e-11:
                raise ValueError(("Minimum witness", case, stage))
            for level in (2, 3):
                j = stage[f"first_le{level}_row"]
                if j and (int(j) not in a["witness_index"] or r[int(j)] > level):
                    raise ValueError(("Entry witness", case, stage))
        totals["witness_coordinates"] += len(a["witness_xyz"])
        totals["events"] += len(r)
        if raw is None:
            continue
        binary = raw/(case+".bin")
        if independent.digest(binary) != row["raw_sha256"]:
            raise ValueError(f"Raw file digest disagrees for {case}")
        n = len(serials)
        dtype = np.dtype([("tags", "<u4", (5,)), ("E", "<f8"), ("x", "<f8", (n, 3))])
        trace = np.memmap(binary, dtype=dtype, mode="r", offset=4*(n+1))
        with binary.open("rb") as f:
            header = np.fromfile(f, "<u4", n+1)
        assert header.tolist() == [n, *a["serials"]]
        for k, field in enumerate(("traj", "step", "phase", "event", "local")):
            assert np.array_equal(trace["tags"][:, k], a[field])
        assert np.array_equal(trace["x"][a["witness_index"]], a["witness_xyz"])
        references = crystal.GetConformer().GetPositions()[list(perms)]
        for i in range(0, len(trace), 256):
            xyz = trace["x"][i:i+256][:, reorder, :]
            d = xyz[:, None, :, :] - references[None, :, :, :]
            val = np.sqrt(np.einsum("bmij,bmij->bm", d, d).min(axis=1) / n)
            expected = r[i:i+len(val)]
            max_raw = max(max_raw, float(np.max(np.abs(val-expected))))
            if not np.allclose(val, expected, rtol=5e-14, atol=1e-10):
                raise ValueError(("Independent complete-coordinate RMSD", case, i))
        # Validate every returned grid minimization against its actual first
        # coordinate and most recently adopted line-search point.
        end = np.flatnonzero(np.isin(a["event"], [4, 6]))
        start = np.flatnonzero((a["event"] == 2) & (a["local"] == 1) & np.isin(a["phase"], [1, 2]))
        for lo, hi in zip(start, end, strict=True):
            assert np.array_equal(trace["tags"][lo,:3], trace["tags"][hi,:3])
            steps = np.flatnonzero(a["event"][lo:hi] == 3) + lo
            expected = steps[-1] if len(steps) else lo
            close = lambda l,h: np.array_equal(trace["x"][l], trace["x"][h])
            if not close(hi, expected):
                if not close(hi, lo):
                    raise ValueError(("Returned coordinates differ from both eligible states", case, lo, hi))
                totals["returned_initial_over_last_adopted"] += 1
            if not close(hi, hi-1):
                totals["returned_different_from_last_evaluation"] += 1
            totals["grid_local_returns"] += 1
        totals["full_coordinate_events_checked"] += len(trace)
    if max_witness > 1e-9:
        raise ValueError(("Independent witness error", max_witness))
    return dict(cases=len(rows), models=by_drug, counts=dict(totals),
                max_witness_abs_diff_A=max_witness, max_full_coordinate_abs_diff_A=max_raw if raw else None)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("results", type=Path)
    p.add_argument("--raw", type=Path)
    args = p.parse_args()
    print(json.dumps(verify(args.work, args.results, args.raw), indent=2))
