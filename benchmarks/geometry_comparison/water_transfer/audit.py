"""Read the raw BRDT pose files and reconstruct assay order and receptor-frame RMSD."""
import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation


def coordinates(line):
    return [float(line[i:i+8]) for i in (30,38,46)]


def protein_alignment(query, receiver, data):
    if query == receiver:
        return np.array(data[query]["mapped"])
    left, right = data[query], data[receiver]
    keys = sorted(k for k in left["ca"].keys() & right["ca"].keys() if all(
        math.dist(entry["ca"][k], entry["center"]) < 14 for entry in (left, right)))
    x, y = (np.asarray([entry["ca"][k] for k in keys]) for entry in (left, right))
    rotation, _ = Rotation.align_vectors(y - y.mean(0), x - x.mean(0))
    src = np.asarray(left["mapped"])
    return rotation.apply((src - x.mean(0)).reshape(-1, 3)).reshape(src.shape) + y.mean(0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path, help="Directory containing cells.csv, modes.csv, inputs and arm subdirectories")
    args = p.parse_args()
    work = args.work
    rows = list(csv.DictReader((work/"cells.csv").open()))
    modes = list(csv.DictReader((work/"modes.csv").open()))
    data = {}
    for code in ("7MRH","7MRG"):
        entry = work/"inputs"/code
        data[code] = json.loads((entry/"metadata.json").read_text())
        data[code]["ca"] = {(l[22:27], l[17:20]):coordinates(l) for l in (entry/(code+".pdb")).read_text().splitlines()
            if l.startswith("ATOM  ") and l[21] == "A" and l[12:16].strip() == "CA" and l[16] in (" ","A")}
    groups = defaultdict(list)
    for r in modes:
        groups[r["query"],r["receiver"],r["arm"],r["seed"]].append(r)
    outcomes = defaultdict(lambda: defaultdict(list))
    energy = defaultdict(list)
    max_error, count = 0.0, 0
    for r in rows:
        q,rec,arm,seed = r["query"],r["receiver"],r["arm"],r["seed"]
        ref = protein_alignment(q,rec,data)
        raw = (work/arm/"poses"/f"{q}-{rec}-{seed}-{arm}.pdbqt").read_bytes()
        if hashlib.sha256(raw).hexdigest() != r["pdbqt_sha256"]:
            raise ValueError("Pool hash mismatch")
        raw_scores, rmsds = [], []
        for model in raw.decode().split("ENDMDL"):
            lines = model.splitlines()
            atoms = {int(l[6:11]):coordinates(l) for l in lines if l.startswith(("ATOM  ","HETATM"))}
            if not atoms:
                continue
            ordered = [atoms[i] for i in data[q]["serial"]]
            rmsds.append(min(math.sqrt(math.fsum(math.dist(a,b)**2 for a,b in zip(ordered,c)) / len(c)) for c in ref))
            raw_scores.append(float(next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3]))
        reported = groups[q,rec,arm,seed]
        assert len(reported) == len(rmsds) == int(r["modes"])
        for score, rmsd, mode in zip(raw_scores,rmsds,reported):
            assert score == float(mode["score"])
            max_error = max(max_error,abs(float(mode["rmsd_A"])-rmsd))
        outcomes[arm][r["kind"]].append((rmsds[0], min(rmsds[:5]), min(rmsds)))
        energy[arm,q,rec].append(raw_scores[0])
        count += len(rmsds)
    assert max_error < 1e-9
    results = {}
    for arm in outcomes:
        groups = outcomes[arm]
        scores = {rec: {q:statistics.median(energy[arm,q,rec]) for q in data} for rec in data}
        results[arm] = dict(hits_top1_top5_any={kind:[sum(x[i]<=2 for x in cells) for i in range(3)]
                                             for kind,cells in groups.items()},
                            medians_by_receiver_query=scores,
                            correct_assay_order_by_receiver={rec:m["7MRH"] < m["7MRG"] for rec,m in scores.items()})
    print(json.dumps(dict(runs=len(rows),modes=count,max_abs_rmsd_difference_A=max_error,arms=results),indent=2))


if __name__=="__main__":
    main()
