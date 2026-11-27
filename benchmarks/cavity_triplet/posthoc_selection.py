"""Crystal-informed diagnostic only: rescue a known proximal candidate at fixed coordinates."""

import argparse
import csv
import gzip
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from generate import write


def main(w):
    with gzip.open(w/"candidate_analysis/candidate_metrics.csv.gz","rt") as f:
        best = {}
        for row in csv.DictReader(f):
            key = row["source"],row["arm"]
            if not key[0].startswith("82I_") or float(row["rmsd_A"]) > 2:
                continue
            if key not in best or int(row["coarse_rank"]) < int(best[key]["coarse_rank"]):
                best[key] = row
    output = []
    for (source,arm),row in sorted(best.items()):
        i = int(row["index"])
        with np.load(w/"candidates"/(source+".npz")) as res, np.load(w/"sources"/(source+".npz")) as data:
            R = Rotation.from_matrix(res[arm+"_rotation"][i])
            x,q = data["xyz"],data["q"]
            origin = R.apply(q[:3]-x.mean(axis=0))+res["center"][i]
            quat = R*Rotation.from_quat(q[[4,5,6,3]])
            output.append(dict(source=source,arm=arm+"_rescue",index=i,
                **dict(zip(("qx","qy","qz","qw","qi","qj","qk"),np.r_[origin,quat.as_quat()[[3,0,1,2]]],strict=True)),
                input_rmsd_A=row["rmsd_A"], coarse_rank=row["coarse_rank"]))
    assert len(output)==4
    write(w/"posthoc_selections.csv",output)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    main(p.parse_args().work)
