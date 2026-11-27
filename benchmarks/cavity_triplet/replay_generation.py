"""Measure the full two-arm prediction driver and require exact replay of frozen candidates."""

import argparse
import csv
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


def run(w, output):
    w = w.resolve()
    with tempfile.TemporaryDirectory() as tmp:
        replica = Path(tmp) / "candidates"
        t0 = time.perf_counter()
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        subprocess.run([sys.executable, str(Path(__file__).with_name("generate.py")), str(w/"sources.csv"), str(replica)],
                       env=os.environ | {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}, check=True)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        meter = dict(wall_s=time.perf_counter()-t0,
                     cpu_s=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
        assert (replica.parent/"selections.csv").read_bytes() == (w/"selections.csv").read_bytes()
        n = 0
        for row in csv.DictReader((w/"sources.csv").open()):
            source = row["source"]
            with np.load(replica/(source+".npz")) as regenerated, np.load(w/"candidates"/(source+".npz")) as frozen:
                assert regenerated.files==frozen.files
                for name in regenerated.files:
                    np.testing.assert_array_equal(regenerated[name],frozen[name])
                    n += 1
        for path in replica.glob("sites*.npy"):
            np.testing.assert_array_equal(np.load(path),np.load(w/"candidates"/path.name))
        output.write_text(json.dumps(dict(**meter, arrays_compared=n,
            numpy=np.__version__, same_selection_csv_bytes=True,
            scope="fresh Python process, imports, all ordinary inputs, shared receptor setup and both arms for 24 sources, file compression/writes; excludes MC and offline evaluation"),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    run(a.work,a.out)
