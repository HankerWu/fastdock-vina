"""Replay the archived BRD9 scalar/batched failure at a single identical selected ID."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
FIELD = ROOT/"benchmarks/cavity_field_rank"
sys.path.insert(0, str(ROOT/"benchmarks/cavity_triplet"))
from pose_runs import invoke


def replay(work, inputs, binary):
    work.mkdir(parents=True, exist_ok=True)
    (work/"raw").mkdir(exist_ok=True)
    c = dict(json.loads((inputs/"inputs.json").read_text())["cells"]["82I"])
    for key in ("ligand","bonds","receptor"):
        c[key] = inputs/c[key]
    costs = []
    with tarfile.open(FIELD/"engine/live_raw.tar.xz") as old, tarfile.open(FIELD/"pose/raw.tar.xz") as pose:
        vectors = json.load(old.extractfile("live_failure/vectors.json"))
        base = np.loadtxt(io.BytesIO(pose.extractfile("banks/82I_91813.bank").read()))
        for key, archive_dir in (("native_scalar", "live_pinned"), ("vectored", "live_failure")):
            reproduced = base.copy()
            reproduced[vectors["candidate"]] = vectors[key]
            path = work/(key+".bank")
            np.savetxt(path, reproduced, fmt="%.17g")
            event = work/"event"
            event.write_text("1 1 0 0 0 1 0 0 0\n")
            trace = work/(key+".inject")
            pool = invoke(work,costs,key,binary,c,91813,env={
                "FD_FIELD_BANK":path, "FD_FIELD_ARM":"true", "FD_RIGID_INTERVENTION":event,
                "FD_RIGID_LOG":trace})
            artifact = archive_dir+"/82I_91813_true.pdbqt"
            expected = old.extractfile(artifact).read()
            assert pool.read_bytes() == expected
            line = next(l.split() for l in trace.read_text().splitlines() if l.startswith("field_rank"))
            assert int(line[3]) == vectors["candidate"]
            costs[-1].update(archive=artifact,field_index=int(line[3]),
                             input_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (work/"roundoff.json").write_text(json.dumps(costs,indent=2)+"\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("work","inputs","binary"):
        p.add_argument(key,type=lambda value:Path(value).resolve())
    args = p.parse_args()
    replay(args.work,args.inputs,args.binary)
