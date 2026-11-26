"""Replay normal-size boxes with the projected explicit gradient and compare full pools."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from run import experiment, measure, reference, write_csv


def rows(fixtures, original, out, binary, workers):
    cells = [r for r in csv.DictReader((original/"confirmation.csv").open())
             if r["arm"] in ("native","cubic")]
    prepared = {p.parent.name:json.loads(p.read_text()) for p in fixtures.glob("*/preparation.json")}
    jobs = [(fixtures,out,binary,prepared[r["query"]],prepared[r["receiver"]],r["arm"],
             int(r["seed"]),int(r["cap"])) for r in cells]
    with ProcessPoolExecutor(workers) as executor:
        for r,new in zip(cells, executor.map(experiment, jobs)):
            q,rec=prepared[r["query"]],prepared[r["receiver"]]
            name=r["pool"]
            oldpath, newpath = original/name, out/new["pool"]
            values=measure(newpath,q["serial"],reference(q,rec))
            yield dict(query=r["query"], receiver=r["receiver"], arm=r["arm"], seed=r["seed"],
                       original_sha256=hashlib.sha256(oldpath.read_bytes()).hexdigest(),
                       corrected_sha256=hashlib.sha256(newpath.read_bytes()).hexdigest(),
                       old_first=r["first"], new_first=values[0], old_best=r["best"], new_best=min(values))


if __name__ == "__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures","original","out","binary"):
        ap.add_argument(name,type=Path)
    ap.add_argument("--jobs",type=int,default=1)
    a=ap.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    data=list(rows(a.fixtures.resolve(),a.original.resolve(),a.out.resolve(),a.binary.resolve(),a.jobs))
    write_csv(a.out/"wall_replay.csv",data)
    print(sum(r["original_sha256"]==r["corrected_sha256"] for r in data),"identical complete pools of",len(data))
    print(sum((float(r["old_first"])<=2)!=(r["new_first"]<=2) for r in data),"top-1 hit status changes")
