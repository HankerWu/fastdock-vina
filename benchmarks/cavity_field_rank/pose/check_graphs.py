"""Compare adjacency, neighbor order, XS typing and torsion flags from both Pro parser paths."""
import argparse
import csv
import io
from pathlib import Path
import subprocess

from run import configs, table, write_json

def check(p, probe, out):
    audit = {}
    allrows = []
    for name,c in configs(p).items():
        graphs = []
        for label in ("inferred","explicit"):
            cmd=[probe,c["ligand"]]+([c["bonds"]] if label=="explicit" else [])
            rows=list(csv.DictReader(io.StringIO(subprocess.check_output(cmd,text=True,stderr=subprocess.DEVNULL))))
            allrows.extend(dict(drug=name,graph=label,**r) for r in rows)
            graphs.append(rows)
        edge=lambda rows:{tuple(sorted((int(r["atom"]),int(r["neighbor"])))) for r in rows}
        a,b=map(edge,graphs)
        expected={tuple(sorted(map(int,l.split()))) for l in c["bonds"].read_text().splitlines()}
        assert expected==b
        audit[name]=dict(inferred_bonds=len(a),explicit_bonds=len(b),missing=sorted(b-a),
            extra=sorted(a-b),identical_adjacency_types_and_rotors=graphs[0]==graphs[1])
    table(out/"graphs.csv",allrows)
    write_json(out/"graphs.json",audit)

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("work","probe","output"): p.add_argument(k,type=Path)
    a=p.parse_args(); check(a.work.resolve(),a.probe.resolve(),a.output.resolve())
