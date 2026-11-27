"""Crystal-informed diagnostic on the four pre-existing BRD9 candidate banks."""
import argparse
from pathlib import Path

from pose_metrics import inspect_fork, rows, scan, table


def main(p):
    cell = scan.config(p/"aurora",p/"brd",p)["82I"]
    detail=[]
    for r in rows(p/"posthoc_selections.csv"):
        stem=r["source"]+"_"+r["arm"]
        _,_,chain=inspect_fork(p,stem,cell,1,1)
        assert abs(float(r["input_rmsd_A"])-chain["selected_A"])<1e-8
        pool=scan.metrics(cell,p/"raw"/(stem+".pdbqt"))
        near=[i for i,(rms,_) in enumerate(pool,1) if rms<=2]
        detail.append(dict(**chain,coarse_rank=r["coarse_rank"],modes=len(pool),
                           top1_A=pool[0][0],oracle_A=min(rms for rms,_ in pool),
                           first_le2_rank=near[0] if near else ""))
    table(p/"posthoc.csv",detail)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work",type=lambda v:Path(v).resolve())
    main(ap.parse_args().work)
