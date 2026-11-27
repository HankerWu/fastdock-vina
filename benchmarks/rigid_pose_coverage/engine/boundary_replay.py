"""Run the registered boundary sensitivity (prediction and post-run evaluation separately)."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from coverage_gate import ROOT, candidates, select

HERE = Path(__file__).resolve().parent
PANEL = HERE.parent / "pose"
ARMS = ("feasible_edge", "coverage", "uniform")
CASE = "3ANR_96107"


def horn_mask_gate(pocket, xyz, serials, ligand, bonds, seed, arm):
    # Both masks are computed on the same actual source. Ship only original SVD transforms.
    sys.path.insert(0, str(ROOT/"benchmarks/fragment_anchor/engine"))
    from verify_banks import horn
    pool = candidates(pocket, xyz, serials)
    anchor = xyz[pool["ligand_triangle"]]
    dest = pocket.sites[pocket.triangles]
    rotation = horn(anchor, dest)
    placed = (xyz-anchor.mean(0)) @ rotation.transpose(0,2,1) + dest.mean(1)[:,None,:]
    svd = pool["eligible"]
    alternate = ((placed >= pocket.center-pocket.half_width) & (placed <= pocket.center+pocket.half_width)).all(axis=(1,2))
    pool["eligible"] = alternate
    ids, _ = select(pool, xyz, seed, arm)
    centers, rots = pool["center"][ids], pool["matched_rotation"][ids]
    poses = (xyz-xyz.mean(0)) @ rots.transpose(0,2,1) + centers[:,None,:]
    shape = pocket._score(poses)
    return dict(triangle_id=ids, center=centers, matched_rotation=rots, matched_score=shape,
                edge_residual_A2=pool["edge_residual_A2"][ids], changed=np.where(svd!=alternate)[0],
                svd_count=svd.sum(), horn_count=alternate.sum())


def predict(work, engine, out):
    sys.path.insert(0, str(PANEL))
    import run_panel as panel
    cfg = panel.old.configs(work)["3ANR"]
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for arm in ARMS:
        dump = out/arm
        product = out/(arm+".pdbqt")
        native = panel.parent.native(cfg, 96107, product)
        env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        cmd = [sys.executable, __file__, "live", "--engine", engine, "--arm", arm, "--dump", dump, "--", *native]
        begin, cpu = time.perf_counter(), resource.getrusage(resource.RUSAGE_CHILDREN)
        with (out/(arm+".log")).open("w") as f:
            subprocess.run(list(map(str,cmd)), stdout=f, stderr=subprocess.STDOUT, check=True, env=env)
        end = resource.getrusage(resource.RUSAGE_CHILDREN)
        rows.append(dict(arm=arm, wall_s=time.perf_counter()-begin,
                         cpu_s=end.ru_utime+end.ru_stime-cpu.ru_utime-cpu.ru_stime))
    (out/"timing.json").write_text(json.dumps(rows, indent=2)+"\n")


def evaluate(work, out, dest):
    sys.path.insert(0, str(PANEL))
    import evaluate as panel_eval
    g=panel_eval.g
    cell,*_ = g.reference_inputs(work)["3ANR"]
    rows = []
    timing = {r["arm"]:r for r in json.loads((out/"timing.json").read_text())}
    for arm in ARMS:
        primary = work/"dump"/(CASE+"_"+arm)
        alt = out/arm
        assert (primary/"source.xyz").read_bytes() == (alt/"source.xyz").read_bytes()
        assert (primary/"source.q").read_bytes() == (alt/"source.q").read_bytes()
        p_original = work/"raw"/(CASE+"_"+arm+".pdbqt")
        p_alternate = out/(arm+".pdbqt")
        with np.load(primary/"bank.npz") as a, np.load(alt/"bank.npz") as b:
            ids, aa, bb = np.intersect1d(a["triangle_id"], b["triangle_id"], return_indices=True)
            for key in ("center", "matched_rotation", "edge_residual_A2", "matched_score"):
                np.testing.assert_array_equal(a[key][aa], b[key][bb])
            overlap, changes = len(ids), " ".join(map(str, b["changed"]))
        for tag, dump, product in (("svd",primary,p_original),("horn_mask",alt,p_alternate)):
            field = next(l.split() for l in (dump/"inject.log").read_text().splitlines() if l.startswith("field_rank "))
            winner=int(field[3])
            scores=np.genfromtxt(dump/"scores.csv", names=True,delimiter=",")["score"]
            assert winner == int(np.argmin(scores))
            with np.load(dump/"bank.npz") as bank:
                tri=int(bank["triangle_id"][winner])
            lines=(dump/"stages.xyz").read_text().splitlines()
            order=list(map(int,lines[0].split()[1:]))
            ix=[order.index(s) for s in cell["serials"]]
            metrics={}
            for stage in ("selected","hunt","second"):
                entry=next(words for l in lines if (words:=l.split())[:3]==[stage,"1","1"])
                coords=np.array(entry[3:],float).reshape(-1,3)[ix]
                metrics[stage+"_A"]=float(g.errors(coords,cell["refs"]))
            es=[float(g.errors(g.xyz(block,cell["serials"]),cell["refs"])) for block in product.read_text().split("ENDMDL")[:-1]]
            rows.append(dict(arm=arm, mask=tag, changed_ids=changes, overlap=overlap, winner_index=winner,
                             triangle_id=tri, coarse_score=scores[winner], **metrics, released=len(es),
                             top1_A=es[0], oracle_A=min(es), released_sha256=hashlib.sha256(product.read_bytes()).hexdigest(),
                             byte_equal=(p_original.read_bytes()==p_alternate.read_bytes()),
                             wall_s=timing[arm]["wall_s"] if tag=="horn_mask" else "",
                             cpu_s=timing[arm]["cpu_s"] if tag=="horn_mask" else ""))
    with dest.open("w", newline="") as f:
        writer=csv.DictWriter(f,rows[0].keys(),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    if sys.argv[1:2] == ["live"]:
        sys.argv.pop(1)
        from live import field_live
        field_live.main(generator=horn_mask_gate, arms=ARMS)
    else:
        p=argparse.ArgumentParser(description=__doc__)
        s=p.add_subparsers(dest="mode",required=True)
        run=s.add_parser("predict")
        audit=s.add_parser("evaluate")
        for cmd in (run, audit):
            cmd.add_argument("work",type=lambda p:Path(p).resolve())
            cmd.add_argument("out",type=lambda p:Path(p).resolve())
        run.add_argument("engine",type=lambda p:Path(p).resolve())
        audit.add_argument("dest",type=Path)
        args=p.parse_args()
        if args.mode=="predict":
            predict(args.work,args.engine,args.out)
        else:
            evaluate(args.work,args.out,args.dest)
