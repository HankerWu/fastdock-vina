"""Read frozen one-pass pose banks, lineage and published pools after the blind runs."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path

import numpy as np


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = Path(__file__).resolve().parent
panel = load("fragment_pose_panel", HERE/"panel.py")
old = load("reference_evaluation", panel.RANK/"pose/evaluate.py")
old.panel.CASES = panel.CASES
old.COMPONENTS.update({"2NNQ":"T4B","3FR2":"8CA"})


def geometry(work):
    return old.reference_inputs(work)


def distribution(rows, arm, group, unit="instance"):
    valid = [r for r in rows if r["n"]]
    if unit == "molecule":
        groups = {r["group_molecule"] for r in valid}
        measures = [{key:np.mean([r[key] for r in valid if r["group_molecule"]==m])
                    for key in ("top1_A","best5_A","oracle_A","regret_A")} for m in sorted(groups)]
    else:
        measures = valid
    x = np.array([r["top1_A"] for r in measures])
    return dict(arm=arm,group=group,unit=unit,n=len(measures),cells=len(rows),empty=len(rows)-len(valid),
                **{k:np.mean([r[k] for r in measures]) if len(x) else "" for k in ("top1_A","best5_A","oracle_A","regret_A")},
                **dict(zip(("q25","median","q75","p90","p95"),np.percentile(x,[25,50,75,90,95]) if len(x) else [""]*5,strict=True)),
                gt5=sum(x>5),gt10=sum(x>10))


def run(p, out):
    out.mkdir(parents=True,exist_ok=True)
    geo = geometry(p)
    runs = {r["stem"]:r for fname in ("control_cost.csv","run_cost.csv") for r in panel.old.csvrows(p/fname)}
    cells,modes,chem,forks,stages,history,dist = [],[],[],[],[],[],[]
    with gzip.open(out/"candidates.csv.gz","wt",newline="") as stream:
        writer = None
        for name, (query,receiver,seeds) in panel.CASES.items():
            cell,mol,src,tree = geo[name]
            for seed in seeds:
                for arm in (*panel.ARMS,"pro","pro-high","official","official-high"):
                    stem = f"{name}_{seed}_{arm}"
                    raw = p/"raw"/stem
                    pool = raw.with_suffix(".pdbqt")
                    assert panel.digest(pool)==runs[stem]["pool_sha256"]
                    k = dict(stem=stem,partition=panel.partition(name,seed),query=query,receiver=receiver,
                             drug=name,group_molecule=query,seed=seed,arm=arm)
                    e = []
                    for rank,block in enumerate(pool.read_text().split("ENDMDL")[:-1],1):
                        x = old.xyz(block,cell["serials"])
                        d = float(old.errors(x,cell["refs"]))
                        score = float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
                        modes.append(dict(**k,rank=rank,rmsd_A=d,score=score))
                        chem.append(dict(**k,rank=rank,**old.quality(x,src,mol,tree)))
                        e.append(d)
                    cells.append(dict(**k,n=len(e),top1_A=e[0] if e else "",best5_A=min(e[:5],default=""),
                                      oracle_A=min(e,default=""),regret_A=e[0]-min(e) if e else ""))
                    if arm not in panel.ARMS:
                        continue
                    b = p/"dump"/stem
                    head, text = (b/"source.xyz").read_text().splitlines()
                    ids = list(map(int,head.split()[1:]))
                    x = np.fromstring(text.split(maxsplit=3)[3],sep=" ").reshape(len(ids),3)
                    order = [ids.index(n) for n in cell["serials"]]
                    with np.load(b/"bank.npz") as bank:
                        poses = (x-x.mean(0)) @ bank["matched_rotation"].transpose(0,2,1) + bank["center"][:,None,:]
                        ds = old.errors(poses[:,order],cell["refs"])
                        ts = np.array(bank["triangle_id"])
                        anchors = " ".join(map(str,bank["anchor_serials"]))
                        component = " ".join(map(str,bank["component_serials"]))
                        residuals = bank["edge_residual_A2"]
                    scores = panel.old.csvrows(b/"scores.csv")
                    assert [int(r["index"]) for r in scores]==list(range(1024))
                    vals = np.array([float(r["score"]) for r in scores])
                    indices = np.argsort(vals,kind="stable")
                    ranks = np.argsort(indices)+1
                    i = int(indices[0])
                    inj = next(l.split() for l in (b/"inject.log").read_text().splitlines() if l.startswith("field_rank "))
                    assert int(inj[3]) == i and int(inj[11]) == len(x)*1024
                    rec, hist, chain = old.inspect_fork(p,stem,cell,1,1)
                    stages.extend(rec)
                    history.extend(hist)
                    observed = (b/"stages.xyz").read_text().splitlines()
                    assert list(map(int,observed[0][1:].split())) == ids
                    root = next(l for l in observed if l.startswith("selected "))
                    np.testing.assert_allclose(poses[i],np.array(root.split()[3:],float).reshape(-1,3),rtol=0,atol=1e-9)
                    forks.append(dict(**k,anchor_serials=anchors,component=component,index=i,
                                      n_le2=sum(ds<=2),generation_best_A=ds.min(),first_le2_rank=min(ranks[ds<=2],default=""),
                                      selected_field_score=vals[i],rank_field_queries=int(inj[11]),
                                      second_A=next((v["rmsd_A"] for v in rec if v["stage"]=="second"),""),
                                      later_first_le2_step=next((v["index"] for v in rec if v["stage"]=="offered" and v["scope"]==1 and v["index"]>1 and v["rmsd_A"]<=2),""),
                                      **{kk:vv for kk,vv in chain.items() if kk!="stem"}))
                    for j in range(len(ds)):
                        row=dict(stem=stem,index=j,triangle_id=int(ts[j]),edge_residual_A2=residuals[j],
                                 rmsd_A=ds[j],field_score=vals[j],field_rank=ranks[j])
                        if writer is None:
                            writer = csv.DictWriter(stream,fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
    for arm in sorted({c["arm"] for c in cells}):
        for group in ("all","previous10","transfer_A1120","transfer_FABP4",*panel.CASES):
            rows = [r for r in cells if r["arm"]==arm and (group=="all" or group in (r["drug"],r["partition"]))]
            dist.append(distribution(rows,arm,group))
            if group=="all":
                dist.append(distribution(rows,arm,group,"molecule"))
    for file, rows in (("cells",cells),("modes",modes),("chemistry",chem),("rank_stages",forks),
                       ("stages",stages),("residency",history),("distributions",dist)):
        panel.table(out/(file+".csv"),rows)
    panel.old.write_json(out/"references.json",dict(
        archives={str(f.relative_to(panel.old.ROOT)):panel.digest(f) for f in
                  (panel.old.LIB["brd"],panel.BENCH/"cubic_grid/human_panel/depositions.tar.xz")},
        serials={name:dict(serials=cell["serials"],equivalent_frames=len(cell["refs"])) for name,(cell,*_) in geo.items()}))


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work","output"): parser.add_argument(key,type=Path)
    args=parser.parse_args();run(args.work.resolve(),args.output.resolve())
