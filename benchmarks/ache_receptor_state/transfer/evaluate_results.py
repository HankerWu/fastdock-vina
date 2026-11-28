"""Read the experimental geometry only after all transfer selections are frozen."""
import argparse
import csv
import json
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from predict import BLIND, RECEIVERS, SEEDS, model_positions, rms, table, save_json, sha, tick
from evaluate import distribution, references


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def geometry(sites, serial, ref):
    xyz=np.stack([sites[i] for i in serial])
    return xyz, min(rms(xyz, mapping) for mapping in ref)


def main(result, blind):
    start=tick()
    if sha((result/"raw.tar.xz").read_bytes()) != json.loads((result/"receipt.json").read_text())["raw_sha256"]:
        raise ValueError("Changed raw archive")
    poses = {}
    with tempfile.TemporaryDirectory() as tmp:
        refs, serial, fits=references(Path(tmp))
    for receiver in RECEIVERS:
        for seed in SEEDS:
            for arm in ("pro24k","pro48k"):
                pool=(blind/"raw"/receiver/f"{arm}_{seed}"/"dry.pdbqt").read_text().split("ENDMDL")[:-1]
                for m,block in enumerate(pool,1):
                    raw=model_positions(block.splitlines())
                    poses[receiver,seed,arm,m]=geometry(raw,serial,refs[receiver])
    job_data=read(result/"jobs.csv")
    coords, detail={},{}
    with tarfile.open(result/"raw.tar.xz") as ar:
        for j in job_data:
            rec,source,seed,m=j["receiver"],j["source"],int(j["seed"]),int(j["mode"])
            key=rec,source,seed,m
            parent=f"{rec}/{source}_{seed}_{m:02}"
            def geom(file, expected):
                data=ar.extractfile(parent+"/"+file).read()
                if sha(data)!=expected:
                    raise ValueError("Coordinate digest differs")
                return geometry(model_positions(data.decode().splitlines()),serial,refs[rec])
            _,din=geom("start.pdbqt",j["start_sha256"])
            _,dsource=poses[source,seed,"pro24k",m]
            xyz,dout=geom("end.pdbqt",j["end_sha256"]) if j["status"]=="ok" else (None,None)
            coords[key]=xyz
            detail[key]=dict(receiver=rec, source=source, source_seed=seed, source_mode=m, status=j["status"],
                             source_rmsd_A=dsource, start_rmsd_A=din, end_rmsd_A=dout,
                             local_change_A=dout-din if dout is not None else None)
    table(result/"geometry.csv",list(detail.values()))
    pool_rows,cells=[],[]
    rows=read(result/"selections.csv")
    for rec in RECEIVERS:
        for seed in SEEDS:
            for arm in ("recipient","cross","same","double"):
                subset=[r for r in rows if (r["receiver"],int(r["seed"]),r["arm"])==(rec,seed,arm)]
                near=[]
                for r in subset:
                    source,t,m=r["source"],int(r["source_seed"]),int(r["source_mode"])
                    if arm=="double":
                        _,d=poses[rec,seed,"pro48k",m]
                        near.append(dict(**r,source_rmsd_A=d,start_rmsd_A=d,end_rmsd_A=d,
                                         nearest_recipient_endpoint_A=None))
                        continue
                    key=rec,source,t,m
                    j=detail[key]
                    own=[xyz for (rr,ss,ts,_),xyz in coords.items() if rr==ss==rec and ts==seed and xyz is not None]
                    nearest=min(rms(coords[key],x) for x in own) if own and coords[key] is not None else None
                    near.append(dict(**r,**{k:j[k] for k in ("source_rmsd_A","start_rmsd_A","end_rmsd_A")},
                                     nearest_recipient_endpoint_A=nearest))
                pool_rows.extend(near)
                returned=sorted((r for r in near if r["decision"]=="retained"),key=lambda r:int(r["rank"]))
                by_rank=[r["end_rmsd_A"] for r in returned]
                def minimum(k):
                    return min((r[k] for r in near if r[k] is not None),default=None)
                ignored={s:sum(r["decision"]==s for r in near) for s in ("failed","window","duplicated","rank_cut")}
                donors=[r for r in near if (r["source"],int(r["source_seed"]))!=(rec,seed)]
                distinct=[r for r in donors if r["end_rmsd_A"] is not None and r["end_rmsd_A"]<2
                          and r["nearest_recipient_endpoint_A"]>=1]
                cells.append(dict(receiver=rec,seed=seed,arm=arm,planned=len(near),published=len(returned),**ignored,
                                  origin_best_A=minimum("source_rmsd_A"),start_best_A=minimum("start_rmsd_A"),
                                  local_best_A=minimum("end_rmsd_A"),top1_A=by_rank[0] if returned else None,
                                  best5_A=min(by_rank[:5]) if returned else None,oracle_A=min(by_rank) if returned else None,
                                  regret_A=by_rank[0]-min(by_rank) if returned else None,
                                  serious=by_rank[0]>5 if returned else None,distinct_donor_sub2=len(distinct),
                                  distinct_donor_sub2_returned=sum(r["decision"]=="retained" for r in distinct),
                                  distinct_donor_sub2_top1=sum(r["rank"]=="1" for r in distinct),
                                  selected_source=returned[0]["source"] if returned else None,
                                  selected_seed=returned[0]["source_seed"] if returned else None,
                                  selected_mode=returned[0]["source_mode"] if returned else None))
    table(result/"stages.csv",pool_rows)
    table(result/"cells.csv",cells)
    paired,groups=[],[]
    index={(r["receiver"],r["arm"],r["seed"]):r for r in cells}
    for rec in RECEIVERS:
        for arm in ("recipient","cross","same","double"):
            selected=[index[rec,arm,seed] for seed in SEEDS]
            for k in ("top1_A","best5_A","oracle_A","regret_A"):
                values=[r[k] for r in selected if r[k] is not None]
                groups.append(dict(receiver=rec,arm=arm,measure=k,empty=len(selected)-len(values),
                                   serious=sum(r["serious"] or False for r in selected),**distribution(values)))
        for seed in SEEDS:
            for arm in ("same","recipient","double"):
                c,b=index[rec,"cross",seed],index[rec,arm,seed]
                paired.append(dict(receiver=rec,seed=seed,arm1="cross",arm0=arm,
                                   **{k:c[k]-b[k] if c[k] is not None and b[k] is not None else None
                                      for k in ("top1_A","best5_A","oracle_A","regret_A")}))
    table(result/"groups.csv",groups)
    table(result/"paired.csv",paired)
    pair_summary=[]
    for rec in RECEIVERS:
        for arm in ("same","recipient","double"):
            for k in ("top1_A","best5_A","oracle_A","regret_A"):
                values=[r[k] for r in paired if r["receiver"]==rec and r["arm0"]==arm and r[k] is not None]
                pair_summary.append(dict(receiver=rec,arm1="cross",arm0=arm,measure=k,complete=len(values),
                                         **distribution(values)))
    table(result/"paired_summary.csv",pair_summary)
    save_json(result/"evaluation.json",dict(protein_fits=fits,total_wall_cpu_s=(tick()-start).tolist()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("result",type=Path)
    p.add_argument("--blind",type=Path,default=BLIND)
    a=p.parse_args()
    main(a.result.resolve(),a.blind)
