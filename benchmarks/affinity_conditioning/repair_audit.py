"""Recompute the historical FKBP coordinate audits with the scoring operator of each binary."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

from replay import BENCH, COVER, RECEPTORS, audited_maps, calculate_components

HIST = BENCH / "receptor_linkage/aurora/internal_loops/rigid_fkbp"
TRIPLET = BENCH / "cavity_triplet/fkbp"
FIELD = BENCH / "cavity_field_rank/fkbp"
FRAGMENT = BENCH / "fragment_anchor/fkbp"
SOURCES = {
    HIST: [HIST/"pools.tar.xz"],
    TRIPLET: [HIST/"pools.tar.xz",TRIPLET/"blind_forks.tar.xz"],
    FIELD: [HIST/"pools.tar.xz",TRIPLET/"blind_forks.tar.xz",FIELD/"field_forks.tar.xz"],
    FRAGMENT: sorted(p for p in (FRAGMENT/"raw").glob("*.tar.xz") if p.name != "replay.tar.xz"),
    COVER: sorted(p for p in (COVER/"raw").glob("*.tar.xz") if p.name != "replay.tar.xz"),
}
PATTERN = re.compile(r"pools/([^/]+)/(1FK[GH])_(\d+).pdbqt$")


def load_pools(archive):
    raw = {}
    with tarfile.open(archive, "r|xz") as tar:
        for member in tar:
            match = PATTERN.fullmatch(member.name)
            if match:
                arm, rec, compound = match.groups()
                raw[(arm,rec,int(compound))] = tar.extractfile(member).read()
    return raw


def run(out):
    timer = time.perf_counter(),time.process_time()
    prep = json.loads((COVER/"results/inputs.json").read_text())
    archives = {a:load_pools(a) for group in SOURCES.values() for a in group}
    audit, summary, directions, cache = [],[],[],{}
    with tempfile.TemporaryDirectory() as tmp:
        maps = {r:audited_maps(r,prep,Path(tmp)) for r in RECEPTORS}
        for source, files in SOURCES.items():
            pools = {k:v for f in files for k,v in archives[f].items()}
            prior = pd.read_csv(source/"results/modes.csv")
            label = source.relative_to(BENCH).as_posix()
            cellrows=[]
            for (arm,rec,compound), group in prior.groupby(["arm","receiver","compound"],sort=False):
                raw = pools[(arm,rec,compound)]
                digest = hashlib.sha256(raw).hexdigest()
                k = (rec,digest)
                if k not in cache:
                    result=[]
                    for piece in raw.decode("ascii").split("ENDMDL"):
                        if not piece.strip():
                            continue
                        lines=piece.strip().splitlines()
                        reported=float(next(x for x in lines if x.startswith("REMARK VINA RESULT")).split()[3])
                        values = [sum(calculate_components(lines,maps[rec],operator=op).values())
                                  for op in ("separate","combined")]
                        result.append((reported,*values))
                    cache[k] = result
                before = group["ad4_from_xyz" if "ad4_from_xyz" in group else "from_xyz"]
                printed = group["reported" if "reported" in group else "score"]
                operator = "separate" if arm.startswith("official") else "combined"
                native_scores=[]
                for rank, ((reported,s,c),archived,native) in enumerate(zip(cache[k],before,printed,strict=True),1):
                    corrected = s if operator=="separate" else c
                    if abs(archived-s)>1e-8 or abs(native-reported)>1e-8:
                        raise ValueError(f"Historical coordinate/pool mismatch: {label} {k} {rank}")
                    native_scores.append(corrected)
                    cellrows.append(dict(panel=label,arm=arm,receiver=rec,compound=compound,rank=rank,
                                         operator=operator,printed=reported,old=archived,corrected=corrected,
                                         delta=corrected-archived))
                if len(native_scores)!=9:
                    raise ValueError("Unexpected number of archived modes")
            table=pd.DataFrame(cellrows)
            audit.extend(cellrows)
            winners=table.groupby(["arm","receiver","compound"])
            flips=sum(g.old.idxmin()!=g.corrected.idxmin() for _,g in winners)
            order=sum(list(g.sort_values(["old","rank"])["rank"])!=list(g.sort_values(["corrected","rank"])["rank"])
                      for _,g in winners)
            first=table[table["rank"]==1].set_index(["arm","receiver","compound"])
            pairs=pd.read_csv(source/"results/pairs.csv")
            for r in pairs.itertuples():
                def margin(col):
                    recs=RECEPTORS if r.receiver=="mean" else (r.receiver,)
                    return np.mean([first.loc[(r.arm,rec,r.weaker),col]-
                                    first.loc[(r.arm,rec,r.stronger),col] for rec in recs])
                old, new, published=map(margin,("old","corrected","printed"))
                if not np.isclose(r.xyz_margin, old, atol=1e-8):
                    raise ValueError("Pre-existing pair margin mismatch")
                directions.append(dict(panel=label,arm=r.arm,receiver=r.receiver,stronger=r.stronger,
                                       weaker=r.weaker,old=old,corrected=new,reported=published))
            subset=[r for r in directions if r["panel"]==label]
            summary.append(dict(panel=label,observations=len(table),pools=len(winners),pairs=len(subset),
                                corrected_operator_modes=int((table.operator=="combined").sum()),
                                max_score_change=float(table.delta.abs().max()),
                                max_old_print_error=float((table.old-table.printed).abs().max()),
                                max_new_print_error=float((table.corrected-table.printed).abs().max()),
                                recomputed_top_index_changes=flips,recomputed_complete_order_changes=order,
                                pair_direction_changes=sum(np.sign(r["old"])!=np.sign(r["corrected"]) for r in subset)))
    out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(summary).to_csv(out/"validation_repair.csv",index=False)
    for label,records in (("modes",audit),("pairs",directions)):
        payload=pd.DataFrame(records).to_csv(index=False).encode()
        (out/f"validation_{label}.csv.gz").write_bytes(gzip.compress(payload,mtime=0))
    costs = dict(unique_pool_receptor_inputs=len(cache),wall_s=time.perf_counter()-timer[0],
                 cpu_s=time.process_time()-timer[1])
    (out/"validation_cost.json").write_text(json.dumps(costs,indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    run(parser.parse_args().out)
