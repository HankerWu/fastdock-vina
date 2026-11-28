"""Write a PDBQT pool for each selection from the registered receiver-score rule."""
import argparse
import csv
import tarfile
from pathlib import Path

from predict import BLIND, RECEIVERS, SEEDS, save_json, sha, table, tick


def export(result, blind):
    start=tick()
    out=result/"published"
    out.mkdir(exist_ok=True)
    rows=list(csv.DictReader((result/"selections.csv").open()))
    costs=list(csv.DictReader((result/"costs.csv").open()))
    evidence, profiles=[],[]
    with tarfile.open(result/"raw.tar.xz") as ar:
        # Random access to LZMA here is avoidable; keep the 108 small ligand poses.
        ends={m.name:ar.extractfile(m).read().decode() for m in ar if m.name.endswith("/end.pdbqt")}
    for rec in RECEIVERS:
        for seed in SEEDS:
            for arm in ("cross","same","recipient","double"):
                t=tick()
                keep=sorted([r for r in rows if r["receiver"]==rec and int(r["seed"])==seed and r["arm"]==arm
                             and r["decision"]=="retained"],key=lambda r:int(r["rank"]))
                path=out/f"{rec}_{arm}_{seed}.pdbqt"
                if arm=="double":
                    path.write_bytes((blind/"raw"/rec/f"pro48k_{seed}"/"dry.pdbqt").read_bytes())
                else:
                    records=[]
                    for row in keep:
                        source=row["source"]
                        tag=f"{source}_{row['source_seed']}_{int(row['source_mode']):02d}"
                        coords=ends[f"{rec}/{tag}/end.pdbqt"]
                        records += [f"MODEL {row['rank']}", f"REMARK RECIPIENT SCORE: {row['score']} kcal/mol",
                                    f"REMARK SOURCE: {tag}", coords.rstrip(), "ENDMDL"]
                    path.write_text("\n".join(records)+("\n" if records else ""))
                wall,cpu=tick()-t
                evidence.append(dict(receiver=rec,seed=seed,arm=arm,modes=len(keep),pdbqt_sha256=sha(path.read_bytes()),
                                     export_wall_s=wall,export_cpu_s=cpu))
                charged=next(c for c in costs if c["receiver"]==rec and int(c["seed"])==seed and c["arm"]==arm)
                profiles.append(dict(receiver=rec,seed=seed,arm=arm,searches=charged["searches"],
                                     local_starts=charged["local_starts"],modes=len(keep),
                                     charged_wall_incl_export_s=float(charged["charged_wall_s"])+wall,
                                     charged_cpu_incl_export_s=float(charged["charged_cpu_s"])+cpu))
    table(result/"published.csv",evidence)
    table(result/"charged_cost.csv",profiles)
    save_json(result/"export_receipt.json",dict(total_wall_cpu_s=(tick()-start).tolist()))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("result",type=Path)
    p.add_argument("--blind",type=Path,default=BLIND)
    a=p.parse_args()
    export(a.result.resolve(),a.blind)
