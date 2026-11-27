"""Time the uninstrumented production defaults; compare their released pools to the count builds."""
import argparse
from pathlib import Path
from assay_runner import RECEPTORS, SEED, sha, table, time_one


def replay(work,exe):
    dest=work/"ordinary_pro"
    dest.mkdir(exist_ok=True)
    rows=[]
    for rec in RECEPTORS:
        for n in range(4,12):
            for arm in ("base","high"):
                stem=dest/f"{arm}_{rec}_{n}"
                out=stem.with_suffix(".pdbqt")
                cmd=[str(exe),"--ligand",f"affinity/{n}.pdbqt","--maps",f"maps/{rec}/receptor",
                     "--scoring","ad4","--cpu","1","--exhaustiveness","3","--max_evals",str(48000 if arm=="high" else 24000),
                     "--num_modes","9","--energy_range","100","--seed",str(SEED),"--out",str(out.relative_to(work))]
                wall,cpu=time_one(cmd,work,stem.with_suffix(".log"))
                rows.append(dict(arm=arm,receiver=rec,compound=n,wall_s=wall,cpu_s=cpu,sha256=sha(out),executable_sha256=sha(exe)))
    table(work/"unmodified_pro.csv",rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("exe",type=Path)
    a=p.parse_args()
    replay(a.work.resolve(),a.exe.resolve())
