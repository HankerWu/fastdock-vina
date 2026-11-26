"""Replay the upstream selective-flexibility pathway at two crystal-diagnosed sites."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"confirmation"))
from run import measure, reference, write_csv

TARGETS = {"FABP4": (("2NNQ","3FR2"),60), "RBP4": (("1RBP","3FMZ"),98)}
SEEDS = (62483, 87911, 120997)


def ligand_only(text):
    lines=[]
    residue=False
    for line in text.splitlines():
        if line.startswith("BEGIN_RES"):
            residue=True
        elif line.startswith("END_RES"):
            residue=False
        elif not residue:
            lines.append(line)
    if residue: raise ValueError("Truncated flexible residue")
    return "\n".join(lines)+"\n"


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("fixtures","out","official"):
        p.add_argument(name,type=Path)
    p.add_argument("--jobs",type=int,default=2)
    a=p.parse_args()
    fixtures,out,official=[p.resolve() for p in (a.fixtures,a.out,a.official)]
    out.mkdir(parents=True,exist_ok=True)
    meta={p.parent.name:json.loads(p.read_text()) for p in fixtures.glob("*/preparation.json")}
    jobs=[]
    for target,(ids,site) in TARGETS.items():
        for pdb in ids:
            with (out/f"{pdb}.prep.log").open("w") as log:
                subprocess.run([sys.executable,"-m","meeko.cli.mk_prepare_receptor","--read_pdb",
                                str(fixtures/pdb/"chain-A.pdb"),"-f",f"A:{site}","-p",str(out/pdb)],
                               stdout=log,stderr=subprocess.STDOUT,check=True)
            def chemical_rows(text):
                return sorted(line[12:26]+line[30:54]+line[70:]
                              for line in text.splitlines() if line.startswith(("ATOM  ","HETATM")))
            rigid=(out/f"{pdb}_rigid.pdbqt").read_text()
            flex=(out/f"{pdb}_flex.pdbqt").read_text()
            if chemical_rows(rigid+"\n"+flex)!=chemical_rows((fixtures/pdb/"receptor.pdbqt").read_text()):
                raise ValueError(f"Prepared atom identities, types, charges or coordinates differ for {pdb}")
        jobs.extend(itertools.product(ids,ids,SEEDS,("rigid","flex")))

    def dock(job):
        q,r,seed,arm=job
        dest=out/f"{q}-{r}-{seed}-{arm}.pdbqt"
        cmd=[official,"--ligand",fixtures/q/"ligand.pdbqt","--seed",str(seed),"--cpu","1",
             "--max_evals","24000","--exhaustiveness","3","--num_modes","9","--energy_range","100","--out",dest]
        if arm=="rigid":
            cmd.extend(["--receptor",fixtures/r/"receptor.pdbqt"])
        else:
            cmd.extend(["--receptor",out/f"{r}_rigid.pdbqt","--flex",out/f"{r}_flex.pdbqt"])
        for axis,c in zip("xyz",meta[r]["center"]):
            cmd.extend(["--center_"+axis,str(c),"--size_"+axis,"24"])
        with dest.with_suffix(".log").open("w") as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        trimmed=dest.with_suffix(".ligand")
        trimmed.write_text(ligand_only(dest.read_text()))
        values=measure(trimmed,meta[q]["serial"],reference(meta[q],meta[r]))
        return [dict(query=q,receiver=r,target=meta[q]["target"],arm=arm,seed=seed,rank=i,rmsd=v)
                for i,v in enumerate(values,1)]
    with ThreadPoolExecutor(a.jobs) as pool:
        rows=[row for group in pool.map(dock,jobs) for row in group]
    write_csv(out/"modes.csv",rows)
    for target,arm in itertools.product(TARGETS,("rigid","flex")):
        selected=[r for r in rows if r["target"]==target and r["arm"]==arm]
        for cross in (True,False):
            hit=[len({(r["query"],r["receiver"],r["seed"]) for r in selected
                      if (r["query"]!=r["receiver"])==cross and r["rank"]<=n and r["rmsd"]<=2}) for n in (1,5,9)]
            print(target,arm,"cross" if cross else "cognate",*hit)


if __name__=="__main__":
    main()
