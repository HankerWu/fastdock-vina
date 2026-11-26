"""Recompute the chain-copy diagnostic directly from the deposited references and raw pools."""
import argparse,csv,hashlib,json,math,statistics,sys
from collections import defaultdict
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from audit import coordinates, protein_alignment


def read_csv(p):
    with p.open() as stream:
        return list(csv.DictReader(stream))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("original",type=Path)
    p.add_argument("copies",type=Path)
    a=p.parse_args()
    mode_rows=read_csv(a.copies/"modes.csv")
    groups=defaultdict(list)
    for r in mode_rows:
        groups[r["chain"],r["query"],r["receiver"],r["seed"]].append(r)
    ct={(r["chain"],r["receiver"],r["seed"]):r for r in read_csv(a.copies/"contacts.csv")}
    energies=defaultdict(list)
    error=0.; total=0; checked=0
    for chain in ("A","B"):
        root=a.original/"inputs" if chain=="A" else a.copies/"B"
        docking=a.original/"dry" if chain=="A" else a.copies/"B"
        dm={}
        for code in ("7MRH","7MRG"):
            d=root/code; metadata=d/"metadata.json"
            dm[code]=json.loads((metadata if metadata.exists() else d/"query_only.json").read_text())
            lines=(d/(code+".pdb")).read_text().splitlines()
            dm[code]["ca"]={(l[22:27],l[17:20]):coordinates(l) for l in lines if l.startswith("ATOM  ")
                and l[21]=="A" and l[12:16].strip()=="CA" and l[16] in (" ","A")}
        for row in read_csv(a.copies/"cells.csv"):
            if row["chain"]!=chain:continue
            q,r,s=row["query"],row["receiver"],row["seed"]
            blob=(docking/"poses"/f"{q}-{r}-{s}-dry.pdbqt").read_bytes()
            assert hashlib.sha256(blob).hexdigest()==row["pdbqt_sha256"]
            maps=protein_alignment(q,r,dm)
            nz=[coordinates(l) for l in (docking/r/"receptor.pdbqt").read_text().splitlines() if
                l.startswith("ATOM") and l[17:20]=="LYS" and int(l[22:26])==60 and l[12:16].strip()=="NZ"][0]
            metrics=[]
            for model in blob.decode().split("ENDMDL"):
                lines=model.splitlines()
                pos={int(l[6:11]):coordinates(l) for l in lines if l.startswith(("ATOM  ","HETATM"))}
                if not pos:continue
                rms=min(math.sqrt(math.fsum(math.dist(pos[i],xyz)**2 for i,xyz in zip(dm[q]["serial"],m))/len(m)) for m in maps)
                energy=float(next(l for l in lines if l.startswith("REMARK VINA RESULT:")).split()[3])
                reference=groups[chain,q,r,s][len(metrics)]
                assert float(reference["score"])==energy
                error=max(error,abs(rms-float(reference["rmsd_A"])))
                metrics.append((rms,min(math.dist(nz,pos[i]) for i in (4,5)) if q=="7MRH" else None,energy))
            assert len(metrics)==int(row["modes"])
            total+=len(metrics); checked+=1
            energies[chain,r,q].append(metrics[0][2])
            if q=="7MRH":
                report=ct[chain,r,s]
                for label,rank in (("first",0),("nearest",min(range(len(metrics)),key=lambda k:metrics[k][0]))):
                    assert rank+1==int(report[label+"_rank"])
                    error=max(error,abs(float(report[label+"_rmsd_A"])-metrics[rank][0]),
                              abs(float(report[label+"_lys60_A"])-metrics[rank][1]))
    for r in read_csv(a.copies/"orders.csv"):
        key=r["chain"],r["receiver"]
        strong,weak=(statistics.median(energies[*key,q]) for q in ("7MRH","7MRG"))
        assert strong==float(r["ZMJ_77nM"]) and weak==float(r["ZMM_875nM"])
        assert str(strong<weak)==r["correct"]
    assert error<1e-8
    print(json.dumps(dict(pools=checked,modes=total,max_geometry_discrepancy_A=error),indent=2))

if __name__=="__main__":
    main()
