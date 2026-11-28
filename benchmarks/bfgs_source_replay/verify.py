"""Independent sparse-table audit from coordinates, event records and archived inputs."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np


def rows(directory, name):
    return list(csv.DictReader((directory/name).open()))


def reference(query, receiver):
    keys = sorted(set(query["ca"]) & set(receiver["ca"]))
    use = [k for k in keys if all(np.sum((np.array(p["ca"][k]) - p["center"])**2) < 14**2
                                  for p in (query, receiver))]
    x, y = [np.array([p["ca"][i] for i in use]) for p in (query, receiver)]
    assert len(x) >= 12
    a, _, b = np.linalg.svd((x-x.mean(0)).T @ (y-y.mean(0)))
    reflection = np.diag([1., 1., np.linalg.det(a @ b)])
    rotation = a @ reflection @ b
    return (np.array(query["mapped"])-x.mean(0)) @ rotation + y.mean(0)


def check(work, tables):
    historical = rows(tables, "locked_proposals.csv")
    added = rows(tables, "candidate_locked.csv")
    lineage = rows(tables, "published_lineage.csv")
    attrition = rows(tables, "path_attrition.csv")
    splits = rows(tables, "first_divergence.csv")
    max_delta = 0.
    n, n_published, n_splits, n_paths, n_same = 0, 0, 0, 0, 0
    for row in historical:
        case, driver = row["case"], row["driver"]
        a, b, _ = case.split("-")
        case_dir = work/case
        inputs = work/"fixtures/reciprocal"
        q,r = (json.loads((inputs/i/"preparation.json").read_text()) for i in (a,b))
        ref = reference(q,r)
        # Use query serials to pick the indices in actual prepared PDBQT record order.
        ligand = (inputs/a/"ligand.pdbqt").read_text()
        atoms = [int(l[6:11]) for l in ligand.splitlines() if l.startswith(("ATOM","HETATM"))
                 and int(l[6:11]) in q["serial"]]
        order = np.array([atoms.index(i) for i in q["serial"]])
        def error(xyz):
            diff = xyz[order] - ref
            return np.sqrt(np.einsum("sij,sij->s", diff, diff)/len(order)).min()
        def compare(value,target):
            nonlocal max_delta
            diff = abs(value - float(target))
            assert diff < 1.e-9, (value,target)
            max_delta = max(max_delta, diff)
        def decode(f):
            return error(np.fromstring(f[8], sep=" ").reshape(-1,3)*0.001)
        traces = {}
        for arm in ("e7031e9","2326430","ed0d42d","0bb827e"):
            words = [l.split("\t") for l in (case_dir/f"{driver}-{arm}.shadow").read_text().splitlines()]
            traces[arm] = words, np.array([decode(f) for f in words])
        profile = next(t for t in added if t["driver"]==driver and t["case"]==case)
        for arm, col, t in (("e7031e9","old",row),("2326430","repaired",row),
                            ("ed0d42d","stable",profile),("0bb827e","candidate",profile)):
            words, values = traces[arm]
            compare(values.min(),t[col+("_best" if t is row else "_best_rmsd")])
            assert (values <= 2).sum() == int(t[col+"_near"])
            assert sum(int(f[4]) for f in words) == int(t[col+"_evals"])
        words1, errors1 = traces["ed0d42d"]
        words2, errors2 = traces["0bb827e"]
        assert all(x[9:11]==y[9:11] for x,y in zip(words1,words2))
        energy = np.array([float(y[2])-float(x[2]) for x,y in zip(words1,words2)])
        compare(np.median(energy),profile["median_paired_energy"])
        compare(np.median(errors2-errors1),profile["median_paired_rmsd"])
        lines = (case_dir/f"{driver}.native").read_text().splitlines()
        native = [f.split("\t")[2:] for f in lines if f.startswith("H\t")]
        same = traces[driver][0]
        assert len(native)==len(same) and all(f==s[2:] for f,s in zip(native,same))
        n_same += len(same)
        n += len(energy)
        individual = []
        for line in lines:
            fields = line.split("\t")
            if fields[0]=="T":
                individual.append([])
            elif fields[0] in ("S","P","H","A","D","M"):
                individual[-1].append(fields)
        for i, events in enumerate(individual,1):
            t = next(z for z in attrition if z["trajectory"]==str(i) and z["case"]==case and z["driver"]==driver)
            for typ, key, geom in (("H","hunts","hunt"),("D","saved_attempts","saved")):
                hits = [decode(f) for f in events if f[0]==typ]
                assert len(hits) == int(t[key])
                assert sum(g<=2 for g in hits) == int(t[geom+"_near"])
                compare(min(hits),t[geom+"_best"])
            n_paths += 1
        prefix = f"reciprocal-{case}-{driver}.pdbqt"
        models = (case_dir/prefix).read_text().split("ENDMDL")[:-1]
        for rank, model in enumerate(models,1):
            t=next(x for x in lineage if x["case"]==case and x["driver"]==driver and int(x["rank"])==rank)
            coords={int(l[6:11]):[float(l[j:j+8]) for j in (30,38,46)] for l in model.splitlines()
                    if l.startswith(("ATOM","HETATM"))}
            compare(error(np.array([coords[i] for i in atoms])),t["published_rmsd"])
            n_published += 1
        if driver=="e7031e9":
            other=(case_dir/"2326430.native").read_text()
            for i,(aa,bb) in enumerate(zip(("\n".join(lines)+"\n").split("T\n")[1:],other.split("T\n")[1:]),1):
                def states(text):
                    for l in text.splitlines():
                        f=l.split("\t")
                        if f[0] in {"S","P","H","A","D","M"}:
                            yield (f[0],int(f[1]),f[2] if f[0]!="P" and len(f)>2 else "",tuple(f[3:4]))
                left,right=next((u,v) for u,v in zip(states(aa),states(bb)) if u!=v)
                t=next(x for x in splits if x["case"]==case and int(x["trajectory"])==i)
                assert left[0]==right[0]==t["state_stage"] and left[1]==right[1]==int(t["state_step"])
                n_splits += 1
    return dict(method="Independent SVD pocket fit and raw event parser; log/coordinate/energy checks",
                historical_driver_paths=n_paths, matched_native_controls=n_same,
                stable_candidate_pairs=n, first_divergences=n_splits, verified_published_ranks=n_published,
                max_abs_numeric_discrepancy=max_delta)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("restored",type=Path)
    p.add_argument("tables",type=Path)
    params=p.parse_args()
    print(json.dumps(check(params.restored,params.tables),indent=2))
