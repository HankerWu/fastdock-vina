"""Second-stage fixed-pool and fixed-assay comparison with the intrinsic torsion prior."""
import argparse
from functools import lru_cache
from pathlib import Path
import tempfile

from affinity import ARMS, HERE, assay_inputs, save, summarize
from bonded import add_bonded_prior
from scorer import Receptor, score_pose
from evaluate import archive_cases, blocks, read


def replay(output):
    first_pass = {}
    for r in read(HERE / "results/affinity_poses.csv"):
        first_pass[(r["target"], r["ligand"], r["receptor"], int(r["seed"]))] = r

    @lru_cache(maxsize=None)
    def bonded(protein):
        return add_bonded_prior(Receptor(protein))

    results = []
    for target, key, measured, rec, seed, prot, top, data in assay_inputs():
        r = first_pass[target,key,rec,seed]
        b = bonded(prot.decode())
        terms = score_pose(b, top.decode(), blocks(data.decode())[0])
        scale, native = float(r["scale"]), float(r["native"])
        result = {**r, "measurement": measured, "rotors": int(r["rotors"]), **{k:float(r[k]) for k in ARMS}}
        result["ensemble"] = native + scale*terms["raw_ensemble"]
        result["rotor_only"] = native + scale*terms["raw_rotor_only"]
        results.append(result)
    sums, _ = summarize(results)
    save(output/"bonded_affinity_summary.csv", [r for r in sums if r["arm"] in ("ensemble","rotor_only")])
    save(output/"bonded_affinity_poses.csv", [dict(target=r["target"],ligand=r["ligand"],receptor=r["receptor"],
        seed=r["seed"],ensemble=r["ensemble"],rotor_only=r["rotor_only"]) for r in results])

    old = {(r["panel"],r["design"],r["query"],r["receiver"],r["seed"],int(r["rank"])):r
           for r in read(HERE.parent/"hydroxyl_partition/results/modes.csv")}
    poses, cells = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for case, sites, pool, table in archive_cases(Path(tmp)):
            receptor = bonded((sites/case["receiver"]/"receptor.pdbqt").read_text())
            ligand = (sites/case["query"]/"ligand.pdbqt").read_text()
            values = []
            for i, (part, ref) in enumerate(zip(blocks(pool),table,strict=True),1):
                id_ = tuple(case[k] for k in ("panel","design","query","receiver","seed"))+(i,)
                original = old.pop(id_)
                terms = score_pose(receptor,ligand,part)
                score = float(original["native"]) + float(original["scale"])*terms["raw_ensemble"]
                row = dict(**case, rank=i,rmsd_A=float(ref["rmsd_A"]), native=float(original["native"]),
                           ensemble=float(original["ensemble"]),bonded=score)
                values.append(row)
            cell = dict(**case,oracle=min(p["rmsd_A"] for p in values))
            for field in ("native","ensemble","bonded"):
                order = sorted(values,key=lambda x:(x[field],x["rank"]))
                cell[field+"_rank"] = order[0]["rank"]
                cell[field+"_top"] = order[0]["rmsd_A"]
                cell[field+"_top5"] = min(p["rmsd_A"] for p in order[:5])
            cells.append(cell)
            poses.extend(values)
    if old:
        raise ValueError("Not all archived poses accounted for")
    save(output/"bonded_pose_cells.csv",cells)
    save(output/"bonded_pose_modes.csv",poses)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    replay(args.output)
