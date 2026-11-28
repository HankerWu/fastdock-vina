"""Matched native Vina scores and local relaxations on three measured AChE structures."""
import argparse
import hashlib
import json
import lzma
import re
import resource
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
OLD = BENCH / "water_rank_envelope/transfer"
SOURCE = BENCH / "receptor_amide_states/transfer/transfer_evidence.tar.xz"
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as chemistry

CODES = {"4EY5": "HUP", "4EY6": "GNT", "4EY7": "E20"}
ARMS = {"pro24k": OLD, "pro48k": OLD / "pro48k", "official24k": OLD / "official24k"}
LABELS = ["Estimated Free Energy of Binding", "Final Intermolecular Energy",
          "Final Total Internal Energy", "Torsional Free Energy", "Unbound System's Energy"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tick():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return np.array([time.perf_counter(), time.process_time() + child.ru_utime + child.ru_stime])


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def xyz(row):
    return np.array([float(row[i:i+8]) for i in (30, 38, 46)])


def model_positions(lines):
    return {int(l[6:11]): xyz(l) for l in lines if l.startswith(("ATOM  ", "HETATM"))}


def move(ligand, coords, rotation, translation):
    moved = []
    for line in ligand.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            pos = coords[int(line[6:11])] @ rotation + translation
            line = line[:30] + "".join(f"{v:8.3f}" for v in pos) + line[54:]
        moved.append(line)
    return "\n".join(moved) + "\n"


def fit(x, y):
    u, _, vt = np.linalg.svd((x - x.mean(0)).T @ (y - y.mean(0)))
    u[:, -1] *= np.sign(np.linalg.det(u @ vt))
    rot = u @ vt
    return rot, y.mean(0) - x.mean(0) @ rot


def pocket(rows, ligand):
    return dict(ca={l[22:27]+l[17:20]: xyz(l) for l in rows
                    if l.startswith("ATOM  ") and l[12:16] == " CA " and l[21] == "A"},
                center=np.mean(list(ligand.values()), axis=0))


def transform(pockets, a, b):
    if a == b:
        return np.eye(3), np.zeros(3), 0., len(pockets[a]["ca"])
    return chemistry.alignment(pockets[a], pockets[b])


def prepare(work):
    clock = tick()
    raw, hashes = work / "raw", {}
    raw.mkdir(parents=True, exist_ok=False)
    for pdb in CODES:
        (raw / pdb).mkdir()
    archive_sha = digest(SOURCE)
    old_meta = json.loads((OLD/"manifest.json").read_text())
    assert archive_sha == old_meta["source_archive_sha256"]
    (work/"4EY5.pdb").write_bytes(lzma.decompress((OLD/"4EY5.pdb.xz").read_bytes()))
    with tarfile.open(SOURCE) as archive:
        for pdb in ["4EY6", "4EY7"]:
            for dest, member in [(work/f"{pdb}.pdb", f"novel/raw/{pdb}.pdb"),
                                 (raw/pdb/"receptor.pdbqt", f"novel/{pdb}/receptor.pdbqt"),
                                 (raw/pdb/"preparation.json", f"novel/{pdb}/preparation.json")]:
                dest.write_bytes(archive.extractfile(member).read())
            previous = json.loads((raw/pdb/"preparation.json").read_text())
            assert digest(work/f"{pdb}.pdb") == previous["pdb_sha256"]
            assert previous["receptor_status"] == "ready"
    for name, path in [("4EY5.pdb",work/"4EY5.pdb"), ("HUP.cif",OLD/"HUP.cif"),
                       ("HUP.pdbqt",OLD/"HUP.pdbqt"), ("receptor.pdbqt",raw/"4EY6/receptor.pdbqt"),
                       ("preparation.json",raw/"4EY6/preparation.json")]:
        assert digest(path) == old_meta["hashes"][name]
    rows = {p: (work/f"{p}.pdb").read_text().splitlines() for p in CODES}
    ligands = {p: chemistry.dep_records(rows[p], code)[0] for p, code in CODES.items()}
    pockets = {p: pocket(rows[p], ligands[p]) for p in CODES}
    graph = chemistry.ligand(OLD/"HUP.cif", [ligands["4EY5"]], work/"HUP.check.pdbqt")
    assert digest(work/"HUP.check.pdbqt") == digest(OLD/"HUP.pdbqt") and len(graph["mapped"]) == 1
    truth = np.array(graph["mapped"][0])
    serials = graph["serial"]
    template = (OLD/"HUP.pdbqt").read_text()
    original = model_positions(template.splitlines())
    rot, tr = fit(np.array([original[i] for i in serials]), truth)
    projected = {i: v @ rot + tr for i,v in original.items()}
    included, guards = [], {}
    closest = lambda l: float(np.min(np.linalg.norm(truth - xyz(l), axis=1)))
    omitted = [l for l in rows["4EY5"] if l.startswith(("ATOM  ", "HETATM")) and
               l[76:78].strip() not in ("H", "D") and l[17:20] not in ("HOH", "DOD") and
               not (l[21]=="A" and (l.startswith("ATOM  ") or l[17:20] == "HUP"))]
    nearest = min(omitted, key=closest)
    prep = chemistry.receptor(rows["4EY5"], raw/"4EY5", "HUP", np.array(list(ligands["4EY5"].values())))
    prep["nearest_omitted_nonwater_A"] = [closest(nearest), nearest]
    guards["4EY5"] = prep
    save(raw/"4EY5/preparation.json", prep)
    starts, ids = {}, {}
    for arm,path in ARMS.items():
        assert digest(path/"dry.pdbqt") == json.loads((path/"manifest.json").read_text())["hashes"]["dry.pdbqt"]
        parts = (path/"dry.pdbqt").read_text().split("ENDMDL")[:-1]
        old_modes = pd.read_csv(path/"modes.csv").query("order == 11")
        ids[arm] = int(old_modes.loc[old_modes.rmsd.idxmin(), "mode"])
        for i, block in enumerate(parts,1):
            starts[f"{arm}_{i:02d}"] = model_positions(block.splitlines())
    assert len(starts) == 26
    targets = {}
    for pdb in CODES:
        if pdb == "4EY5" and (closest(nearest) < 8 or (prep["closest_repaired_ligand_distance"] is not None
                                                       and prep["closest_repaired_ligand_distance"] < 8)):
            guards["4EY5"]["excluded"] = "local omitted or repaired nonwater heavy atom"
            continue
        included.append(pdb)
        r,t,err,n = transform(pockets,"4EY6",pdb)
        q,v,qerr,qn = transform(pockets,"4EY5",pdb)
        target = truth @ q + v
        rc,tc,_,_ = transform(pockets,"4EY5","4EY6")
        via = (truth @ rc + tc) @ r + t
        targets[pdb] = dict(truth=target.tolist(), center=pockets[pdb]["center"].tolist(),
                            from_4ey6_rotation=r.tolist(), from_4ey6_translation=t.tolist(),
                            from_4ey5_rotation=q.tolist(), from_4ey5_translation=v.tolist(),
                            fit_from_4ey6=[err,n],fit_from_4ey5=[qerr,qn],
                            fit_cycle_rmsd_A=float(np.sqrt(np.mean(np.sum((via-target)**2,axis=1)))))
        for key, coords in dict(starts, projection=projected).items():
            job = raw / pdb / key
            job.mkdir()
            matrix, translation = (q,v) if key=="projection" else (r,t)
            (job/"input.pdbqt").write_text(move(template, coords, matrix, translation))
    for pdb in CODES:
        for p in (work/f"{pdb}.pdb", raw/pdb/"receptor.pdbqt",raw/pdb/"preparation.json"):
            hashes[str(p.relative_to(work))] = digest(p)
    info = dict(sources_archive_sha256=archive_sha,prior_pools={a:digest(p/"dry.pdbqt") for a,p in ARMS.items()},
                hup_cif_sha256=digest(OLD/"HUP.cif"),hup_pdbqt_sha256=digest(OLD/"HUP.pdbqt"),
                hashes=hashes,included=included,cognate_guard=guards,targets=targets,serials=serials,
                oracle_mode=ids, preparation_wall_cpu_s=(tick()-clock).tolist())
    save(work/"inputs.json",info)


def call(exe, args, job, name, costs):
    start = tick()
    with (job/f"{name}.log").open("w") as f:
        p = subprocess.run([str(exe), *args], cwd=job, stdout=f, stderr=subprocess.STDOUT)
    wall, cpu = tick()-start
    costs.append(dict(receiver=job.parent.name,start=job.name,task=name,wall_s=wall,cpu_s=cpu,
                      exit_code=p.returncode))
    if p.returncode:
        raise ValueError(f"{name} failed with code {p.returncode}")


def score(log):
    text = log.read_text()
    return [float(re.search(re.escape(s)+r"\s*:\s*(-?[\d.]+)", text)[1]) for s in LABELS]


def geom(path, target, center, serials):
    positions = model_positions(path.read_text().splitlines())
    xyz_all = np.stack(list(positions.values()))
    if np.abs(xyz_all-center).max() >= 18:
        raise ValueError(f"Pose outside diagnostic box: {path}")
    x = np.array([positions[i] for i in serials])
    return float(np.sqrt(np.mean(np.sum((x-target)**2,axis=1))))


def run(work, binary, observer):
    meta = json.loads((work/"inputs.json").read_text())
    records, costs, failed = [], [], []
    for pdb in meta["included"]:
        m = meta["targets"][pdb]
        ref = np.array(m["truth"])
        frame = ["--receptor","../receptor.pdbqt","--cpu","1","--seed","1103501"]
        for axis, x in zip("xyz",m["center"]):
            frame += ["--center_"+axis,str(x),"--size_"+axis,"36"]
        for job in sorted((work/"raw"/pdb).iterdir()):
            if not job.is_dir():
                continue
            try:
                before = geom(job/"input.pdbqt",ref,m["center"],meta["serials"])
                args = frame + ["--ligand","input.pdbqt"]
                call(binary,args+["--score_only","--unbound_energy","0"],job,"initial",costs)
                call(binary,args+["--local_only","--out","output.pdbqt"],job,"local",costs)
                call(binary,frame+["--ligand","output.pdbqt","--score_only","--unbound_energy","0"],job,"final",costs)
                call(observer,args+["--local_only","--out","observed.pdbqt"],job,"observer",costs)
                if digest(job/"output.pdbqt") != digest(job/"observed.pdbqt"):
                    raise ValueError("observer changes production coordinates")
                observations = re.findall(r"REFINE_OBSERVE\s+(\S+)\s+(\S+)\s+(\d+)", (job/"observer.log").read_text())
                if not observations:
                    raise ValueError("no observed gradients")
                after = geom(job/"output.pdbqt",ref,m["center"],meta["serials"])
                row=dict(receiver=pdb,start=job.name,rmsd_in=before,rmsd_out=after,
                         objective_in=float(observations[0][0]),objective_out=float(observations[-1][1]),
                         gradient_evals=sum(int(a[2]) for a in observations),
                         input_sha256=digest(job/"input.pdbqt"),output_sha256=digest(job/"output.pdbqt"))
                for name,stage in (("initial","in"),("final","out")):
                    for key,val in zip(("score","inter","intra","torsion","unbound"),score(job/f"{name}.log")):
                        row[key+"_"+stage]=val
                records.append(row)
            except (ValueError,FileNotFoundError,IndexError,TypeError) as e:
                failed.append(dict(receiver=pdb,start=job.name,error=str(e)))
    df=pd.DataFrame(records)
    df.to_csv(work/"endpoints.csv",index=False)
    pd.DataFrame(costs).to_csv(work/"costs.csv",index=False)
    save(work/"run.json",dict(binary_sha256=digest(binary),observer_sha256=digest(observer),failures=failed,
                             planned=27*len(meta["included"]),complete=len(df)))
    results=[]
    for pdb in meta["included"]:
        group=df[df.receiver==pdb].set_index("start")
        if "projection" not in group.index:
            continue
        proj=group.loc["projection"]
        for arm,oracle in meta["oracle_mode"].items():
            pool=group.loc[group.index.str.startswith(arm+"_")]
            if f"{arm}_{oracle:02d}" not in pool.index or f"{arm}_01" not in pool.index:
                continue
            a,b=pool.loc[f"{arm}_01"],pool.loc[f"{arm}_{oracle:02d}"]
            for stage in ("in","out"):
                score_col, rmsd_col="score_"+stage,"rmsd_"+stage
                winner=pool[score_col].idxmin()
                best=pool[rmsd_col].min()
                results.append(dict(receiver=pdb,arm=arm,stage=stage,n=len(pool),oracle_id=oracle,
                                   delta_oracle_original=b[score_col]-a[score_col],winner=winner,
                                   winner_rmsd=pool.loc[winner,rmsd_col],pool_oracle=best,
                                   regret=pool.loc[winner,rmsd_col]-best,projected_rmsd=proj[rmsd_col],
                                   projected_score=proj[score_col],best_pool_score=pool[score_col].min(),
                                   projection_minus_best=proj[score_col]-pool[score_col].min()))
    pd.DataFrame(results).to_csv(work/"contrasts.csv",index=False)
    same_state_control(work, df)


def same_state_control(work, df):
    rows = []
    state = df[df.receiver == "4EY6"].set_index("start")
    for arm, source in ARMS.items():
        blocks = (source/"dry.pdbqt").read_text().split("ENDMDL")[:-1]
        for rank, block in enumerate(blocks, 1):
            key = f"{arm}_{rank:02d}"
            if key not in state.index:
                continue
            original_inter = float(re.search(r"REMARK INTER:\s*([-0-9.]+)", block)[1])
            original_score = float(re.search(r"REMARK VINA RESULT:\s*([-0-9.]+)", block)[1])
            rows.append(dict(arm=arm, original_rank=rank, printed_native=original_score,
                             printed_inter=original_inter, replay_inter=state.loc[key, "inter_in"],
                             replay_score=state.loc[key,"score_in"]))
    table = pd.DataFrame(rows)
    table["replay_rank"] = table.groupby("arm").replay_score.rank(method="first").astype(int)
    table["replay_minus_printed_inter"] = table.replay_inter - table.printed_inter
    table.to_csv(work/"same_state.csv",index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=["prepare","run"])
    parser.add_argument("work",type=Path)
    parser.add_argument("--pro",type=Path)
    parser.add_argument("--observer",type=Path)
    a=parser.parse_args()
    if a.step=="prepare":
        prepare(a.work.resolve())
    else:
        if not (a.pro and a.observer):
            parser.error("run needs --pro and --observer")
        run(a.work.resolve(),a.pro.resolve(),a.observer.resolve())
