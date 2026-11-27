"""Run matched native AD4 local search on two published starts and the feasible 5 intermediate."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import resource
import subprocess
import tarfile
import tempfile
import time

import numpy as np

from manifold import csv_out
from replay import COVER, HUMAN, LEGACY, calculate_components, files, models, pool_sets, read_map
from verify import Ligand

HERE = Path(__file__).resolve().parent
KEY = ("panel", "arm", "seed", "receiver")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def prep(work):
    meta = json.loads((COVER/"results/inputs.json").read_text())
    lig_names = [f"affinity/{n}.pdbqt" for n in (5,6)]
    data = files(HUMAN/"affinity-pools.tar.xz",
                 meta["archives"]["benchmarks/cubic_grid/human_panel/affinity-pools.tar.xz"],
                 {name: meta["files"][name] for name in lig_names})
    ligands, templates, maps = {}, {}, {}
    (work/"evidence").mkdir()
    for n in (5,6):
        template = data[f"affinity/{n}.pdbqt"].decode()
        lig = Ligand(template)
        edges = "".join(f"{lig.order[i]} {lig.order[j]}\n" for i,j in lig.bonds)
        assert sha(edges.encode()) == meta["files"][f"affinity/{n}.bonds"]
        (work/"evidence"/f"{n}.bonds").write_text(edges)
        ligands[n], templates[n] = lig, template
    for rec in ("1FKG","1FKH"):
        names = {name: val for name,val in meta["files"].items() if name.startswith(f"maps/{rec}/")}
        archive = LEGACY/f"maps-{rec}.tar.xz"
        payload = files(archive, meta["archives"][str(archive.relative_to(HUMAN.parents[2]))], names)
        maps[rec]={}
        for name, content in payload.items():
            dest = work/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(content)
            maps[rec][dest.name.split(".")[1]] = read_map(dest)
    return ligands, templates, maps


def coordinates(template, ligand, xyz):
    sites = dict(zip(ligand.order, xyz, strict=True))
    lines = []
    for line in template.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            s = int(line[6:11])
            line = line[:30]+"".join(f"{v:8.3f}" for v in sites[s])+line[54:]
        lines.append(line)
    return "\n".join(lines)+"\n"


def invoke(exe, args, work, log):
    t0 = time.perf_counter()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    run = subprocess.run([str(exe), *args], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    wall = time.perf_counter()-t0
    log.write_bytes(run.stdout)
    return run.returncode, wall, cpu.ru_utime+cpu.ru_stime-usage.ru_utime-usage.ru_stime


def energy_log(text):
    names = ("Estimated Free Energy of Binding", "(1) Final Intermolecular Energy",
             "(2) Final Total Internal Energy", "(3) Torsional Free Energy", "(4) Unbound System")
    columns = []
    for label in names:
        values = [float(line.split(":",1)[1].split()[0]) for line in text.splitlines() if line.startswith(label)]
        if len(values) != 2:
            raise ValueError(f"Expected before/after for {label}")
        columns.append(values)
    return [dict(zip(("S", "inter", "intra", "torsion", "unbound"), e, strict=True))
            for e in zip(*columns, strict=True)]


def run_leg(work, name, ligand, number, text, receptor, maps, binaries, capture_rejection):
    native, observer, graph = binaries
    folder = work/"evidence"/name
    folder.mkdir()
    inp, out = folder/"input.pdbqt", folder/"output.pdbqt"
    inp.write_text(text)
    bonds = work/"evidence"/f"{number}.bonds"
    normal = subprocess.check_output([graph, inp])
    explicit = subprocess.check_output([graph, inp, bonds])
    (folder/"inferred.graph").write_bytes(normal)
    (folder/"explicit.graph").write_bytes(explicit)
    if normal != explicit:
        return dict(job=name, receiver=receptor, ligand=number, status="graph_mismatch",
                    input_sha256=sha(inp.read_bytes()), inferred_sha256=sha(normal), explicit_sha256=sha(explicit))
    opts = ["--ligand",str(inp.relative_to(work)),"--maps",f"maps/{receptor}/receptor",
            "--scoring","ad4","--local_only","--cpu","1","--seed","42","--verbosity","2",
            "--out",str(out.relative_to(work))]
    (folder/"command.json").write_text(json.dumps(opts)+"\n")
    if capture_rejection:
        status, _, _ = invoke(native, [*opts, "--ligand_bonds", str(bonds.relative_to(work))],
                              work, work/"evidence/explicit_cli_rejection.log")
        if status == 0 or "ligand_bonds requires" not in (work/"evidence/explicit_cli_rejection.log").read_text():
            raise ValueError("Unexpected change in rigid AD4 explicit-graph CLI support")
    status, wall, cpu = invoke(native, opts, work, folder/"native.log")
    if status:
        raise ValueError(f"Production search failed with {status}: {(folder/'native.log').read_text()}")
    obs = folder/"observer.pdbqt"
    obs_opts = opts[:-1]+[str(obs.relative_to(work))]
    status, obs_wall, obs_cpu = invoke(observer, obs_opts, work, folder/"observer.log")
    if status or obs.read_bytes() != out.read_bytes():
        raise ValueError(f"Inconsistent observed and unmodified AD4 coordinates in {name}")
    obs.unlink()  # identical to saved output
    log = (folder/"native.log").read_text()
    endpoints = energy_log(log)
    obs_log=(folder/"observer.log").read_text()
    if energy_log(obs_log) != endpoints:
        raise ValueError("Observer changed a reported energy")
    observed = [line.split("REFINE_OBSERVE ", 1)[1].split() for line in obs_log.splitlines()
                if "REFINE_OBSERVE " in line]
    if not observed:
        raise ValueError(f"Missing true BFGS objective: {obs_log}")
    e0, e1 = float(observed[0][0]), float(observed[-1][1])
    calls = sum(int(row[2]) for row in observed)
    result = dict(job=name, receiver=receptor, ligand=number, status="ok", atoms=len(ligand.order),
                  graph_sha256=sha(normal), input_sha256=sha(inp.read_bytes()), output_sha256=sha(out.read_bytes()),
                  true_before=e0, true_after=e1, observed_gradient_calls=calls, local_attempts=len(observed),
                  steps_per_attempt=int(re.search(r"Number of local optimization steps: (\d+)", log)[1]),
                  wall_s=wall, cpu_s=cpu, observer_wall_s=obs_wall, observer_cpu_s=obs_cpu)
    xyz = []
    for prefix, raw, en in zip(("before","after"), (text, out.read_text()), endpoints, strict=True):
        lines = raw.splitlines()
        x=ligand.coords(lines)
        ligand.check_center(x)
        xyz.append(x)
        terms = calculate_components(lines, maps, operator="combined", uncapped_terms=True)
        if abs(en["unbound"]-en["intra"]) > 1e-10 or abs(en["S"]-en["inter"]-en["torsion"]) > 0.0011:
            raise ValueError(f"Printed energy identity failed: {name}")
        result.update({prefix+"_"+k:v for k,v in en.items()})
        result.update({prefix+"_map_"+k:v for k,v in terms.items()})
        result[prefix+"_S_from_xyz"] = sum(terms.values())
    result["fixed_serial_drift_A"] = np.linalg.norm(xyz[1]-xyz[0])/np.sqrt(len(ligand.order))
    return result


def rectangular(rows):
    columns = dict.fromkeys(k for r in rows for k in r)
    return [dict.fromkeys(columns, "") | r for r in rows]


def consensus_rows(cells):
    pooled = {}
    for cell in cells:
        key = tuple(cell[k] for k in ("panel", "arm", "seed"))
        pooled.setdefault(key, []).append(cell)
    consensus = []
    for key, group in pooled.items():
        assert {c["receiver"] for c in group} == {"1FKG", "1FKH"} and len(group) == 2
        row = dict(zip(("panel", "arm", "seed"), key, strict=True))
        row["status"] = "ok" if all(c["status"] == "ok" for c in group) else "graph_mismatch"
        if row["status"] == "ok":
            for when in ("before", "after"):
                for term in ("substitution", "placement", "net"):
                    col = when+"_"+term
                    row[col] = sum(c[col] for c in group)/2
        consensus.append(row)
    return rectangular(consensus)


def main(native, observer, graph, source, out):
    tic, cpu0 = time.perf_counter(), time.process_time()
    subprocess0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    out.mkdir(exist_ok=True, parents=True)
    with (source/"manifold_poses.csv").open() as f:
        poses = list(csv.DictReader(f))
    geom = np.load(source/"manifold_xyz.npz")
    starts = {tuple(r[k] for k in KEY): geom["xyz"][i] for i,r in enumerate(poses) if r["rank"]=="1"}
    assert len(starts)==26 and len(poses)==len(geom["xyz"])
    results, indices, cells = {}, [], []
    with tempfile.TemporaryDirectory(prefix="fd_refine_") as temp:
        work = Path(temp)
        ligs, templates, maps = prep(work)
        np.testing.assert_array_equal(ligs[5].order, geom["pdbqt_serials"])
        try:
            for key, pair in pool_sets():
                panel, arm, seed, rec = key
                pool_key = (panel, arm, str(seed), rec)
                legs = {}
                original = {n: ligs[n].coords(models(data)[0]) for n,data in zip((5,6),pair,strict=True)}
                for label, n, x in (("native5",5,original[5]), ("feasible5",5,starts[pool_key]), ("native6",6,original[6])):
                    text = coordinates(templates[n],ligs[n],x)
                    job = f"{rec}_{n}_"+sha(text.encode())[:20]
                    if job not in results:
                        results[job] = run_leg(work,job,ligs[n],n,text,rec,maps[rec],(native,observer,graph), not results)
                    row = dict(zip(KEY,key,strict=True))
                    indices.append(row | dict(leg=label,job=job))
                    legs[label] = results[job]
                c = row.copy()
                c["status"] = "ok" if all(leg["status"] == "ok" for leg in legs.values()) else "graph_mismatch"
                if c["status"] != "ok":
                    cells.append(c)
                    continue
                symmetry = np.array(ligs[5].graph.GetSubstructMatches(ligs[5].graph, uniquify=False, useChirality=True))
                for endpoint in ("before","after"):
                    file = "input.pdbqt" if endpoint == "before" else "output.pdbqt"
                    left, right = [ligs[5].coords((work/"evidence"/legs[k]["job"]/file).read_text().splitlines())
                                   for k in ("native5","feasible5")]
                    c[endpoint+"_two5_frame_rmsd_A"] = float(np.sqrt(((left[None]-right[symmetry])**2).sum(axis=2).mean(axis=1).min()))
                    scores = [legs[k][endpoint+"_S_from_xyz"] for k in ("native5","feasible5","native6")]
                    c.update({endpoint+"_substitution":scores[2]-scores[1],
                              endpoint+"_placement":scores[1]-scores[0], endpoint+"_net":scores[2]-scores[0]})
                cells.append(c)
        finally:
            with tarfile.open(out/"refine_raw.tar.xz","w:xz") as tar:
                tar.add(work/"evidence",arcname="evidence")
    assert len(indices)==78
    ok = [r for r in results.values() if r["status"] == "ok"]
    csv_out(out/"refine_jobs.csv",rectangular(results.values()))
    csv_out(out/"refine_legs.csv",indices)
    csv_out(out/"refine_cells.csv",rectangular(cells))
    csv_out(out/"refine_consensus.csv", consensus_rows(cells))
    costs = {k: sum(r[k] for r in ok) for k in ("wall_s","cpu_s","observer_wall_s","observer_cpu_s")}
    subcpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    info = dict(sha256_native=sha(native.read_bytes()),sha256_observer=sha(observer.read_bytes()),sha256_graph=sha(graph.read_bytes()),
                sha256_manifold=sha((source/"manifold_xyz.npz").read_bytes()), sha256_raw=sha((out/"refine_raw.tar.xz").read_bytes()),
                n_cells=len(cells), requested_legs=len(indices), unique_jobs=len(results), successful=len(ok),
                calls=sum(r["observed_gradient_calls"] for r in ok),
                max_before_map_abs_kcal_mol=max(abs(r["before_S"]-r["before_S_from_xyz"]) for r in ok),
                max_after_map_abs_kcal_mol=max(abs(r["after_S"]-r["after_S_from_xyz"]) for r in ok),
                full_wall_s=time.perf_counter()-tic, full_cpu_s=time.process_time()-cpu0+subcpu.ru_utime+subcpu.ru_stime-subprocess0.ru_utime-subprocess0.ru_stime,
                **costs)
    (out/"refine_receipt.json").write_text(json.dumps(info,indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("native","observer","graph","source","out"):
        parser.add_argument(arg,type=lambda s:Path(s).resolve())
    a=parser.parse_args()
    main(a.native,a.observer,a.graph,a.source,a.out)
