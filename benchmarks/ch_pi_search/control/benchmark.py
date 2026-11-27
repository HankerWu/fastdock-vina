"""Official Vina version and increased native-search budget controls for the GAA factorial."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify import geometry, point
from study import digest, timed_process

KEY = ["query", "receptor", "seed"]
OFFICIAL = "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"


def run(gaa, original, pro, official, out):
    if "FASTDOCK_RESEARCH_CHPI" in os.environ:
        raise ValueError("Unexpected research grid in native controls")
    out.mkdir(parents=True, exist_ok=True)
    pools = out / "pools"
    pools.mkdir(exist_ok=True)
    engines = {"official": official, "pro": pro}
    versions = {k: {"sha256": digest(v), "version": subprocess.check_output([v, "--version"], text=True).strip()}
                for k, v in engines.items()}
    if versions["official"]["sha256"] != OFFICIAL or "1.2.7" not in versions["official"]["version"]:
        raise ValueError("Official binary differs from prior GBA1 comparator")
    if versions["pro"]["sha256"] != json.loads((original / "binary.json").read_text())["sha256"]:
        raise ValueError("Pro binary differs from the four-arm comparator")
    (out / "binaries.json").write_text(json.dumps(versions, indent=2) + "\n")
    jobs = pd.read_csv(original / "costs.csv").query("cap==12000 and search=='native'")
    if len(jobs) != 12 or jobs.duplicated(KEY).any():
        raise ValueError("Expected twelve paired conditions")
    completed = []
    control_arms = [("pro", 48000), ("official", 12000), ("official", 24000)]
    for j, row in enumerate(jobs.itertuples()):
        start = j % len(control_arms)
        for engine, cap in control_arms[start:] + control_arms[:start]:
            args = shlex.split(row.args)
            args[args.index("--max_evals") + 1] = str(cap)
            if engine == "official":
                for option in ("--pose_selection", "--refine_candidates"):
                    i = args.index(option)
                    del args[i:i+2]
                args.remove("--search_counters")
            label = f"{row.query}_{row.receptor}_{row.seed}_{engine}_{cap}"
            pose = pools / (label + ".pdbqt")
            log = pools / (label + ".log")
            archive_args = shlex.join(args)
            args[args.index("--out") + 1] = str(pose)
            timing = timed_process([str(engines[engine]), *args], gaa, log)
            traces = re.findall(r"hunt_evals=(\d+).*?second_evals=(\d+)", log.read_text())
            explicit = re.findall(r"explicit_evals=(\d+)", log.read_text())
            if engine == "pro" and (len(traces) != 4 or len(explicit) != 1):
                raise ValueError(f"Missing search counters in {label}")
            completed.append(dict(query=row.query, receptor=row.receptor, seed=row.seed, engine=engine, cap=cap,
                                  label=label, args=archive_args, sha256=digest(pose), **timing,
                                  search_evals=sum(int(a)+int(b) for a,b in traces) if traces else "",
                                  explicit_evals=sum(map(int, explicit)) if explicit else ""))
    pd.DataFrame(completed).to_csv(out / "costs.csv", index=False)
    with tarfile.open(out / "outputs.tar.xz", "w:xz") as archive:
        for path in sorted(pools.glob("*")):
            info = archive.gettarinfo(path, arcname=path.name)
            info.mtime = info.uid = info.gid = 0
            info.uname = info.gname = ""
            with path.open("rb") as data:
                archive.addfile(info, data)


def assess(gaa, original, out):
    costs = pd.read_csv(out / "costs.csv")
    modes, cells = [], []
    with tarfile.open(out / "outputs.tar.xz") as archive:
        for r in costs.itertuples():
            data = archive.extractfile(r.label + ".pdbqt").read()
            if hashlib.sha256(data).hexdigest() != r.sha256:
                raise ValueError("Changed pose output")
            refs, order = geometry(gaa, "GAA", r.query, r.receptor)
            case = dict(query=r.query, receptor=r.receptor, seed=r.seed, engine=r.engine, cap=r.cap)
            pool = []
            for block in data.decode().split("ENDMDL"):
                if not block.strip():
                    continue
                rank = int(block.strip().splitlines()[0].split()[1])
                xyz = {int(l[6:11]): point(l) for l in block.splitlines() if l.startswith(("ATOM  ", "HETATM"))}
                rmsd = np.sqrt(np.square(np.array([xyz[i] for i in order]) - refs).sum(axis=-1).mean(axis=-1)).min()
                score = float(next(l.split()[3] for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")))
                pool.append(dict(**case, rank=rank, rmsd=rmsd, score=score))
            if not pool or [p["rank"] for p in pool] != list(range(1,len(pool)+1)):
                raise ValueError(f"Empty or malformed output {r.label}")
            modes.extend(pool)
            errors = [p["rmsd"] for p in pool]
            cells.append(dict(**case, modes=len(pool), top1=errors[0], best5=min(errors[:5]),
                              pool_min=min(errors), regret=errors[0]-min(errors)))
    mode_table = pd.DataFrame(modes)
    cell_table = pd.DataFrame(cells)
    mode_table.to_csv(out / "modes.csv", index=False)
    cell_table.to_csv(out / "cells.csv", index=False)
    orig_cells = pd.read_csv(original / "cells.csv").query("search==rerank").copy()
    orig_cells["engine"] = orig_cells.search.map(dict(native="pro", field="joint"))
    orig_cells.rename(columns={"output_count":"modes"}, inplace=True)
    all_cells = pd.concat([orig_cells[cell_table.columns], cell_table])
    orig_cost = pd.read_csv(original / "costs.csv")
    orig_cost["engine"] = orig_cost.search.map(dict(native="pro",field="joint"))
    orig_cost["ranking_s"] = np.where(orig_cost.engine=="joint", orig_cost.analyze_pool_s, 0.)
    costs["ranking_s"] = 0.
    all_cost = pd.concat([orig_cost, costs], ignore_index=True)
    summary = []
    for (engine, cap), group in all_cells.groupby(["engine", "cap"]):
        price = all_cost[(all_cost.engine == engine) & (all_cost.cap == cap)]
        row = dict(engine=engine, cap=cap, jobs=len(group), modes_total=group.modes.sum(),
                   wall_s=price.wall_s.sum(), pool_analysis_s=price.ranking_s.sum(),
                   cpu_s=price.cpu_s.sum(), median_peak_rss_kib=price.peak_rss_kib.median(),
                   search_evals=price.search_evals.sum(min_count=1), explicit_evals=price.explicit_evals.sum(min_count=1))
        for metric in ["top1","best5","pool_min","regret"]:
            x = group[metric]
            row.update({metric+"_"+stat: v for stat, v in dict(
                mean=x.mean(), median=x.median(), q25=x.quantile(.25), q75=x.quantile(.75),
                p90=x.quantile(.9), p95=x.quantile(.95)).items()})
        summary.append(row)
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False)
    pairs = []
    for cap in (12000,24000):
        joint = all_cells[(all_cells.engine=="joint") & (all_cells.cap==cap)].set_index(KEY)
        for engine, budget in (("pro",48000),("official",cap)):
            base = all_cells[(all_cells.engine==engine)&(all_cells.cap==budget)].set_index(KEY)
            delta = joint[["top1","pool_min","regret"]]-base[["top1","pool_min","regret"]]
            if len(delta)!=12 or delta.isna().any().any():
                raise ValueError("Unpaired control")
            for ix, row in delta.iterrows():
                pairs.append(dict(zip(KEY, ix)) | dict(joint_cap=cap, comparator=engine, comparator_cap=budget) | row.to_dict())
    pd.DataFrame(pairs).to_csv(out / "paired_differences.csv", index=False)


def verify_native(gaa, out, executable):
    successful = 0
    if "FASTDOCK_RESEARCH_CHPI" in os.environ:
        raise ValueError("Unexpected research grid in stable Pro control")
    with TemporaryDirectory() as tmp:
        pose = Path(tmp) / "pose.pdbqt"
        log = Path(tmp) / "run.log"
        for row in pd.read_csv(out / "costs.csv").query("engine == 'pro'").itertuples():
            args = shlex.split(row.args)
            args[args.index("--out") + 1] = str(pose)
            timed_process([str(executable), *args], gaa, log)
            if digest(pose) != row.sha256:
                raise ValueError(f"Unpatched release output differs at {row.label}")
            successful += 1
    if successful != 12:
        raise ValueError("Expected twelve increased-budget native comparisons")
    data = dict(paired_jobs=successful, byte_identical=successful, binary_sha256=digest(executable))
    (out / "native_parity.json").write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("gaa","original","out"):
        p.add_argument(key,type=lambda s:Path(s).resolve())
    p.add_argument("--binaries", nargs=2, metavar=("PRO", "OFFICIAL"), type=lambda s:Path(s).resolve())
    p.add_argument("--verify-native", type=lambda s:Path(s).resolve(),
                   help="Replay increased-budget results using a separate unpatched stable Pro binary")
    a=p.parse_args()
    if a.binaries:
        run(a.gaa, a.original, *a.binaries, a.out)
    assess(a.gaa, a.original, a.out)
    if a.verify_native:
        verify_native(a.gaa, a.out, a.verify_native)
