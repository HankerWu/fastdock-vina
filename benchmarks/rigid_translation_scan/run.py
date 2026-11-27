"""Run the registered lattice and equal-count controls, then evaluate held-out labels."""

import argparse
import csv
import itertools
import json
import re
import resource
import subprocess
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

from scan import read_ligand, receptor_maps, translate_search, write_start

HERE = Path(__file__).resolve().parent

spec = spec_from_file_location(
    "paired_existing", HERE.parent / "paired_crossdocking/run.py"
)
paired = module_from_spec(spec)
spec.loader.exec_module(paired)
CASES = {
    "5C7N": ["5EPS", "5EPR", "5FG5"],
    "4OD0": ["3WKE", "4X6Y", "3ANT"],
    "1A28": ["1SQN"],
    "1SQN": ["1A28"],
    "2P54": ["1K7L"],
    "1K7L": ["2P54"],
}


def save(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def rms(coords, refs):
    return np.sqrt(np.square(refs - coords).sum(axis=-1).mean(axis=-1)).min()


def parse_pose(path, serials):
    lines = path.read_text().splitlines()
    atoms = {
        int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
        for l in lines
        if l.startswith(("ATOM", "HETATM"))
    }
    return np.array([atoms[i] for i in serials])


def cpu():
    ch = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + ch.ru_utime + ch.ru_stime


def main():
    p = argparse.ArgumentParser()
    for name in ("dev", "hold", "work", "probe", "vina"):
        p.add_argument("--" + name, required=True, type=Path)
    p.add_argument("--seed", type=int, nargs="+", default=[43103, 43117])
    p.add_argument("--receivers", nargs="+", default=list(CASES))
    p.add_argument(
        "--arms",
        nargs="+",
        choices=["lattice", "random"],
        default=["lattice", "random"],
    )
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    devmeta = json.loads((a.dev / "metadata.json").read_text())
    poses = []
    costs = []
    all_ranks = []
    cells = []
    for rec in a.receivers:
        dev = rec in ("5C7N", "4OD0")
        base = a.dev if dev else a.hold
        meta = (
            devmeta
            if dev
            else {
                entry: json.loads((base / entry / "preparation.json").read_text())
                for entry in CASES[rec] + [rec]
            }
        )
        center = np.array(meta[rec]["center"])
        start = time.perf_counter()
        c = cpu()
        fields, origin, spacing = receptor_maps(
            a.vina, base / rec / "receptor.pdbqt", center, a.work / "maps" / rec
        )
        costs.append(
            {
                "receiver": rec,
                "query": "",
                "seed": "",
                "arm": "receiver_maps",
                "probes": "",
                "lookups": "",
                "refines": "",
                "gradients": "",
                "wall_s": time.perf_counter() - start,
                "cpu_s": cpu() - c,
            }
        )
        for query in CASES[rec]:
            start = time.perf_counter()
            c = cpu()
            ligand_path = base / query / "ligand.pdbqt"
            ligand = read_ligand(a.probe, ligand_path)
            costs.append(
                {
                    "receiver": rec,
                    "query": query,
                    "seed": "",
                    "arm": "ligand_types",
                    "probes": "",
                    "lookups": "",
                    "refines": "",
                    "gradients": "",
                    "wall_s": time.perf_counter() - start,
                    "cpu_s": cpu() - c,
                }
            )
            for seed, arm in itertools.product(a.seed, a.arms):
                start = time.perf_counter()
                c = cpu()
                dest = a.work / "runs" / f"{rec}_{query}_{seed}_{arm}"
                dest.mkdir(parents=True, exist_ok=True)
                proposals, starts, probes = translate_search(
                    fields, origin, spacing, center, ligand, seed, arm
                )
                paths = []
                for k, (_e, i, coords) in enumerate(starts):
                    path = dest / f"rank{k + 1}.pdbqt"
                    write_start(path, *ligand[:2], coords)
                    paths.append(path)
                command = [
                    str(a.probe),
                    "refine",
                    str(base / rec / "receptor.pdbqt"),
                    str(ligand_path),
                    *map(str, center),
                    *map(str, paths),
                ]
                output = subprocess.run(
                    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
                )
                (dest / "refine.log").write_text(output.stdout)
                if output.returncode:
                    raise RuntimeError(str(dest) + " " + output.stdout[-1200:])
                refined = [
                    l.split()
                    for l in output.stdout.splitlines()
                    if l.startswith("result ")
                ]
                counts = list(
                    map(int, re.findall(r"gradient_evals (\d+)", output.stdout))
                )
                if len(counts) != len(paths) or len(refined) != len(paths):
                    raise ValueError("incomplete refinements")
                costs.append(
                    {
                        "receiver": rec,
                        "query": query,
                        "seed": seed,
                        "arm": arm,
                        "probes": probes,
                        "lookups": probes * len(ligand[3]),
                        "refines": len(paths),
                        "gradients": sum(counts),
                        "wall_s": time.perf_counter() - start,
                        "cpu_s": cpu() - c,
                    }
                )
                # Experimental coordinates enter only after every start was scored and refined.
                rotation, shift, _, _ = paired.alignment(meta[query], meta[rec])
                refs = np.asarray(meta[query]["mapped"]) @ rotation + shift
                serials = meta[query]["serial"]
                pool = []
                for (e, i, _), path, line, n in zip(starts, paths, refined, counts):
                    before = parse_pose(path, serials)
                    after = parse_pose(Path(str(path) + ".refined"), serials)
                    pool.append(
                        {
                            "receiver": rec,
                            "query": query,
                            "seed": seed,
                            "arm": arm,
                            "slot": path.stem,
                            "orientation": i,
                            "coarse_energy": e,
                            "explicit_before": float(line[2]),
                            "explicit_after": float(line[3]),
                            "before_rmsd": rms(before, refs),
                            "after_rmsd": rms(after, refs),
                            "gradients": n,
                        }
                    )
                # Native-style minimum spacing on indexed chemical heavy atoms; rank and keep original coordinates.
                kept = []
                kept_xyz = []
                for row in sorted(pool, key=lambda r: r["explicit_after"]):
                    xyz = parse_pose(dest / (row["slot"] + ".pdbqt.refined"), serials)
                    if all(
                        np.square(xyz - prior).sum(axis=1).mean() >= 1
                        for prior in kept_xyz
                    ):
                        kept.append(row)
                        kept_xyz.append(xyz)
                for idx, row in enumerate(kept):
                    all_ranks.append(dict(**row, rank=idx + 1))
                poses.extend(pool)
                cells.append(
                    {
                        "receiver": rec,
                        "query": query,
                        "stage": "development" if dev else "transfer",
                        "seed": seed,
                        "arm": arm,
                        "reported": len(kept),
                        "top": kept[0]["after_rmsd"],
                        "top5": min(r["after_rmsd"] for r in kept[:5]),
                        "oracle": min(r["after_rmsd"] for r in kept),
                        "best_before": min(r["before_rmsd"] for r in pool),
                    }
                )
            for name, rows in [
                ("poses", poses),
                ("ranked", all_ranks),
                ("cost", costs),
                ("cells", cells),
            ]:
                save(a.work / (name + ".csv"), rows)
            print(rec, query, flush=True)


if __name__ == "__main__":
    main()
