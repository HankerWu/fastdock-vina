"""Paired fixed-budget native Late docking from a repeated or diversified input."""

import argparse
import csv
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

SEEDS = (3141593, 2718281, 1732051, 2236067)


def pools(work, data, seed):
    for rec in data["receptors"]:
        for state, slot in [
            ("base", 0),
            ("base", 1),
            ("alternate", 1),
            ("sham", 1),
            ("first", 0),
            ("first", 1),
        ]:
            path = (
                work / "pools" / data["target"] / rec / f"{seed}_{slot}_{state}.pdbqt"
            )
            yield rec, state, slot, path


def dock(args):
    binary, folder, rec, state, seed, output, max_evals = args
    output.parent.mkdir(exist_ok=True, parents=True)
    tmp = output.with_suffix(".pending")
    cmd = [
        str(binary),
        "--receptor",
        str(folder / rec / "receptor.pdbqt"),
        "--ligand",
        str(folder / f"{state}.pdbqt"),
        "--config",
        str(folder / rec / "box.txt"),
        "--cpu",
        "1",
        "--seed",
        str(seed),
        "--exhaustiveness",
        "3",
        "--max_evals",
        str(max_evals),
        "--pose_selection",
        "late",
        "--num_modes",
        "9",
        "--refine_candidates",
        "12",
        "--energy_range",
        "100",
        "--verbosity",
        "0",
        "--out",
        str(tmp),
    ]
    out = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if out.returncode or not tmp.exists() or "ENDMDL" not in tmp.read_text():
        raise ValueError(f"{output}: {out.stdout} {out.stderr}")
    tmp.replace(output)


def models(path, data, rec, state):
    ref = np.array(data["receptors"][rec]["reference"])
    permutations = np.array(data["permutations"])
    result = []
    for index, model in enumerate(path.read_text().split("ENDMDL")[:-1], 1):
        positions = {}
        atom_map = {}
        for line in model.splitlines():
            if line.startswith("REMARK INTER + INTRA:"):
                energy = float(line.split()[-1])
            elif line.startswith("REMARK VINA RESULT:"):
                naive = float(line.split()[3])
            elif line.startswith("REMARK INDEX MAP "):
                v = list(map(int, line.split()[3:]))
                atom_map.update(zip(v[::2], v[1::2], strict=True))
            elif line.startswith(("ATOM", "HETATM")):
                positions[int(line[6:11])] = [
                    float(line[k : k + 8]) for k in (30, 38, 46)
                ]
        xyz = np.array([positions[atom_map[i + 1]] for i in range(len(ref))])
        diff = xyz - ref[permutations]
        rmsd = float(np.sqrt((diff**2).sum(-1).mean(-1)).min())
        result.append(
            {
                "model": index,
                "state": state,
                "common_energy": energy,
                "native_score": naive,
                "rmsd": rmsd,
            }
        )
    return result


def write_csv(out, rows):
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def evaluate(work, items, out):
    out.mkdir(parents=True, exist_ok=True)
    details = []
    cases = []
    for _, data in items:
        for seed in SEEDS:
            arms = {}
            for rec, state, slot, path in pools(work, data, seed):
                if state != "base" and not data["states"][state + "_distinct"]:
                    path = path.with_name(f"{seed}_{slot}_base.pdbqt")
                entries = models(path, data, rec, state)
                for entry in entries:
                    details.append(
                        dict(
                            target=data["target"],
                            receptor=rec,
                            seed=seed,
                            slot=slot,
                            **entry,
                        )
                    )
                arms.setdefault(rec, {})[slot, state] = entries
            for rec, parts in arms.items():
                info = data["receptors"][rec]
                for arm, label in [
                    ("repeat", "base"),
                    ("ring", "alternate"),
                    ("sham", "sham"),
                    ("first", "first"),
                ]:
                    merged = (
                        parts[0, "first" if arm == "first" else "base"]
                        + parts[1, label]
                    )
                    top = min(merged, key=lambda v: v["common_energy"])
                    best = min(v["rmsd"] for v in merged)
                    naive = min(merged, key=lambda v: v["native_score"])
                    energies = data["states"]["energies_MMFF94s"]
                    costs = {
                        "base": 0,
                        "first": data["states"]["first_energy_MMFF94s"] - energies[0],
                        "alternate": energies[1] - energies[0],
                        "sham": data["states"]["sham_energy_MMFF94s"] - energies[0],
                    }
                    gas = min(
                        merged, key=lambda v: v["common_energy"] + costs[v["state"]]
                    )
                    cases.append(
                        {
                            "target": data["target"],
                            "receptor": rec,
                            "kind": info["kind"],
                            "seed": seed,
                            "arm": arm,
                            "ring_distinct": data["states"]["alternate_distinct"],
                            "top_state": top["state"],
                            "top_rmsd": top["rmsd"],
                            "naive_rmsd": naive["rmsd"],
                            "gas_offset_rmsd": gas["rmsd"],
                            "best_rmsd": best,
                            "top_hit": int(top["rmsd"] <= 2.0),
                            "pool_hit": int(best <= 2.0),
                            "n_poses": len(merged),
                        }
                    )
    write_csv(out / "poses.csv", details)
    write_csv(out / "cases.csv", cases)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--target", nargs="*")
    parser.add_argument(
        "--evaluate-existing",
        action="store_true",
        help="Evaluate the complete archived pools without new docking",
    )
    parser.add_argument("--max-evals", type=int, default=12000)
    args = parser.parse_args()
    if not args.evaluate_existing and args.binary is None:
        parser.error("--binary is required when running new docking searches")
    tasks = []
    items = []
    for f in sorted((args.work / "prepared").glob("*/panel.json")):
        data = json.loads(f.read_text())
        folder = f.parent
        if args.target and data["target"] not in args.target:
            continue
        items.append((folder, data))
        for seed in SEEDS:
            for rec, state, slot, path in pools(args.work, data, seed):
                if state != "base" and not data["states"][state + "_distinct"]:
                    continue
                if args.evaluate_existing:
                    continue
                tasks.append(
                    (
                        args.binary.resolve(),
                        folder,
                        rec,
                        state,
                        seed + 1009 * slot,
                        path,
                        args.max_evals,
                    )
                )
    if not items:
        raise ValueError("No prepared panels")
    with ThreadPoolExecutor(args.workers) as executor:
        list(executor.map(dock, tasks))
    evaluate(args.work, items, args.out)
    print("Docked", len(tasks), "processes across", len(items), "targets.")
