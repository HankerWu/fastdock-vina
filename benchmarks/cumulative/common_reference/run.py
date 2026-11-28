"""Apply a common reference to both blind seeds, preserving all pool failures."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
POSE = HERE.parents[1] / "geometry_comparison/upstream/four_versions"
sys.path.insert(0, str(POSE))
import run_panel as four

STAGES = {"cognate", "reciprocal", "external"}
TOTAL = re.compile(r"^Estimated Free Energy of Binding\s*:\s*([\d.\-]+)", re.M)
INTER = re.compile(r"^\(1\) Final Intermolecular Energy\s*:\s*([\d.\-]+)", re.M)
INTRA = re.compile(r"^\(2\) Final Total Internal Energy\s*:\s*([\d.\-]+)", re.M)


def csvrows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def read_models(source, ranks, cell, work):
    pool = source / (cell["key"] + ".pdbqt")
    if four.digest(pool) != cell["pool_sha256"]:
        raise ValueError(f"changed pose file {pool}")
    models = []
    for block in pool.read_text().split("ENDMDL"):
        lines = block.strip().splitlines()
        if not lines:
            continue
        rank = int(lines[0].split()[1])
        path = work / f"{cell['seed']}-{rank}.pdbqt"
        path.write_text("\n".join(lines[1:]) + "\n")
        r = ranks[cell["key"], rank]

        def energy(prefix):
            return float(
                next(
                    line.split(":", 1)[1].split()[0]
                    for line in lines
                    if line.startswith(prefix)
                )
            )

        printed = energy("REMARK VINA RESULT:")
        if printed != float(r["score"]):
            raise ValueError("released score mismatch")
        models.append(
            dict(
                seed=int(cell["seed"]),
                rank=rank,
                rmsd=float(r["rmsd"]),
                file=path.name,
                original=printed,
                archived=energy("REMARK INTER + INTRA:"),
                unbound=energy("REMARK UNBOUND:"),
                sha256=four.digest(path),
            )
        )
    if len(models) != int(cell["modes"]):
        raise ValueError("incomplete pose pool")
    return models


def batch_scores(text):
    parsed = {}
    for block in text.split("Scoring ligand: ")[1:]:
        name = Path(block.splitlines()[0]).name
        found = [r.search(block) for r in (TOTAL, INTER, INTRA)]
        if all(found):
            parsed[name] = [float(v.group(1)) for v in found]
    return parsed


def score_pair(group, cells, ranks, raw, fixtures, result, binary, probe):
    cohort, arm, query, receiver = group
    work = result / ("-".join(group))
    work.mkdir()
    info = dict(
        cohort=cohort,
        arm=arm,
        query=query,
        receiver=receiver,
        target=cells[0]["target"],
        search_cpu_s=sum(float(c["cpu_s"]) for c in cells),
        search_wall_s=sum(float(c["wall_s"]) for c in cells),
    )
    if sorted(int(c["seed"]) for c in cells) != [85231, 85237]:
        raise ValueError("missing scheduled seed row")
    if any(c["status"] != "ok" for c in cells):
        return dict(info, status="unavailable_search"), [], []
    models = []
    for cell in sorted(cells, key=lambda c: c["seed"]):
        models.extend(read_models(raw, ranks, cell, work))
    meta_file = "preparation.json" if cohort == "reciprocal" else "metadata.json"
    center = json.loads((fixtures / receiver / meta_file).read_text())["center"]
    receptor, ligand = (
        fixtures / receiver / "receptor.pdbqt",
        fixtures / query / "ligand.pdbqt",
    )
    files = [str(work / m["file"]) for m in models]
    cmd = [
        str(binary),
        "--receptor",
        str(receptor),
        "--cpu",
        "1",
        "--seed",
        "85231",
        "--score_only",
        "--spacing",
        "0.375",
        "--dir",
        str(work),
    ]
    for axis, v in zip("xyz", center):
        cmd += ["--center_" + axis, str(v), "--size_" + axis, "24"]
    stages = [
        ("batch", cmd + ["--unbound_energy", "0", "--batch"] + files),
        (
            "fixed",
            [str(probe), str(receptor), str(ligand)] + list(map(str, center)) + files,
        ),
    ]
    if cohort == "external":
        stages.append(("minus3", cmd + ["--unbound_energy", "-3", "--batch"] + files))
    costs = []
    for stage, args in stages:
        log = work / (stage + ".log")
        cost = four.execute(args, log, 180)
        costs.append(dict(info, stage=stage, **cost))
        output = log.read_text()
        if stage == "fixed":
            for line in output.splitlines():
                if line.startswith("POSE "):
                    _, idx, energy, inter, intra = line.split()
                    models[int(idx)].update(
                        fixed=float(energy),
                        fixed_inter=float(inter),
                        fixed_intra=float(intra),
                    )
        else:
            by_name = batch_scores(output)
            for m in models:
                if m["file"] in by_name:
                    energy, inter, intra = by_name[m["file"]]
                    m.update(
                        {
                            stage: energy,
                            stage + "_inter": inter,
                            stage + "_intra": intra,
                        }
                    )
    for m in models:
        m.update({k: info[k] for k in ("arm", "cohort", "target", "query", "receiver")})
    info.update(modes=len(models), oracle=min(m["rmsd"] for m in models))
    for name in ["original", "archived"] + [s for s, _ in stages]:
        if any(name not in m for m in models):
            continue
        minimum = min(m[name] for m in models)
        ties = [m for m in models if m[name] == minimum]
        winner = ties[0]
        info.update(
            {
                name + "_winner": winner["file"],
                name + "_ties": len(ties),
                name + "_A": winner["rmsd"],
                name + "_regret": winner["rmsd"] - info["oracle"],
            }
        )
    if "batch_A" in info:
        info["max_within_seed_shift_spread"] = max(
            max(ds) - min(ds)
            for ds in (
                [m["batch"] - m["original"] for m in models if m["seed"] == s]
                for s in (85231, 85237)
            )
        )
        info["within_seed_changed_top"] = sum(
            min((m for m in models if m["seed"] == s), key=lambda m: m["batch"])["rank"]
            != 1
            for s in (85231, 85237)
        )
    return (
        dict(
            info,
            status=(
                "ok" if all(c["status"] == "ok" for c in costs) else "partial_rescore"
            ),
        ),
        models,
        costs,
    )


def write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, keys, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "fixtures", "result", "binary", "probe"):
        p.add_argument(name, type=Path)
    p.add_argument(
        "--fixture-subdirectories",
        action="store_true",
        help="look up fixtures by cohort",
    )
    a = p.parse_args()
    src, out = a.source.resolve(), a.result.resolve()
    out.mkdir(parents=True, exist_ok=False)
    pairs = defaultdict(list)
    for cell in csvrows(src / "cells.csv"):
        if cell["cohort"] in STAGES:
            pairs[
                tuple(cell[k] for k in ("cohort", "arm", "query", "receiver"))
            ].append(cell)
    ranks = {(r["key"], int(r["rank"])): r for r in csvrows(src / "ranks.csv")}
    summaries, scores, costs = [], [], []
    for key, pair in sorted(pairs.items()):
        fx = a.fixtures / key[0] if a.fixture_subdirectories else a.fixtures
        s, m, c = score_pair(
            key,
            pair,
            ranks,
            src / "raw",
            fx.resolve(),
            out,
            a.binary.resolve(),
            a.probe.resolve(),
        )
        summaries.append(s)
        scores.extend(m)
        costs.extend(c)
    write(out / "pairs.csv", summaries)
    write(out / "scores.csv", scores)
    write(out / "costs.csv", costs)
    (out / "binaries.json").write_text(
        json.dumps({p.name: four.digest(p) for p in (a.binary, a.probe)}, indent=2)
        + "\n"
    )


if __name__ == "__main__":
    main()
