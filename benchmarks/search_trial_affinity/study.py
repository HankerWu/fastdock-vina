"""Run or independently audit the four registered arms, reusing tree-identical raw baselines."""

import argparse
import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "cumulative/affinity"
sys.path.insert(0, str(SOURCE))
from run import digest, prepare
from evaluate import analyze, write
from verify import independent
from scipy.stats import spearmanr, kendalltau


def setup(work):
    work.mkdir(parents=True, exist_ok=False)
    prepare(work)
    historical = SOURCE / "results"
    hashes = {
        line.split()[1].lstrip("*"): line.split()[0]
        for line in (historical / "SHA256SUMS").read_text().splitlines()
    }
    hashes["common_cost.json"] = (
        "e4b0e1d4445e13623983aa004d2a122b92f797b3c89fe8e8e385267536186e4b"
    )
    for name in ("raw.tar.xz", "common_cost.json"):
        assert digest(historical / name) == hashes[name]
    with tarfile.open(historical / "raw.tar.xz") as archive:
        for m in archive:
            if any(
                m.name.startswith("outputs/" + e + "/")
                for e in ("official", "parent", "candidate")
            ):
                archive.extract(m, work, filter="data")
    (work / "outputs/candidate").rename(work / "outputs/stable")
    prior = json.loads((historical / "common_cost.json").read_text())
    prior["attempts"] = [
        dict(a, engine="stable" if a["engine"] == "candidate" else a["engine"])
        for a in prior["attempts"]
        if a["engine"] != "stable"
    ]
    prior["final"]["stable"] = prior["final"].pop("candidate")
    (work / "baseline_cost.json").write_text(json.dumps(prior, indent=2) + "\n")


def experiments(work, binary, jobs):
    for stages in (["nominal", "low", "matched"], ["default"]):
        subprocess.run(
            [
                sys.executable,
                SOURCE / "run.py",
                work,
                "--binary",
                binary,
                "--engine",
                "candidate",
                "--stage",
                *stages,
                "--jobs",
                str(jobs),
            ],
            check=True,
        )
    subprocess.run(
        [
            sys.executable,
            SOURCE / "common_cost.py",
            work,
            "--candidate",
            binary,
            "--replace-candidate",
            work / "baseline_cost.json",
        ],
        check=True,
    )


def chemical_series(out):
    groups = defaultdict(list)
    with (out / "ligands.csv").open() as f:
        for r in csv.DictReader(f):
            if r["stage"] == "nominal":
                groups[
                    tuple(
                        r[k]
                        for k in ("target", "assay", "receptor", "engine", "scaffold")
                    )
                ].append(r)
    rows = []
    for key, available in sorted(groups.items()):
        valid = [r for r in available if r["score"]]
        x = [float(r["value_nM"]) for r in valid]
        y = [float(r["score"]) for r in valid]
        eligible = [
            (j, k)
            for j, k in combinations(range(len(x)), 2)
            if max(x[j], x[k]) >= 3 * min(x[j], x[k])
        ]
        credit = sum(
            0.5 if y[j] == y[k] else float((x[j] < x[k]) == (y[j] < y[k]))
            for j, k in eligible
        )
        diverse = len(set(x)) > 1 and len(set(y)) > 1
        rows.append(
            dict(
                zip(("target", "assay", "receptor", "engine", "scaffold"), key),
                planned=len(available),
                observed=len(valid),
                pairs_3x=len(eligible),
                correct_3x=credit,
                rho=spearmanr(x, y).statistic if diverse else None,
                tau_b=kendalltau(x, y).statistic if diverse else None,
            )
        )
    write(out / "scaffolds.csv", rows)


def analyze_and_audit(work, out):
    analyze(work, out)
    chemical_series(out)
    audit = independent(work, out)
    audit["new_search_receipts"] = len(
        list((work / "outputs/candidate").glob("*/*.json"))
    )
    audit["prepared_input_hashes"] = json.loads((work / "inputs.json").read_text())[
        "files"
    ]
    (out / "verification.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: v for k, v in audit.items() if not isinstance(v, dict)}))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("step", choices=["setup", "run", "collect", "replay"])
    p.add_argument("work", type=Path)
    p.add_argument("--binary", type=Path)
    p.add_argument("--jobs", type=int, default=2)
    p.add_argument("--artifacts", type=Path, default=HERE / "results")
    args = p.parse_args()
    work, out = args.work.resolve(), args.artifacts.resolve()
    if args.step == "setup":
        setup(work)
    elif args.step == "run":
        experiments(work, args.binary.resolve(strict=True), args.jobs)
    elif args.step == "collect":
        analyze_and_audit(work, out)
        for name in ("common_cost.json", "inputs.json"):
            shutil.copyfile(work / name, out / name)
        with tarfile.open(out / "raw_candidate.tar.xz", "w:xz") as tar:
            tar.add(work / "outputs/candidate", arcname="outputs/candidate")
    else:
        setup(work)
        with tarfile.open(out / "raw_candidate.tar.xz") as tar:
            tar.extractall(work, filter="data")
        shutil.copyfile(out / "common_cost.json", work / "common_cost.json")
        fresh = work / "recomputed"
        analyze_and_audit(work, fresh)
        for file in fresh.iterdir():
            assert file.read_bytes() == (out / file.name).read_bytes(), file.name


if __name__ == "__main__":
    main()
