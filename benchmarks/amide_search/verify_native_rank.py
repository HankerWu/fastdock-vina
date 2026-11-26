"""Check production --rank_amide against 80 matched single-process hCA II docking pools.

First run guidance.py dock and score. Then run this script with
BINARY EXTRACTED_FIXTURES GUIDANCE_WORK NEW_OUTPUT_DIRECTORY [--jobs 4].
The native reference uses each original Vina process's score and unbound energy.
"""

import argparse
import csv
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from guidance import SEARCH_BARRIER, pools
from rank_amide_pool import rank
from run import dock, sources
from select_zinc_pose import read_models


def fields(lines, prefix):
    [record] = [
        line[len(prefix) :].split() for line in lines if line.startswith(prefix)
    ]
    if not all(math.isfinite(float(x)) for x in record):
        raise ValueError(f"Nonfinite {prefix} in MODEL")
    return record


def stable_records(lines):
    # Match both coordinates and the three-decimal native energy components.
    return [
        line
        for line in lines
        if line.startswith(
            (
                "ATOM ",
                "HETATM ",
                "REMARK INTER",
                "REMARK INTRA:",
                "REMARK UNBOUND:",
            )
        )
    ]


def check_one(binary, reference, output, args, sitefile, labels):
    dock(binary, (output, args), rank_amide=sitefile)
    native_text = reference.read_text()
    models = read_models(native_text.splitlines())
    ranked = read_models(output.read_text().splitlines())
    python = rank(args[1].read_text(), native_text, SEARCH_BARRIER, "all-arms")
    expected = {r["model"]: r for r in python}
    if list(ranked) != list(range(1, len(models) + 1)):
        raise ValueError(f"Model count or numbering changed in {output}")

    order, scores, errors, score_errors = [], [], [], []
    for lines in ranked.values():
        printed_native, _, _ = map(float, fields(lines, "REMARK VINA RESULT:"))
        native, native_rank = fields(lines, "REMARK VINA NATIVE:")
        native, native_rank = float(native), int(native_rank)
        score, _, _ = map(float, fields(lines, "REMARK VINA AMIDE RESULT:"))
        [penalty] = map(float, fields(lines, "REMARK VINA AMIDE PENALTY:"))
        r = expected[native_rank]
        if stable_records(lines) != stable_records(models[native_rank]):
            raise ValueError(
                f"Coordinates or native energy components changed in {output}"
            )
        if (
            printed_native != r["native"]
            or abs(native - r["native"]) > 0.000501
            or abs(native + penalty - score) > 1e-9
        ):
            raise ValueError(f"Amide and native energy records disagree in {output}")
        # Python reads coordinates rounded to 0.001 A and energies rounded to
        # 0.001 kcal/mol. C++ ranks before either value is rounded for PDBQT.
        errors.append(abs(penalty - r["amide_penalty"]))
        score_errors.append(abs(score - r["score"]))
        order.append(native_rank)
        scores.append(score)
    if sorted(order) != sorted(models) or scores != sorted(scores):
        raise ValueError(
            f"Native model membership or full-precision sort failed in {output}"
        )
    if max(errors) > 0.025:
        raise ValueError(f"Large torsion-score discrepancy in {output}: {max(errors)}")

    top, python_top = order[0], python[0]["model"]
    if expected[top]["score"] - expected[python_top]["score"] > 2 * (
        max(errors) + 0.000501
    ):
        raise ValueError(f"First model disagrees beyond rounding tolerance in {output}")
    rec, lig, state, seed, slot = reference.stem.split("_")
    hit = [labels[rec, lig, state, seed, slot, n] for n in models]
    return {
        "pool": reference.stem,
        "models": len(models),
        "native_hit": int(hit[0]),
        "ranked_hit": int(labels[rec, lig, state, seed, slot, top]),
        "oracle": int(any(hit)),
        "cross": rec != lig,
        "top_matches_python": (top == python_top),
        "order_matches_python": (order == [r["model"] for r in python]),
        "max_penalty_error": max(errors),
        "max_score_error": max(score_errors),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("binary", "base", "guidance_work", "out"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    a = parser.parse_args()
    work = a.guidance_work
    results = work / "results" / "modes.csv"
    if not results.is_file():
        results = work / "results-ca" / "modes.csv"
    with results.open() as source:
        labels = {
            (
                r["rec"],
                r["ligand"],
                r["state"],
                r["seed"],
                r["slot"],
                int(r["model"]),
            ): (int(r["hand"]) == 1 and float(r["rmsd"]) <= 2.0)
            for r in csv.DictReader(source)
            if r["panel"] == "ca" and r["search"] == "baseline"
        }
    a.out.mkdir(parents=True, exist_ok=True)
    tasks = []
    for case in sources(a.base):
        panel, _, lig, *_ = case
        if panel != "ca":
            continue
        for (state, _, _, search), ref, args, _ in pools(case, work):
            if search == "baseline":
                sites = work / "inputs" / f"ca_{lig}_{state}.amide"
                dest = a.out / ref.name
                if dest.resolve() == ref.resolve():
                    parser.error("Output must be separate from the baseline pool")
                tasks.append((ref, dest, args, sites))
    with ThreadPoolExecutor(a.jobs) as executor:
        rows = list(
            executor.map(lambda t: check_one(a.binary.resolve(), *t, labels), tasks)
        )
    if len(rows) != 80 or sum(r["models"] for r in rows) != len(labels):
        raise ValueError("Incomplete hCA II production conformance panel")
    with (a.out / "production.csv").open("w", newline="") as sink:
        writer = csv.DictWriter(sink, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for scope, group in (("all", rows), ("cross", [r for r in rows if r["cross"]])):
        print(
            scope,
            "pools",
            len(group),
            "native/amide/oracle joint hits",
            *(sum(r[k] for r in group) for k in ("native_hit", "ranked_hit", "oracle")),
        )
    print("Strict coordinate and native energy-component matches:", len(labels))
    print(
        "Pool winners / complete orders matching Python on rounded coordinates:",
        *(
            sum(r[k] for r in rows)
            for k in ("top_matches_python", "order_matches_python")
        ),
    )
    print(
        "Largest absolute penalty and ranking score errors (kcal/mol):",
        *(max(r[k] for r in rows) for k in ("max_penalty_error", "max_score_error")),
    )


if __name__ == "__main__":
    main()
