"""Exercise real cross-compound contact ranking through the production CLI."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from evaluate import GEOM, HERE, atoms, distinct_template, input_data, native_components, profile, shape, tanimoto


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def coordinates(block):
    return "\n".join(l for l in block.splitlines() if l.startswith(("ATOM  ", "HETATM")))


def remark(block, name):
    prefix = "REMARK VINA " + name + ":"
    return next(l[len(prefix):].split() for l in block.splitlines() if l.startswith(prefix))


def check(binary, baseline, directory, q, r, seed, cap, selection, refinement, panel):
    source = "brdt" if panel == "brdt" else "old"
    _, ligand, _ = input_data(source, q)
    meta, _, receptor = input_data(source, r)
    rec, donor = distinct_template(source, q, r)
    donor_profile = profile(rec, donor)
    (directory / "ligand").write_text(ligand)
    (directory / "receptor").write_text(receptor)
    (directory / "donor").write_text("".join(f"{e} {x} {y} {z}\n" for e, (x,y,z) in zip(donor[1], donor[0])))

    def dock(engine, stage=None, window=100, modes=9):
        cmd = [str(engine), "--ligand", "ligand", "--receptor", "receptor", "--out", "poses",
               "--cpu", "1", "--seed", str(seed), "--max_evals", str(cap), "--exhaustiveness", "3",
               "--num_modes", str(modes), "--energy_range", str(window), "--search_counters"]
        if selection == "late":
            cmd += ["--pose_selection", "late"]
        if refinement != "native":
            cmd += ["--refinement_optimizer", refinement]
        for axis, value in zip("xyz", meta["center"]):
            cmd += ["--center_" + axis, str(value), "--size_" + axis, "24"]
        if stage:
            cmd += ["--pose_template", "donor", "--template_stage", stage]
        run = subprocess.run(cmd, cwd=directory, capture_output=True, text=True, check=True, timeout=90)
        blocks = [s for s in (directory / "poses").read_text().split("ENDMDL") if s.strip()]
        accounting = [l for l in run.stdout.splitlines() if l.startswith(("Search trajectory:", "Search counters:"))]
        assert len(accounting) == 4
        if stage:
            counts = run.stdout.split("Template pose scoring: ", 1)[1].splitlines()[0]
            prefix = f"stage={stage} selection_evals=0 refinement_evals=0 rank_evals="
            assert counts.startswith(prefix) and int(counts[len(prefix):]) >= len(blocks)
            if window == 100:
                assert int(counts[len(prefix):]) == len(blocks)
        return blocks, accounting

    native, counters = dock(binary)
    if baseline:
        assert (native, counters) == dock(baseline)
    native_coords = [coordinates(p) for p in native]
    native_scores = [native_components(p) for p in native]
    data = {"query": q, "receiver": r, "seed": seed, "cap": cap, "selection": selection,
            "refinement": refinement, "native_pool_sha256": digest("ENDMDL".join(native) + "ENDMDL\n"),
            "search_accounting": counters, "modes": len(native)}
    case = (q, r, str(seed))
    if panel == "late_fixed":
        with (GEOM / "upstream/reciprocal/cells.csv").open() as stream:
            prior = next(row for row in csv.DictReader(stream)
                         if (row["query"], row["receiver"], row["seed"], row["arm"]) == (*case, "late"))
    else:
        with (HERE / "cells.csv").open() as stream:
            prior = next(row for row in csv.DictReader(stream)
                         if (row["query"], row["receiver"], row["seed"], row["stage"]) == (*case, panel))
    data["matches_archived_native_pool_sha256"] = data["native_pool_sha256"] == prior["pool_sha256"]
    data["native_coordinate_sha256"] = digest("\n".join(native_coords))
    stages = ("rank", "contact_rank") if (selection, refinement) == ("late", "native") else ("contact_rank",)
    for stage in stages:
        blocks, accounting = dock(binary, stage)
        if stage == "rank" and baseline:
            assert (blocks, accounting) == dock(baseline, stage)
        assert counters == accounting
        order, overlaps, diffs, predicted_order = [], [], [], []
        for block in blocks:
            actual_stage, = remark(block, "TEMPLATE STAGE")
            score = float(remark(block, "TEMPLATE RESULT")[0])
            native_score, rank = remark(block, "NATIVE")
            rank = int(rank)
            native_score = float(native_score)
            overlap = float(remark(block, "TEMPLATE OVERLAP")[0])
            parsed = atoms(block)
            if stage == "rank":
                oracle = shape(parsed, donor)
            else:
                both = [tanimoto(a,b) for a,b in zip(profile(rec, parsed), donor_profile)]
                active = [n for n in both if n is not None]
                oracle = sum(active)/len(active) if active else 0
            rounded_score, divisor = native_scores[rank-1]
            assert actual_stage == stage and coordinates(block) == native_coords[rank-1]
            assert abs(rounded_score-native_score) <= 0.000501
            assert abs(overlap-oracle) < 0.0004  # 0.001 A PDBQT output precision
            assert abs(score-(native_score-6*overlap/divisor)) < 1e-12
            assert abs(score-(native_score-6*oracle/divisor)) < 0.0024
            order.append(rank)
            predicted_order.append((native_score-6*oracle/divisor, rank))
            overlaps.append(abs(overlap-oracle))
            diffs.append(abs(score-(native_score-6*oracle/divisor)))
        assert sorted(order) == list(range(1,len(native)+1))
        assert order == [rank for _, rank in sorted(predicted_order)]
        pool_digest = digest("\n".join(text for _, text in sorted(zip(order, map(coordinates, blocks)))))
        assert pool_digest == data["native_coordinate_sha256"]
        data[stage] = dict(native_ranks=order, native_order_coordinate_sha256=pool_digest,
                           max_overlap_delta=max(overlaps), max_rank_score_delta=max(diffs))
    if panel != "late_fixed" and data["matches_archived_native_pool_sha256"]:
        assert data["contact_rank"]["native_ranks"][0] == int(prior["contact_rank"])
        data["archived_native_top1_rmsd_A"] = float(prior["native"])
        data["archived_contact_top1_rmsd_A"] = float(prior["contact"])
    if selection == "native":
        small, small_counts = dock(binary, modes=2)
        ranked, counts = dock(binary, "contact_rank", modes=2)
        assert len(small) == len(ranked) == 2 and small_counts == counts
        assert sorted(map(coordinates, small)) == sorted(map(coordinates, ranked))
        data["two_mode_cap_coordinate_parity"] = True
    # Place a nontrivial energy window between native affinities, independent
    # of either template. Reranking must retain exactly the eligible coordinates.
    index = len(native)//2
    window = (native_scores[index-1][0] + native_scores[index][0])/2 - native_scores[0][0]
    control, _ = dock(binary, window=window)
    guided, _ = dock(binary, "contact_rank", window)
    assert 0 < len(control) < len(native)
    assert sorted(map(coordinates, control)) == sorted(map(coordinates, guided))
    data["window_kcal_per_mol"], data["window_modes"] = window, len(control)
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary", type=Path)
    p.add_argument("--baseline", type=Path, help="build predating contact_rank for byte parity on Late and rank")
    args = p.parse_args()
    # The three replayed successes below are named regression cases for the
    # existing frozen-panel analysis. They do not provide independent efficacy estimates.
    cases = [
        ("4GHI", "5TBM", 85231, 12000, "late", "native", "late_fixed"),
        ("2NNQ", "3FR2", 85231, 12000, "late", "native", "late_fixed"),
        ("1FKG", "1FKH", 93563, 70000, "native", "native", "per_complex"),
        ("1FKH", "1FKG", 93529, 66000, "native", "native", "per_complex"),
        ("7MRG", "7MRH", 95701, 120000, "late", "geometry", "brdt"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        results = [check(args.binary.resolve(), args.baseline.resolve() if args.baseline else None, Path(tmp), *case)
                   for case in cases]
    print(json.dumps(results, indent=2))
