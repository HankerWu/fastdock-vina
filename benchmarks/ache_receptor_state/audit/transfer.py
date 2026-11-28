"""Independently reconstruct the two-source selection and all transfer coordinates."""

import argparse
import json
import re
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from audit import STUDY, Frame, all_coordinates, compare, counters, rmsd, score_log, sha

SEEDS = (1103501, 1103503, 1103507)
RECEIVERS = ("4EY6", "4EY7")


def audit_transfer(out, archive):
    wall, cpu = time.perf_counter(), time.process_time()
    blind = STUDY / "blind/results/raw"
    out.mkdir(parents=True, exist_ok=True)
    jobs, endpoints, source_records = {}, {}, {}
    for rec in RECEIVERS:
        for seed in SEEDS:
            pools = (blind / rec / f"pro24k_{seed}/dry.pdbqt").read_text().split("ENDMDL")[:-1]
            for rank, block in enumerate(pools, 1):
                source_records[rec, seed, rank] = block
    with tempfile.TemporaryDirectory() as tmp, tarfile.open(archive) as tar:
        frame = Frame(Path(tmp))
        members = {m.name: m for m in tar.getmembers() if m.isfile()}
        for rec in RECEIVERS:
            for (donor, seed, rank), orig in source_records.items():
                prefix = f"{rec}/{donor}_{seed}_{rank:02d}"

                def get(file, prefix=prefix):
                    return tar.extractfile(members[prefix + "/" + file]).read()

                x = get("start.pdbqt").decode()
                y = get("end.pdbqt").decode()
                assert get("end.pdbqt") == get("observer.pdbqt")
                obs = re.findall(r"REFINE_OBSERVE\s+(\S+)\s+(\S+)\s+(\d+)", get("observer.log").decode())
                assert obs
                err = np.max(np.abs(all_coordinates(x) - frame.moves[donor, rec](all_coordinates(orig))))
                if err > 0.000501:
                    raise ValueError(f"Coordinate transfer mismatch {prefix}: {err}")
                xyz = frame.heavy(y)
                starts = frame.heavy(x)
                key = rec, donor, seed, rank
                jobs[key] = {
                    "receiver": rec,
                    "source": donor,
                    "source_seed": seed,
                    "source_mode": rank,
                    "score": score_log(get("score.log").decode()),
                    "start_A": rmsd(starts, frame.targets[rec]),
                    "end_A": rmsd(xyz, frame.targets[rec]),
                    "calls": sum(int(o[2]) for o in obs),
                    "start_sha256": sha(get("start.pdbqt")),
                    "end_sha256": sha(get("end.pdbqt")),
                    "source_rmsd_A": rmsd(frame.heavy(orig), frame.targets[donor]),
                    "max_transform_component_A": float(err),
                }
                endpoints[key] = xyz
    outcomes, cells = [], []

    def search_calls(host, method, seed):
        counts = counters((blind / host / f"{method}_{seed}/dry.log").read_text())
        return sum(v for k, v in counts.items() if k.endswith("evals"))

    for rec in RECEIVERS:
        other = next(r for r in RECEIVERS if r != rec)
        for si, seed in enumerate(SEEDS):
            nextseed = SEEDS[(si + 1) % len(SEEDS)]
            arms = {
                "recipient": [(rec, seed)],
                "cross": [(rec, seed), (other, seed)],
                "same": [(rec, seed), (rec, nextseed)],
            }
            for arm, sources in arms.items():
                candidates = []
                for donor, source_seed in sources:
                    candidates.extend(k for k in jobs if k[:3] == (rec, donor, source_seed))
                # Python's stable sort retains recipient-first and input rank on ties.
                ranked = sorted(candidates, key=lambda k: jobs[k]["score"])
                keep = []
                for k in ranked:
                    row = jobs[k]
                    duplicate = any(rmsd(endpoints[k], endpoints[j]) < 1.0 for j in keep)
                    if row["score"] > jobs[ranked[0]]["score"] + 100.0:
                        reason = "window"
                    elif duplicate:
                        reason = "duplicated"
                    elif len(keep) >= 9:
                        reason = "rank_cut"
                    else:
                        reason = "retained"
                        keep.append(k)
                    outcomes.append(
                        {
                            "receiver": rec,
                            "seed": seed,
                            "arm": arm,
                            "source": k[1],
                            "source_seed": k[2],
                            "source_mode": k[3],
                            "rank": len(keep) if reason == "retained" else np.nan,
                            "decision": reason,
                        }
                    )
                own = [endpoints[k] for k in ranked if k[1:3] == (rec, seed)]
                distinct = [
                    k
                    for k in ranked
                    if k[1:3] != (rec, seed)
                    and jobs[k]["end_A"] < 2.0
                    and all(rmsd(endpoints[k], coord) >= 1.0 for coord in own)
                ]
                ds = [jobs[k]["end_A"] for k in keep]
                cells.append(
                    {
                        "receiver": rec,
                        "seed": seed,
                        "arm": arm,
                        "top1_A": ds[0],
                        "best5_A": min(ds[:5]),
                        "oracle_A": min(ds),
                        "regret_A": ds[0] - min(ds),
                        "retained": len(ds),
                        "distinct_donor_sub2": len(distinct),
                        "distinct_donor_sub2_returned": sum(k in keep for k in distinct),
                        "distinct_donor_sub2_top1": int(keep[0] in distinct),
                        "published": len(ds),
                        "search_gradient_evals": sum(search_calls(a, "pro24k", b) for a, b in sources),
                        "local_gradient_evals": sum(jobs[k]["calls"] for k in ranked),
                        "local_best_A": min(jobs[k]["end_A"] for k in ranked),
                        "origin_best_A": min(jobs[k]["source_rmsd_A"] for k in ranked),
                        "start_best_A": min(jobs[k]["start_A"] for k in ranked),
                    }
                )
            pool = (blind / rec / f"pro48k_{seed}/dry.pdbqt").read_text().split("ENDMDL")[:-1]
            ds = [rmsd(frame.heavy(m), frame.targets[rec]) for m in pool]
            for i in range(len(ds)):
                outcomes.append(
                    {
                        "receiver": rec,
                        "seed": seed,
                        "arm": "double",
                        "source": rec,
                        "source_seed": seed,
                        "source_mode": i + 1,
                        "rank": i + 1,
                        "decision": "retained",
                    }
                )
            cells.append(
                {
                    "receiver": rec,
                    "seed": seed,
                    "arm": "double",
                    "top1_A": ds[0],
                    "best5_A": min(ds[:5]),
                    "oracle_A": min(ds),
                    "regret_A": ds[0] - min(ds),
                    "retained": len(ds),
                    "published": len(ds),
                    "distinct_donor_sub2": 0,
                    "distinct_donor_sub2_returned": 0,
                    "distinct_donor_sub2_top1": 0,
                    "search_gradient_evals": search_calls(rec, "pro48k", seed),
                    "local_gradient_evals": 0,
                    "local_best_A": min(ds),
                    "origin_best_A": min(ds),
                    "start_best_A": min(ds),
                }
            )
    audited = pd.DataFrame(outcomes)
    results = pd.DataFrame(cells)
    return (
        jobs,
        audited,
        results,
        {
            "archive_sha256": sha(archive.read_bytes()),
            "distinct_local_jobs": len(jobs),
            "observed_gradient_calls": sum(v["calls"] for v in jobs.values()),
            "max_transform_component_error_A": max(
                v["max_transform_component_A"] for v in jobs.values()
            ),
            "computed_cells": len(results),
            "reconstruction_wall_s": time.perf_counter() - wall,
            "reconstruction_cpu_s": time.process_time() - cpu,
        },
    )


def main(path, out):
    wall, cpu = time.perf_counter(), time.process_time()
    archive = path / "raw.tar.xz"
    jobs, selections, cells, receipt = audit_transfer(out, archive)
    primary_jobs = pd.read_csv(path / "jobs.csv")
    for r in primary_jobs.itertuples():
        row = jobs[r.receiver, r.source, r.seed, r.mode]
        assert row["score"] == r.score and row["calls"] == r.gradient_calls
        assert row["start_sha256"] == r.start_sha256 and row["end_sha256"] == r.end_sha256
    assert len(jobs) == len(primary_jobs)
    # Comparison to primary implementation is performed after independent filtering and measurement.
    primary = pd.read_csv(path / "selections.csv")
    keys = ["receiver", "seed", "arm", "source", "source_seed", "source_mode"]
    published = primary[primary.arm.isin(selections.arm.unique())].set_index(keys).sort_index()
    rerun = selections.set_index(keys).sort_index()
    if not published.index.is_unique or not published.index.equals(rerun.index):
        raise ValueError("Source membership differs")
    kept = rerun.decision == "retained"
    assert (kept == (published.decision == "retained")).all()
    assert np.array_equal(rerun.loc[kept, "rank"], published.loc[kept, "rank"])
    # Inspect export files consumed by downstream callers, independently of the primary manifest.
    export_modes = 0
    for (receiver, seed, arm), group in selections.groupby(["receiver", "seed", "arm"]):
        path_out = path / "published" / f"{receiver}_{arm}_{seed}.pdbqt"
        winners = group[group.decision == "retained"].sort_values("rank")
        if arm == "double":
            native = STUDY / "blind/results/raw" / receiver / f"pro48k_{seed}" / "dry.pdbqt"
            assert path_out.read_bytes() == native.read_bytes()
        blocks = path_out.read_text().split("ENDMDL")[:-1]
        assert len(blocks) == len(winners)
        for block, item in zip(blocks, winners.itertuples(), strict=True):
            if arm != "double":
                text = block.strip().splitlines()
                assert text[0] == f"MODEL {int(item.rank)}"
                assert (
                    float(text[1].split()[3])
                    == jobs[receiver, item.source, item.source_seed, item.source_mode]["score"]
                )
                coordinates = ("\n".join(text[3:]) + "\n").encode()
                assert (
                    sha(coordinates)
                    == jobs[receiver, item.source, item.source_seed, item.source_mode]["end_sha256"]
                )
            export_modes += 1
    receipt["verified_export_modes"] = export_modes
    receipt["primary_membership_rank_checked"] = int(kept.sum())
    receipt["candidate_memberships_checked"] = len(rerun)
    separated = pd.concat(
        [published.decision.rename("published"), rerun.decision.rename("independent")], axis=1
    )
    receipt["nonretained_label_differences"] = int(
        (separated["published"] != separated["independent"]).sum()
    )
    assert receipt["nonretained_label_differences"] == 0
    k = ["receiver", "seed", "arm"]
    receipt["primary_metric_max_error"] = compare(
        cells,
        path / "cells.csv",
        k,
        [
            "top1_A",
            "best5_A",
            "oracle_A",
            "regret_A",
            "published",
            "origin_best_A",
            "start_best_A",
            "local_best_A",
            "distinct_donor_sub2",
            "distinct_donor_sub2_returned",
            "distinct_donor_sub2_top1",
        ],
    )
    receipt["charged_gradient_evaluation_max_error"] = compare(
        cells, path / "costs.csv", k, ["search_gradient_evals", "local_gradient_evals"]
    )
    cells.drop(columns=["retained"]).to_csv(out / "transfer_cells.csv", index=False)
    receipt.update(audit_wall_s=time.perf_counter() - wall, audit_cpu_s=time.process_time() - cpu)
    (out / "transfer_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("primary", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    main(a.primary, a.output)
