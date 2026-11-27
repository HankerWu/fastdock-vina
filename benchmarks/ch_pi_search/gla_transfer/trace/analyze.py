"""Attribute lost DGJ pose coverage at observed optimizer and output checkpoints."""

import argparse
import csv
import json
import lzma
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from verify import geometry

STAGES = (
    "hunt",
    "accepted",
    "second",
    "trajectory",
    "merged",
    "selected",
    "refined",
    "retained",
    "published",
)


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def signature(event):
    return tuple(event[3].ravel())


def stages_for(path, ref, chemical_order):
    lines = path.read_text().splitlines()
    (atom_line,) = [l for l in lines if l.startswith("atoms ")]
    serials = [int(n) for n in atom_line.split()[1:]]
    order = [serials.index(i) for i in chemical_order]
    groups = defaultdict(list)
    for l in lines:
        fields = l.split()
        if fields[0] == "atoms":
            continue
        stage, trajectory, index, energy = fields[:4]
        xyz = np.asarray(fields[4:], dtype=float).reshape(len(order), 3)[order]
        rmsd = float(np.sqrt(((xyz - ref) ** 2).sum(-1).mean(-1)).min())
        groups[stage].append((int(trajectory), int(index), float(energy), xyz, rmsd))
    return groups


def analyze(work, runs, out):
    out.mkdir(parents=True, exist_ok=True)
    modes = read_csv(HERE.parent / "results/modes.csv")
    costs = read_csv(HERE.parent / "results/costs.csv")
    cells, summaries, lineage, witnesses, censored = [], [], [], [], []
    event_rows = []
    for case in costs:
        tag = "_".join(
            case[k] for k in ("panel", "query", "receptor", "seed", "cap", "search")
        )
        ref, chemical_order = geometry(work, "GLA", case["query"], case["receptor"])
        events = stages_for(runs / (tag + ".trace"), ref, chemical_order)
        identity = {k: case[k] for k in ("receptor", "seed", "search")}
        log = (runs / (tag + ".log")).read_text()
        for tid, hunts, seconds in re.findall(
            r"Search trajectory: index=(\d+).*?hunt_calls=(\d+).*?second_calls=(\d+)",
            log,
        ):
            trajectory = int(tid)
            assert sum(e[0] == trajectory for e in events["hunt"]) == int(hunts)
            assert sum(e[0] == trajectory for e in events["second"]) == int(seconds)
        for tid, accepted in re.findall(
            r"Hunt chain: index=(\d+).*?accepted=(\d+)", log
        ):
            assert sum(e[0] == int(tid) for e in events["accepted"]) == int(accepted)

        # Current Late concatenates all trajectory returns, with no intervening filter.
        def count_xyz(stage):
            return Counter(signature(e) for e in events[stage])

        assert count_xyz("trajectory") == count_xyz("merged")
        assert not count_xyz("trajectory") - count_xyz("second")
        assert not count_xyz("selected") - count_xyz("merged")
        assert not count_xyz("retained") - count_xyz("refined")
        assert [x[:2] for x in events["selected"]] == [x[:2] for x in events["refined"]]
        for stage in ("hunt", "accepted", "second"):
            assert len({e[:2] for e in events[stage]}) == len(events[stage])
        by_id = {
            stage: {e[:2]: e for e in events[stage]}
            for stage in ("hunt", "accepted", "second")
        }
        for event_id, e in by_id["accepted"].items():
            assert signature(by_id["hunt"][event_id]) == signature(e)
        assert set(by_id["second"]) <= set(by_id["accepted"])
        rows = [r for r in modes if all(r[k] == v for k, v in identity.items())]
        assert len(rows) == len(events["retained"])
        max_diff = max(
            abs(float(r["rmsd"]) - e[4]) for r, e in zip(rows, events["retained"])
        )
        assert max_diff <= 0.001
        # Reconstruct the actual energy-ordered 1 A distinct-pose release in all 12 runs.
        by_energy = sorted(events["refined"], key=lambda e: e[2])
        chosen = []
        for e in by_energy:
            if len(chosen) == 9:
                break
            if all(
                np.sqrt(((e[3] - prev[3]) ** 2).sum(-1).mean()) >= 1.0
                for prev in chosen
            ):
                chosen.append(e)
        assert [signature(e) for e in chosen] == [
            signature(e) for e in events["retained"]
        ]
        for rank, e in enumerate(by_energy, 1):
            if e[4] > 2 or signature(e) in {signature(r) for r in chosen}:
                continue
            nearest = min(
                np.sqrt(((e[3] - prev[3]) ** 2).sum(-1).mean()) for prev in chosen
            )
            censored.append(
                dict(
                    **identity,
                    refinement_input=e[1],
                    rmsd_A=e[4],
                    native_raw_score=e[2],
                    rank_by_score=rank,
                    refinement_inputs=len(by_energy),
                    last_retained_score=chosen[-1][2],
                    nearest_retained_A=nearest,
                )
            )
        observed = {}
        for stage in STAGES:
            if stage == "published":
                distances = [float(r["rmsd"]) for r in rows]
            else:
                stage_events = events[stage]
                best = min(stage_events, key=lambda e: e[4])
                witnesses.append(
                    dict(
                        **identity,
                        stage=stage,
                        trajectory=best[0],
                        index=best[1],
                        score=best[2],
                        rmsd=best[4],
                        chemical_serial_order=chemical_order,
                        xyz=best[3].tolist(),
                    )
                )
                distances = [e[4] for e in stage_events]
                for traj, idx, energy, coords, rmsd in stage_events:
                    event_rows.append((tag, stage, traj, idx, energy, rmsd))
            observed[stage] = sum(d <= 2 for d in distances)
            summaries.append(
                dict(
                    **identity,
                    stage=stage,
                    count=len(distances),
                    best_A=min(distances),
                    within_2A=observed[stage],
                )
            )
        stages_hit = [s for s in STAGES if observed[s]]
        losses = [
            a + "->" + b
            for a, b in zip(STAGES, STAGES[1:])
            if observed[a] and not observed[b]
        ]
        cells.append(
            dict(
                **identity,
                first_observed_hit=stages_hit[0] if stages_hit else "none",
                last_observed_hit=stages_hit[-1] if stages_hit else "none",
                transitions_to_zero="|".join(losses) or "none",
                retained_print_max_delta_A=max_diff,
            )
        )
        h, a, s = (by_id[t] for t in ("hunt", "accepted", "second"))
        accepted_hits = [i for i, e in a.items() if e[4] <= 2]
        retained_by_traj = defaultdict(set)
        for e in events["trajectory"]:
            retained_by_traj[e[0]].add(signature(e))
        lineage.append(
            dict(
                **identity,
                near_hunt_rejected=sum(e[4] <= 2 and i not in a for i, e in h.items()),
                near_accepted_without_second=sum(i not in s for i in accepted_hits),
                near_accepted_second_moved_far=sum(
                    i in s and s[i][4] > 2 for i in accepted_hits
                ),
                near_accepted_second_stayed_near=sum(
                    i in s and s[i][4] <= 2 for i in accepted_hits
                ),
                far_accepted_second_became_near=sum(
                    e[4] <= 2 and a[i][4] > 2 for i, e in s.items()
                ),
                near_second_coord_absent_from_same_trajectory_return=sum(
                    e[4] <= 2 and signature(e) not in retained_by_traj[e[0]]
                    for e in s.values()
                ),
                near_selected_refined_far=sum(
                    before[4] <= 2 and after[4] > 2
                    for before, after in zip(events["selected"], events["refined"])
                ),
                far_selected_refined_near=sum(
                    before[4] > 2 and after[4] <= 2
                    for before, after in zip(events["selected"], events["refined"])
                ),
            )
        )
    write_table(out / "stages.csv", summaries)
    write_table(out / "cells.csv", cells)
    write_table(out / "lineage.csv", lineage)
    write_table(out / "lost_at_release.csv", censored)
    with lzma.open(out / "events.csv.xz", "wt") as stream:
        w = csv.writer(stream, lineterminator="\n")
        w.writerow(("run", "stage", "trajectory", "index", "stage_score", "rmsd_A"))
        w.writerows(event_rows)
    with lzma.open(out / "best_coordinates.json.xz", "wt") as stream:
        json.dump(witnesses, stream, separators=(",", ":"), allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "runs", "output"):
        p.add_argument(name, type=lambda s: Path(s).resolve())
    a = p.parse_args()
    analyze(a.work, a.runs, a.output)
