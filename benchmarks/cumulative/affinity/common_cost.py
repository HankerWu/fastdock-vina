"""Use measured CPU alone to target a reachable budget across actual engines."""

import argparse
import json
import math
import os
from pathlib import Path
from run import ENGINES, cells, digest, run_cell


def interpolate(anchors, target):
    """Select the closest lower/upper observations, handling noisy ties in CPU."""
    lower = [(cpu, cap) for cap, cpu in anchors if cpu <= target]
    upper = [(cpu, cap) for cap, cpu in anchors if cpu >= target]
    if not lower or not upper:
        return min(anchors, key=lambda pair: abs(pair[1] - target))[0]
    low, lo_cap = max(lower)
    high, hi_cap = min(upper)
    if low == high:
        return lo_cap
    estimate = lo_cap + (target - low) / (high - low) * (hi_cap - lo_cap)
    return min(384000, max(3000, round(estimate / 100) * 100))


def measure(name, caps, trials, binaries, identities, cycle, attempts):
    if not caps:
        return {}
    for engine in caps:
        (Path("outputs") / engine / name).mkdir(parents=True, exist_ok=True)
    times = {engine: 0.0 for engine in caps}
    for idx, job in enumerate(trials):
        offset = (idx + cycle) % len(ENGINES)
        for engine in ENGINES[offset:] + ENGINES[:offset]:
            if engine not in caps:
                continue
            run_cell(job, engine, binaries[engine], name, caps[engine], identities[engine])
            state, row, seed = job
            path = Path("outputs") / engine / name / f"{state}_{row['compound']}_{seed}.json"
            record = json.loads(path.read_text())
            if not record["valid"]:
                raise RuntimeError(f"Failed calibration cell {path}")
            times[engine] += record["user_s"] + record["sys_s"]
    for e in caps:
        attempts.append(dict(engine=e, stage=name, cap=caps[e], cpu_s=times[e]))
    return times


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument(
        "--replace-candidate",
        type=Path,
        metavar="PREVIOUS_TIMING_JSON",
        help="keep existing baseline timings and use their registered target for the successor",
    )
    for e in ENGINES:
        parser.add_argument("--" + e, type=Path)
    args = parser.parse_args()
    previous = json.loads(args.replace_candidate.read_text()) if args.replace_candidate else None
    engines = ("candidate",) if previous else ENGINES
    bins = {e: getattr(args, e).resolve(strict=True) for e in engines}
    identities = {e: digest(p) for e, p in bins.items()}
    os.chdir(args.work)
    cpu = max(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    trials = list(cells("low"))
    attempts = [r for r in previous["attempts"] if r["engine"] != "candidate"] if previous else []

    def batch(name, caps, cycle):
        return measure(name, caps, trials, bins, identities, cycle, attempts)

    low = batch("cpu3k", dict.fromkeys(engines, 3000), 0)
    high = batch("cpu96k", dict.fromkeys(engines, 96000), 1)
    target = previous["target_s"] if previous else math.ceil(1.2 * max(low.values()) * 10) / 10
    batch("cpu384k", {e: 384000 for e in engines if high[e] < target}, 2)

    def cap(e):
        return interpolate([(a["cap"], a["cpu_s"]) for a in attempts if a["engine"] == e], target)

    confirmation = batch("confirmation1", {e: cap(e) for e in engines}, 3)
    retry = {e: cap(e) for e in engines if abs(confirmation[e] / target - 1) > 0.15}
    batch("confirmation2", retry, 4)
    final = dict(previous["final"]) if previous else {}
    final.update({e: ("confirmation2" if e in retry else "confirmation1") for e in engines})
    Path("common_cost.json").write_text(
        json.dumps(dict(logical_cpu=cpu, target_s=target, attempts=attempts, final=final), indent=2)
        + "\n"
    )


if __name__ == "__main__":
    main()
