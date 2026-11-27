"""Enumerate every ordering attainable with one ligand-independent two-state weight."""
import argparse
import csv
import itertools
import math
from pathlib import Path

import numpy as np

RT = 0.00198720425864 * 298.15
STATES = ("2J4Z", "2X6E")


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    with path.open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)


def signed_log_difference(near, blind):
    """Sign and log magnitude of exp(-near/RT) minus exp(-blind/RT)."""
    if near == blind:
        return 0, -math.inf
    sign = 1 if near < blind else -1
    log_abs = -min(near, blind) / RT + math.log(-math.expm1(-abs(near - blind) / RT))
    return sign, log_abs


def summarize(path):
    src = read(path / "state_scores.csv")
    energies = {(r["origin"], r["ligand"], r["pose"], r["state"], r["registration"],
                 int(r["nodes_per_scale"])): r for r in src}
    if len(energies) != len(src):
        raise ValueError("Duplicate state row")
    centers = read(path / "centers.csv")
    pairs, bounds, geometry = [], [], []
    frames = sorted({r["registration"] for r in src})
    counts = sorted({int(r["nodes_per_scale"]) for r in src})
    for frame, n, field in itertools.product(frames, counts, ("center", "window")):
        preferences = []
        rmsds = []
        for row in centers:
            origin, drug, pose = row["origin"], row["ligand"], row["pose"]
            if pose == "witness":
                continue
            near = []
            blind = []
            for state in STATES:
                key = (state, "pocket" if state == origin else frame, n)
                near.append(float(energies[origin, drug, "witness", *key][field]))
                blind.append(float(energies[origin, drug, pose, *key][field]))
            A, B = (signed_log_difference(x, y) for x, y in zip(near, blind, strict=True))
            cross = A[1] - B[1] if A[0] * B[0] < 0 else None
            preferences.append((A[0], B[0], cross))
            def half(v):
                return -RT * (np.logaddexp(-v[0]/RT, -v[1]/RT) - math.log(2))
            margin = half(blind) - half(near)
            near_rms = next(float(r["rmsd_A"]) for r in centers if
                            (r["origin"], r["ligand"], r["pose"]) == (origin, drug, "witness"))
            blind_rms = float(row["rmsd_A"])
            if near_rms >= blind_rms:
                raise ValueError("Witness label inconsistent with experimental geometry")
            rmsds.append((blind_rms, near_rms, near_rms if margin > 0 else blind_rms))
            pairs.append(dict(registration=frame, nodes=n, model=field, origin=origin, ligand=drug,
                              seed=pose, near_2J4Z=near[0], blind_2J4Z=blind[0],
                              near_2X6E=near[1], blind_2X6E=blind[1], margin_half=margin,
                              sign_2J4Z=A[0], sign_2X6E=B[0], crossover_log_odds_2X6E=cross))
        roots = sorted({c for _, _, c in preferences if c is not None})
        # The margin is affine in p. A strict win at a crossing also holds
        # immediately on at least one side; a tie never adds a strict win.
        probes = [(-math.inf, -math.inf, -math.inf), (math.inf, math.inf, math.inf)]
        for lo, hi in itertools.pairwise([-math.inf, *roots, math.inf]):
            if lo == -math.inf and hi == math.inf:
                x = 0.
            elif lo == -math.inf:
                x = hi - 1.
            elif hi == math.inf:
                x = lo + 1.
            else:
                x = (lo + hi) / 2.
            probes.append((lo, hi, x))

        def wins(a, b, crossing, x):
            if x == -math.inf:
                return a > 0
            if x == math.inf:
                return b > 0
            if crossing is None:
                return max(a,b) > 0
            return (a > 0 and x < crossing) or (b > 0 and x > crossing)

        scores = [(sum(wins(*pref, x) for pref in preferences), lo, hi) for lo, hi, x in probes]
        best = max(k for k, _, _ in scores)
        bounds.append(dict(registration=frame, nodes=n, model=field, pairs=len(preferences),
                           wins_2J4Z=scores[0][0], wins_2X6E=scores[1][0],
                           wins_home=sum(p[f"sign_{p['origin']}"] > 0 for p in pairs[-len(preferences):]),
                           wins_half=sum(wins(*p, 0) for p in preferences),
                           max_shared=best, max_per_pair=sum(max(a,b)>0 for a,b,_ in preferences),
                           best_log_odds_intervals=" | ".join(f"({lo:.9g},{hi:.9g})"
                              for score,lo,hi in scores if score==best)))
        rr = np.asarray(rmsds)
        for label, arr in (("blind_center", rr[:,0]), ("oracle", rr[:,1]), ("equal_mix", rr[:,2]),
                           ("equal_minus_blind", rr[:,2]-rr[:,0]), ("equal_regret", rr[:,2]-rr[:,1])):
            q1,med,q3,p90,p95 = np.percentile(arr,[25,50,75,90,95])
            geometry.append(dict(registration=frame,nodes=n,model=field,metric=label,n=len(arr),
                                 mean=np.mean(arr),median=med,q25=q1,q75=q3,p90=p90,p95=p95,
                                 max=np.max(arr)))
    write(path / "pairs.csv", pairs)
    write(path / "bounds.csv", bounds)
    write(path / "conditional_geometry.csv", geometry)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    summarize(parser.parse_args().results)
