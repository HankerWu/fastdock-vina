"""Rebuild all independent scores and check the switched pairs from the PDBQT archives."""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from replay import KEY, RT, TOL, lower_bound, recover_factors


def load_archived_model(repo):
    path = repo / "benchmarks/mobile_water"
    sys.path.insert(0, str(path))
    from mobile_water import Site
    spec = importlib.util.spec_from_file_location("mobile_water_replay", path / "replay.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, Site


def main(repo, inputs, pairs_path, output):
    archive, site_class = load_archived_model(repo)
    results = repo / "benchmarks/soft_water_interaction/results"
    poses = pd.read_csv(results / "poses.csv")
    source = pd.read_csv(results / "changed_winner_sites.csv")
    pairs = pd.read_csv(pairs_path).set_index(KEY)
    n_pools = n_poses = n_changed = 0
    max_error = max_site_error = 0.0
    robust = []
    for panel, section in poses.groupby("panel"):
        work = inputs / panel
        meta = json.loads((work / "metadata.json").read_text())
        for receiver, rec in section.groupby("receiver"):
            lines = archive.pdb_lines(work / "raw" / f"{receiver}.pdb")
            protein = cKDTree([archive.xyz(l) for l in lines
                               if l.startswith("ATOM") and l[21] == "A"])
            sites = {name: site_class.prepare(x, b, protein)
                     for name, x, b in archive.cross_sites(lines, protein,
                                                          np.array(meta[receiver]["center"]))}
            for key, pool in rec.groupby(KEY):
                panel, _, cap, query, receiver, seed, search = key
                path = work / f"runs_{cap}/{query}_{receiver}_{seed}_{search}.pdbqt"
                ref = pool.set_index("model")
                per_site, dry = {}, {}
                for model, score, tree in archive.models(path):
                    per_site[model] = np.array([s.evaluate(tree)[1] for s in sites.values()])
                    dry[model] = score
                    max_error = max(max_error, abs(score - ref.loc[model].dry),
                                    abs(per_site[model].sum() - ref.loc[model].water_penalty))
                if set(dry) != set(ref.index):
                    raise ValueError("Coordinate and archived pools differ")
                picked = min(dry, key=lambda i: (dry[i] + per_site[i].sum(), i))
                if picked != ref.hydrated.idxmin():
                    raise ValueError("Coordinate winner disagrees with archive")
                n_poses += len(dry)
                n_pools += 1
                if picked == 1:
                    continue
                n_changed += 1
                tab = source
                for col, val in zip(KEY, key, strict=True):
                    tab = tab[tab[col] == val]
                u = np.stack([per_site[1], per_site[picked]])
                expected = tab.set_index("site").loc[list(sites), [
                    "penalty_original", "penalty_winner"]].to_numpy().T
                max_site_error = max(max_site_error, np.abs(expected - u).max())
                a = recover_factors(u)
                record = pairs.loc[key]
                # Each printed dry score was rounded to 0.001, so the true gap can
                # differ by up to 0.001. For either claim test the adverse endpoint.
                gap = dry[picked] - dry[1]
                if record.result == "forced":
                    decision = lower_bound(a, gap + .001)["lower"] > TOL
                elif record.result == "reversible":
                    decision = lower_bound(a, gap - .001)["feasible_margin"] < -TOL
                else:
                    decision = False
                robust.append(bool(decision))
    if n_changed != len(pairs) or max(max_error, max_site_error) > 1e-10:
        raise ArithmeticError("Archived scores and regenerated coordinates disagree")
    report = dict(pools=n_pools, poses=n_poses, changed=n_changed,
                  robust_to_dry_rounding=sum(robust), max_score_component_error=max_error,
                  max_single_site_error=max_site_error)
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("inputs", type=Path)
    p.add_argument("pairs", type=Path)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    main(args.repo, args.inputs, args.pairs, args.output)
