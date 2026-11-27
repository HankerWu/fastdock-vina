"""Replay archived coordinates; compute bounds without optimizing occupancies from labels."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import time

import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree

from envelope import assess, check_enumeration, RT

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH / "mobile_water"))
from mobile_water import Site, accessible, penalty
from replay import cross_sites, models, pdb_lines, xyz

KEY = ["panel", "target", "cap", "query", "receiver", "seed", "search"]
ARCHIVES = {"historical": ["paired_crossdocking/fixtures.tar.xz",
                            "paired_crossdocking/poses-3000.tar.xz",
                            "paired_crossdocking/poses-12000.tar.xz"],
            "transfer": ["torsion_proposals/fixtures/inputs.tar.xz",
                         "soft_water_interaction/transfer_outputs.tar.xz"]}


def availability(site, ligand):
    if ligand.query(site.center)[0] >= site.radius + 2.6:
        return 1.
    return float(np.clip(site.weights @ accessible(ligand.query(site.samples)[0]), 0., 1.))


def replay(work, order, poses, old_cells, outputs):
    modes, cells, factors = outputs
    max_error = 0.
    for (panel, receiver), group in poses.groupby(["panel", "receiver"]):
        base = work / panel
        meta = json.loads((base / "metadata.json").read_text())
        lines = pdb_lines(base / "raw" / f"{receiver}.pdb")
        protein = cKDTree([xyz(line) for line in lines if line.startswith("ATOM") and line[21] == "A"])
        records = cross_sites(lines, protein, np.asarray(meta[receiver]["center"]))
        sites = [Site.prepare(position, b, protein, order=order) for _, position, b in records]
        for key, table in group.groupby(KEY):
            table = table.sort_values("model")
            _, _, cap, query, _, seed, search = key
            pool = base / f"runs_{cap}/{query}_{receiver}_{seed}_{search}.pdbqt"
            coords = list(models(pool))
            ids = np.array([mode for mode, _, _ in coords])
            energy = np.array([energy for _, energy, _ in coords])
            if not (np.array_equal(ids, table.model) and np.array_equal(energy, table.dry)
                    and (len(sites) == table.n_sites).all()):
                raise ValueError(f"Coordinates, scores or sites differ at {key}")
            a = np.array([[availability(site, tree) for site in sites] for _, _, tree in coords])
            rmsd = table.rmsd.to_numpy()
            hydrated = energy + penalty(a).sum(axis=1)
            im = int(np.argmin(hydrated))
            if order == 11:
                max_error = max(max_error, np.max(np.abs(hydrated - table.hydrated)))
                if ids[im] != old_cells.loc[(*key, "water"), "model"] or max_error > 1e-10:
                    raise ValueError(f"Historical score mismatch {key}")
            details, lower, upper, mask, winners = assess(energy, a, rmsd)
            io = int(np.argmin(rmsd))
            prefix = dict(order=order, **dict(zip(KEY, key, strict=True)))
            for i, (rival, margin) in enumerate(details):
                modes.append(dict(**prefix, model=ids[i], rmsd=rmsd[i], dry=energy[i],
                                  independent=hydrated[i], dominated=rival is not None,
                                  rival=ids[rival] if rival is not None else None,
                                  strongest_pair_margin=margin))
                for (label, _, _), value in zip(records, a[i], strict=True):
                    factors.append(dict(**prefix, model=ids[i], site=label, availability=value))
            dry = rmsd[0]
            cells.append(dict(**prefix, modes=len(ids), sites=len(sites),
                              dry_top=dry, water_top=rmsd[im], water_model=ids[im],
                              best5_dry=min(rmsd[:5]), oracle=rmsd[io], oracle_id=ids[io],
                              dry_regret=dry-rmsd[io], oracle_dominated=details[io][0] is not None,
                              eliminated=sum(j is not None for j, _ in details),
                              best_possible_lower=lower, witness_upper=upper,
                              lower_gap=lower-rmsd[io], closed=abs(upper-lower) < 1e-10,
                              occupied_sites="|".join(record[0] for record, use in zip(records, mask) if use),
                              plausible_winners="|".join(str(x) for x in ids[winners]),
                              pool_sha256=hashlib.sha256(pool.read_bytes()).hexdigest()))
    return max_error


def main(out):
    wall, cpu = time.perf_counter(), time.process_time()
    source = BENCH / "soft_water_interaction/results"
    poses = pd.read_csv(source / "poses.csv", float_precision="round_trip")
    old_cells = pd.read_csv(source / "cells.csv").set_index(KEY + ["ranker"])
    modes, cells, factors = [], [], []
    files = [source / "poses.csv", source / "cells.csv"]
    checks = check_enumeration()
    costs = {}
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for panel, archives in ARCHIVES.items():
            for relative in archives:
                archive = BENCH / relative
                files.append(archive)
                with tarfile.open(archive) as tar:
                    tar.extractall(work / panel, filter="data")
        for order in (11, 17):
            w, c = time.perf_counter(), time.process_time()
            error = replay(work, order, poses, old_cells, (modes, cells, factors))
            costs[order] = dict(wall_s=time.perf_counter()-w, cpu_s=time.process_time()-c,
                                max_historical_score_error=error if order == 11 else None)
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in (("modes", modes), ("factors", factors), ("cells", cells)):
        pd.DataFrame(rows).to_csv(out / (name + ".csv.gz"), index=False,
                                  compression={"method": "gzip", "mtime": 0})
    stats = []
    df = pd.DataFrame(cells)
    for cols in (["order"], ["order", "panel"], ["order", "target"]):
        for key, block in df.groupby(cols):
            for col in ["dry_top", "water_top", "best5_dry", "oracle", "dry_regret",
                        "best_possible_lower", "witness_upper", "lower_gap"]:
                x = block[col]
                quartiles = x.quantile([.25, .5, .75, .9, .95]).to_numpy()
                stats.append(dict(**dict(zip(cols, key)), metric=col, n=len(x),
                                  mean=x.mean(), **dict(zip(("q1","median","q3","p90","p95"), quartiles)),
                                  maximum=x.max()))
    pd.DataFrame(stats).to_csv(out / "distributions.csv", index=False)
    info = dict(enumerated_pair_checks=checks, orders=costs, numpy=np.__version__,
                scipy=scipy.__version__, pandas=pd.__version__, python=sys.version.split()[0],
                pool_count=len(cells)//2, mode_count=len(modes)//2,
                input_sha256={str(p.relative_to(BENCH)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                wall_s=time.perf_counter()-wall, cpu_s=time.process_time()-cpu)
    (out / "manifest.json").write_text(json.dumps(info, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
