"""Compare coarse and explicit potentials before and after refinement on identical stored coordinates."""

import argparse
import csv
import itertools
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from run import parse_pose, save
from scan import read_ligand, receptor_maps


def main():
    p = argparse.ArgumentParser()
    for f in ("dev", "full", "trace", "work", "vina", "probe"):
        p.add_argument("--" + f, type=Path, required=True)
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    meta = json.loads((a.dev / "metadata.json").read_text())
    rows = []
    with (a.full / "ranked.csv").open() as f:
        ranked = list(csv.DictReader(f))
    for query, rec in [("5EPS", "5C7N"), ("3WKE", "4OD0")]:
        center = np.array(meta[rec]["center"])
        fields, origin, spacing = receptor_maps(
            a.vina, a.dev / rec / "receptor.pdbqt", center, a.work / "maps" / rec
        )
        inp = a.dev / query / "ligand.pdbqt"
        _, lines, _, heavy, types = read_ligand(a.probe, inp)
        serials = [lines[i][1] for i in heavy]
        for seed in [43103, 43117] if query == "5EPS" else [43117]:
            select = [
                r
                for r in ranked
                if r["query"] == query
                and r["arm"] == "lattice"
                and r["seed"] == str(seed)
                and r["rank"] == "1"
            ]
            (selected,) = select
            names = {
                "blind_winner": a.full
                / "runs"
                / f"{rec}_{query}_{seed}_lattice"
                / (selected["slot"] + ".pdbqt"),
                "best_close_start": a.trace
                / f"{query}_{rec}_{seed}_lattice"
                / "best_energy_within2.pdbqt",
            }
            for role, initial in names.items():
                for state, path in [
                    ("start", initial),
                    ("refined", Path(str(initial) + ".refined")),
                ]:
                    xyz = parse_pose(path, serials)
                    indices = (xyz - origin) / spacing
                    grid = 0.0
                    scalar = 0.0
                    for typ, idx in zip(types, indices):
                        field = fields[typ]
                        grid += float(
                            map_coordinates(
                                field, idx.reshape(3, 1), order=1, prefilter=False
                            )[0]
                        )
                        lo = idx.astype(int)
                        frac = idx - lo
                        for off in itertools.product([0, 1], repeat=3):
                            weight = np.prod(np.where(off, frac, 1 - frac))
                            scalar += weight * field[tuple(lo + off)]
                    if abs(grid - scalar) > 1e-8:
                        raise ValueError("trilinear disagreement")
                    cmd = [
                        str(a.probe),
                        "fixed",
                        str(a.dev / rec / "receptor.pdbqt"),
                        str(inp),
                        *map(str, center),
                        str(path),
                    ]
                    output = subprocess.check_output(cmd, text=True)
                    (uncapped,) = [
                        l.split()
                        for l in output.splitlines()
                        if l.startswith("uncapped ")
                    ]
                    (caps,) = [
                        l.split() for l in output.splitlines() if l.startswith("fixed ")
                    ]
                    unc_int, unc_intra = map(float, uncapped[1:])
                    cap_int, cap_intra = map(float, caps[2:])
                    rows.append(
                        {
                            "query": query,
                            "receiver": rec,
                            "seed": seed,
                            "arm": "lattice",
                            "role": role,
                            "state": state,
                            "coarse_inter": grid,
                            "explicit_native_inter": cap_int,
                            "explicit_native_intra": cap_intra,
                            "explicit_native_total": cap_int + cap_intra,
                            "uncapped_inter": unc_int,
                            "uncapped_intra": unc_intra,
                            "uncapped_total": unc_int + unc_intra,
                            "interp_diff_vs_uncapped": grid - unc_int,
                        }
                    )
    save(a.work / "fixed_scores.csv", rows)


if __name__ == "__main__":
    main()
