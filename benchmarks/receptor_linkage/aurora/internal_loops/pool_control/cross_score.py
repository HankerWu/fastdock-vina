"""Re-evaluate every released pose on every restored and original rigid receiver."""

import argparse
import hashlib
import io
import itertools
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay
from summarize import read, save

ARMS = tuple(replay.PROTEINS)


def collect(work, out, drug):
    ligand, _, serials, refs = replay.parameters(work, drug)
    rows = []
    known = {
        (r["arm"], r["seed"]): r
        for fold in ("base", "templates")
        for r in read(work / fold / "cells.csv")
        if r["drug"] == drug and r["arm"] in ARMS
    }
    for source, seed in itertools.product(ARMS, replay.geom.SEEDS):
        origin = "base" if source == "dry" else "templates"
        pool = work / origin / "blind" / f"{source}_{drug}_{seed}.pdbqt"
        assert replay.digest(pool) == known[source, str(seed)]["pool_sha256"]
        coord = replay.geom.mode_coords(pool, serials)
        blocks = [b for b in pool.read_text().split("ENDMDL") if b.strip()]
        for rank, (b, xyz) in enumerate(zip(blocks, coord, strict=True), 1):
            label = f"{drug}_{seed}_{source}_{rank}"
            path = out / "poses" / (label + ".pdbqt")
            content = (
                "\n".join(
                    l for l in b.strip().splitlines() if not l.startswith("MODEL")
                )
                + "\n"
            )
            path.write_text(content)
            coord_text = "\n".join(
                l[6:11] + l[30:54]
                for l in content.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            )
            rows.append(
                {
                    "drug": drug,
                    "seed": seed,
                    "source": source,
                    "rank": rank,
                    "label": label,
                    "coord_sha256": hashlib.sha256(coord_text.encode()).hexdigest(),
                    "rmsd_A": replay.mapped_rmsd(xyz, refs),
                    "xyz": xyz,
                    "path": path,
                }
            )
    return ligand, rows


def cross(work, probe, out):
    out.mkdir(parents=True, exist_ok=True)
    (out / "poses").mkdir(exist_ok=True)
    raw = out / "raw_scores"
    raw.mkdir(exist_ok=True)
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    results, costs = [], []
    dry = (work / replay.PROTEINS["dry"]).read_bytes()
    for drug in replay.geom.REFERENCES:
        ligand, rows = collect(work, out, drug)
        # Vina rounds 26 Å up to 70 voxels at 0.375 Å. The explicit receptor
        # then uses the actual coordinates (no clamping) for every frozen heavy atom.
        half_span = np.ceil(26 / 0.375) * 0.375 / 2
        assert max(np.max(np.abs(r["xyz"] - center)) for r in rows) <= half_span
        paths = "".join(str(r["path"]) + "\n" for r in rows)
        dry_scores = None
        for arm, ref in replay.PROTEINS.items():
            receptor = work / ref
            contents = receptor.read_bytes()
            assert contents.startswith(dry)
            added = [
                replay.geom.xyz(l)
                for l in contents[len(dry) :].decode().splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
                and l.split()[-1] not in ("H", "HD")
            ]
            stamp = time.monotonic()
            run = subprocess.run(
                [probe, receptor, ligand, *map(str, center)],
                input=paths,
                text=True,
                capture_output=True,
                check=True,
            )
            costs.append(
                {
                    "drug": drug,
                    "receptor": arm,
                    "pose_count": len(rows),
                    "wall_s": time.monotonic() - stamp,
                }
            )
            (raw / f"{drug}_{arm}.csv").write_text(
                "zero_score,inter,intra\n" + run.stdout
            )
            energies = np.loadtxt(io.StringIO(run.stdout), delimiter=",", ndmin=2)
            assert len(rows) == len(energies) and np.isfinite(energies).all()
            if arm == "dry":
                dry_scores = energies
            for row, energy, baseline in zip(rows, energies, dry_scores, strict=True):
                dist = float(cdist(row["xyz"], added).min()) if added else None
                delta = float(energy[0] - baseline[0])
                assert energy[2] == baseline[2]
                if dist is not None and dist >= 8:
                    assert abs(delta) < 1e-12
                attrs = {k: v for k, v in row.items() if k not in ("xyz", "path")}
                results.append(
                    dict(
                        **attrs,
                        receptor=arm,
                        zero_score=energy[0],
                        inter=energy[1],
                        intra=energy[2],
                        delta_vs_dry=delta,
                        nearest_added_A=dist,
                    )
                )
    save(out / "scores.csv", results)
    save(out / "costs.csv", costs)
    pools = []
    for drug, seed in itertools.product(replay.geom.REFERENCES, replay.geom.SEEDS):
        matches = [r for r in results if r["drug"] == drug and r["seed"] == seed]
        for source in (*ARMS, "union"):
            selected = [
                r for r in matches if source == "union" or r["source"] == source
            ]
            dry_winners = set()
            dry_gap = 0.0
            for arm in ARMS:
                pool = [r for r in selected if r["receptor"] == arm]
                value = min(r["zero_score"] for r in pool)
                winners = [r for r in pool if r["zero_score"] == value]
                hashes = {r["coord_sha256"] for r in winners}
                if arm == "dry":
                    dry_winners = hashes
                    alternatives = [
                        r["zero_score"] for r in pool if r["coord_sha256"] not in hashes
                    ]
                    dry_gap = min(alternatives) - value if alternatives else 0.0
                deltas = [r["delta_vs_dry"] for r in pool]
                oscillation = max(deltas) - min(deltas)
                pools.append(
                    {
                        "drug": drug,
                        "seed": seed,
                        "source": source,
                        "receptor": arm,
                        "n": len(pool),
                        "n_distinct_coords": len({r["coord_sha256"] for r in pool}),
                        "min_score": value,
                        "dry_gap_to_next_coordinate": dry_gap,
                        "score_change_range": oscillation,
                        "gap_certifies_stability": (
                            len(dry_winners) == 1 and dry_gap > oscillation
                        ),
                        "winners": "|".join(r["label"] for r in winners),
                        "n_winner_coords": len(hashes),
                        "same_winner_set": (hashes == dry_winners),
                        "disjoint_from_dry": hashes.isdisjoint(dry_winners),
                        "winner_rmsd_min_A": min(r["rmsd_A"] for r in winners),
                        "winner_rmsd_max_A": max(r["rmsd_A"] for r in winners),
                        "oracle_A": min(r["rmsd_A"] for r in pool),
                        "original_top_is_winner": any(r["rank"] == 1 for r in winners)
                        if source != "union"
                        else "",
                    }
                )
    save(out / "pools.csv", pools)
    summary = []
    for arm in ARMS[1:]:
        for drug in ("both", *replay.geom.REFERENCES):
            obs = [
                r
                for r in results
                if r["receptor"] == arm and drug in ("both", r["drug"])
            ]
            cells = [
                r for r in pools if r["receptor"] == arm and drug in ("both", r["drug"])
            ]
            individual = [r for r in cells if r["source"] != "union"]
            merged = [r for r in cells if r["source"] == "union"]
            summary.append(
                {
                    "receptor": arm,
                    "drug": drug,
                    "modes": len(obs),
                    "inside_8A": sum(r["nearest_added_A"] < 8 for r in obs),
                    "nonzero_delta": sum(r["delta_vs_dry"] != 0 for r in obs),
                    "nearest_A": min(r["nearest_added_A"] for r in obs),
                    "min_delta": min(r["delta_vs_dry"] for r in obs),
                    "max_delta": max(r["delta_vs_dry"] for r in obs),
                    "frozen_pools": len(individual),
                    "gap_certified_pools": sum(
                        r["gap_certifies_stability"] for r in individual
                    ),
                    "frozen_winner_changes": sum(
                        not r["same_winner_set"] for r in individual
                    ),
                    "union_pools": len(merged),
                    "gap_certified_unions": sum(
                        r["gap_certifies_stability"] for r in merged
                    ),
                    "union_winner_changes": sum(
                        not r["same_winner_set"] for r in merged
                    ),
                }
            )
    save(out / "summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for a in ("work", "probe", "out"):
        p.add_argument(a, type=Path)
    a = p.parse_args()
    cross(a.work.resolve(), a.probe.resolve(), a.out.resolve())
