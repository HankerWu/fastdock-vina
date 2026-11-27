"""Independent SciPy reconstruction of AD4 energies; cross-score identical released poses."""

import argparse
import csv
import itertools
from pathlib import Path
import sys

import numpy as np

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "cubic_grid/fkbp_ad4"))
from analyze import read_map
from run_official import table


def cap_energy(x, cap):
    out = x.copy()
    pos = x > 0
    out[pos] = cap * x[pos] / (cap + x[pos])
    return out


def components(model, fields):
    lines = model.splitlines()
    atoms = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    xyz = np.array([[float(l[i : i + 8]) for i in (30, 38, 46)] for l in atoms])
    q = np.array([float(l[70:76]) for l in atoms])
    kinds = np.array([l[77:].strip() for l in atoms])
    lo = [x[0] for x in fields["e"].grid]
    hi = [x[-1] for x in fields["e"].grid]
    clamp = np.clip(xyz, lo, hi)
    miss = np.abs(xyz - clamp).sum(axis=1)
    type_field = np.empty(len(atoms))
    for t in set(kinds):
        mask = kinds == t
        type_field[mask] = fields[t](clamp[mask])
    e = fields["e"](clamp)
    d = fields["d"](clamp)
    torsion = 0.2983 * int(next(l.split()[1] for l in lines if l.startswith("TORSDOF")))
    raw = type_field + q * e + abs(q) * d
    old = (
        cap_energy(type_field, 1000)
        + q * cap_energy(e, 1000)
        + abs(q) * cap_energy(d, 1000)
    )
    return dict(
        parent=float(sum(old + 1e6 * miss * (1 + q + abs(q))) + torsion),
        wall=float(sum(old + 1e6 * miss) + torsion),
        joint=float(sum(cap_energy(raw, 1000) + 1e6 * miss) + torsion),
        uncapped=float(sum(raw) + torsion),
        outside_A=max(miss),
    )


def rounded_bounds(model, fields, arm):
    """Bound a trilinear map over the rounding interval, including lattice planes."""
    lines = model.splitlines()
    torsion = 0.2983 * int(next(l.split()[1] for l in lines if l.startswith("TORSDOF")))
    total_lo = total_hi = torsion
    for line in (l for l in lines if l.startswith(("ATOM", "HETATM"))):
        xyz = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
        q = float(line[70:76])
        axes = fields["e"].grid
        left, right = xyz - 0.00050001, xyz + 0.00050001
        lo = np.array([g[0] for g in axes])
        hi = np.array([g[-1] for g in axes])
        mesh = [
            np.unique(
                np.r_[
                    np.clip(left[j], lo[j], hi[j]),
                    np.clip(right[j], lo[j], hi[j]),
                    g[(g > left[j]) & (g < right[j])],
                ]
            )
            for j, g in enumerate(axes)
        ]
        vertices = np.array(list(itertools.product(*mesh)))
        t = fields[line[77:].strip()](vertices)
        e = fields["e"](vertices)
        d = fields["d"](vertices)
        if arm == "joint":
            energy = cap_energy(t + q * e + abs(q) * d, 1000)
            atom_lo, atom_hi = energy.min(), energy.max()
        else:
            terms = [
                cap_energy(t, 1000),
                q * cap_energy(e, 1000),
                abs(q) * cap_energy(d, 1000),
            ]
            atom_lo, atom_hi = sum(x.min() for x in terms), sum(x.max() for x in terms)
        slope = 1e6 * (1 + q + abs(q) if arm in ("official", "parent") else 1)
        minimum = np.maximum.reduce([lo - right, left - hi, np.zeros(3)])
        maximum = np.maximum.reduce([lo - left, right - hi, np.zeros(3)])
        total_lo += atom_lo + slope * minimum.sum()
        total_hi += atom_hi + slope * maximum.sum()
    return total_lo, total_hi


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work", "panel", "out"):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    modes = {}
    for path in sorted(a.panel.glob("*_models.csv")):
        for r in csv.DictReader(path.open()):
            tag = "-".join(r[n] for n in ("query", "receiver", "seed", "cap", "arm"))
            modes[tag, int(r["rank"])] = r
    fields = {}
    rows = []
    rankings = []
    for path in sorted((a.panel / "runs").glob("*.pdbqt")):
        query, rec, seed, cap, arm = path.stem.split("-")
        if arm == "vina":
            continue
        if rec not in fields:
            maps = a.work / "maps" / rec
            fields[rec] = {
                f.name.split(".")[1]: read_map(f) for f in maps.glob("receptor.*.map")
            }
        pool = []
        for rank, model in enumerate(
            (x for x in path.read_text().split("ENDMDL") if x.strip()), 1
        ):
            comp = components(model, fields[rec])
            mode = modes[path.stem, rank]
            native = comp["parent" if arm == "official" else arm]
            row = dict(
                query=query,
                receiver=rec,
                seed=seed,
                cap=cap,
                arm=arm,
                rank=rank,
                rmsd=float(mode["rmsd"]),
                reported=float(mode["score"]),
                error=native - float(mode["score"]),
                **comp,
            )
            # Within 0.02 kcal/mol direct replay suffices. For steep repulsions and
            # walls, bound the full ±0.0005 Å coordinate box before judging a mismatch.
            if abs(row["error"]) > 0.02:
                low, high = rounded_bounds(model, fields[rec], arm)
                if not low - 0.000501 <= row["reported"] <= high + 0.000501:
                    raise ValueError((path, rank, low, row["reported"], high))
                row.update(rounding_low=low, rounding_high=high)
            else:
                row.update(rounding_low="", rounding_high="")
            rows.append(row)
            pool.append(row)
        if arm in ("parent", "joint"):
            best_old = min(pool, key=lambda r: r["parent"])
            best_new = min(pool, key=lambda r: r["joint"])
            rankings.append(
                dict(
                    query=query,
                    receiver=rec,
                    seed=seed,
                    cap=cap,
                    generated=arm,
                    parent_rank=best_old["rank"],
                    joint_rank=best_new["rank"],
                    parent_top1_rmsd=best_old["rmsd"],
                    joint_top1_rmsd=best_new["rmsd"],
                    pool_min=min(r["rmsd"] for r in pool),
                )
            )
    table(a.out / "fixed_pool.csv", rows)
    table(a.out / "fixed_pool_ranks.csv", rankings)


if __name__ == "__main__":
    main()
