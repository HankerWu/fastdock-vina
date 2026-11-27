"""Measure crystal-mate compatibility of frozen VHL coordinates in two lattices."""
import argparse
import itertools
import json
import re
import sys
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

STUDY = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(STUDY))
from evaluate import PDBS, write_csv

CLASH = 2.3  # Å, as in benchmarks/crystal_environment/packing.py
PADDING = 18.0


def registration(source, receiver):
    keys = sorted(source["ca"].keys() & receiver["ca"].keys())
    a, b = [np.array([m["ca"][k] for k in keys]) for m in (source, receiver)]
    u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    rot = (u * [1, 1, np.linalg.det(u @ vt)]) @ vt
    shift = b.mean(0) - a.mean(0) @ rot
    assert len(keys) == 55 and np.linalg.det(rot) > 0.999999
    return rot, shift


def heavy_modes(path, serials):
    for block in path.read_text().split("ENDMDL"):
        if not block.strip():
            continue
        pos = {int(line[6:11]): [float(line[i:i+8]) for i in (30, 38, 46)]
               for line in block.splitlines() if line.startswith(("ATOM  ", "HETATM"))}
        score = re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", block)
        yield np.array([pos[i] for i in serials]), float(score[1]) if score else ""


def inputs(work, meta):
    """In the original receiver frame, under the original symmetry-aware labels."""
    rows = []
    for group, paths in (("blind", (work / "runs").glob("*.pdbqt")),
                         ("explicit_graph", work.glob("*_*_*_*.pdbqt"))):
        for path in sorted(paths):
            q, r, seed, arm = path.stem.split("_", 3)
            for rank, (xyz, score) in enumerate(heavy_modes(path, meta[q]["serial"]), 1):
                rows.append((group, q, r, seed, arm, rank, score, xyz))
    for q, r in itertools.permutations(PDBS):
        for group, suffix in [("witness", "fixed_local"), ("witness_fit", "fit")]:
            xyz, _ = next(heavy_modes(work / "witness" / f"{q}_{suffix}.pdbqt", meta[q]["serial"]))
            rows.append((group, q, r, "", "fixed", 1, "", xyz))
        rot, shift = registration(meta[q], meta[r])
        rows.append(("crystal_transfer", q, r, "", "native_C", 1, "",
                     np.array(meta[q]["crystal"][0]) @ rot + shift))
        rows.append(("deposited", q, q, "", "native_C", 1, "", np.array(meta[q]["crystal"][0])))
    return rows


class Crystal:
    """Enumerate symmetry images and translations intersecting the observed coordinate box."""
    def __init__(self, path, query):
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()
        atoms, records, host = [], [], []
        for chain in structure[0]:
            for res in chain:
                if res.entity_type != gemmi.EntityType.Polymer:
                    continue
                for a in res:
                    if a.occ <= 0 or a.element.is_hydrogen:
                        continue
                    atoms.append(list(a.pos))
                    records.append((chain.name, str(res.seqid), res.name, a.name,
                                    a.altloc.strip("\x00 "), round(a.occ, 3)))
                    host.append(chain.name in ("A", "B", "C"))
        ortho = np.asarray(structure.cell.orth.mat)
        frac = np.linalg.inv(ortho)
        query_f = query @ frac.T
        atoms_f = np.array(atoms) @ frac.T
        margin_f = PADDING * np.linalg.norm(frac, axis=1)
        lo, hi = query.min(axis=0) - PADDING, query.max(axis=0) + PADDING
        coords = dict(asu=[], image=[])
        labels = dict(asu=[], image=[])
        for op in structure.find_spacegroup().operations():
            a_f = (atoms_f @ np.array(op.rot).T + np.array(op.tran)) / op.DEN
            lower = np.ceil(query_f.min(0) - margin_f - a_f.max(0)).astype(int)
            upper = np.floor(query_f.max(0) + margin_f - a_f.min(0)).astype(int)
            for t in itertools.product(*(range(l, h+1) for l, h in zip(lower, upper))):
                xyz = (a_f + t) @ ortho.T
                mask = ((xyz >= lo) & (xyz <= hi)).all(axis=1)
                identity = op.triplet() == "x,y,z" and not any(t)
                part = "asu" if identity else "image"
                if identity:
                    mask &= ~np.array(host)
                ids = np.flatnonzero(mask)
                coords[part].extend(xyz[ids])
                image = op.triplet() + "|" + ",".join(map(str, t))
                labels[part].extend((*records[i], image) for i in ids)
        self.labels = labels
        self.trees = {k: cKDTree(v) for k, v in coords.items() if v}

    def measure(self, xyz, serials):
        ligand_tree = cKDTree(xyz)
        stats = {}
        for part, tree in self.trees.items():
            dist, index = tree.query(xyz)
            i = int(np.argmin(dist))
            closest = dist[i]
            assert closest < PADDING, "enlarge the image box before claiming this distance"
            chain, number, res, name, alt, occ, image = self.labels[part][index[i]]
            pairs = ligand_tree.sparse_distance_matrix(tree, CLASH, output_type="ndarray")
            # Alternate atom positions count as separate deposited records.
            stats.update({
                f"{part}_min_A": float(closest),
                f"{part}_pairs_lt_2p3": int((pairs["v"] < CLASH).sum()),
                f"{part}_ligand_serial": serials[i],
                f"{part}_atom": f"{chain}:{res}{number}:{name}",
                f"{part}_altloc": alt, f"{part}_occupancy": occ, f"{part}_op": image,
            })
        stats["min_A"] = min(stats[f"{part}_min_A"] for part in self.trees)
        stats["clashes"] = sum(stats[f"{part}_pairs_lt_2p3"] for part in self.trees)
        return stats


def main(work, out):
    meta = {p: json.loads((work / p / "metadata.json").read_text()) for p in PDBS}
    frames = {(q, r): registration(meta[q], meta[r]) for q, r in itertools.permutations(PDBS)}
    cases = inputs(work, meta)
    by_lattice = {p: [] for p in PDBS}
    for group, q, r, seed, arm, rank, score, xyz in cases:
        ref = np.asarray(meta[q]["crystal"])
        if q != r:
            rot, shift = frames[q, r]
            ref = ref @ rot + shift
        rmsd = float(np.sqrt(np.square(xyz - ref).sum(axis=2).mean(axis=1).min()))
        row = dict(group=group, query=q, receiver=r, seed=seed, arm=arm, rank=rank,
                   rmsd_A=rmsd, printed_score=score)
        by_lattice[r].append((dict(row, frame="receiver", lattice=r), xyz))
        if q != r:
            rot, shift = frames[r, q]
            by_lattice[q].append((dict(row, frame="query", lattice=q), xyz @ rot + shift))
    all_rows = []
    for p, data in by_lattice.items():
        crystal = Crystal(work / "raw" / f"{p}.pdb", np.concatenate([xyz for _, xyz in data]))
        for row, xyz in data:
            all_rows.append(dict(**row, **crystal.measure(xyz, meta[row["query"]]["serial"])))
    cells = []
    def key(r):
        return tuple(r[k] for k in ("group", "query", "receiver", "seed", "arm", "frame"))
    blind = [r for r in all_rows if r["group"] in ("blind", "explicit_graph")]
    for _, pool in itertools.groupby(sorted(blind, key=key), key):
        ordered = sorted(pool, key=lambda r: r["rank"])
        eligible = [r for r in ordered if not r["clashes"]]
        first = ordered[0]
        cells.append({**{k: first[k] for k in ("group", "query", "receiver", "seed", "arm", "frame")},
                      "modes": len(ordered), "compatible": len(eligible),
                      "native_A": first["rmsd_A"], "native_min_A": first["min_A"],
                      "native_clashes": first["clashes"],
                      "compatible_rank": eligible[0]["rank"] if eligible else "",
                      "compatible_top_A": eligible[0]["rmsd_A"] if eligible else "",
                      "oracle_A": min(r["rmsd_A"] for r in ordered),
                      "compatible_oracle_A": min((r["rmsd_A"] for r in eligible), default="")})
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "coordinates.csv", all_rows)
    write_csv(out / "pools.csv", cells)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    main(args.work, args.out)
