"""Cross-check all distance/count rows via Gemmi's neighbor search and direct distances."""
import argparse
import csv
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial.distance import cdist


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def positions(xyz):
    return [gemmi.Position(*x) for x in xyz]


def transform(xyz, moving, fixed):
    keys = sorted(moving["ca"].keys() & fixed["ca"].keys())
    a, b = [[gemmi.Position(*m["ca"][k]) for k in keys] for m in (moving, fixed)]
    tr = gemmi.superpose_positions(b, a).transform
    return np.array([list(tr.apply(p)) for p in positions(xyz)])


def coords(row, work, meta, structs):
    q, r = row["query"], row["receiver"]
    group = row["group"]
    if group in ("deposited", "crystal_transfer"):
        ligand = next(res for res in structs[q][0]["C"] if res.name == meta[q]["component"])
        xyz = np.array([list(a.pos) for a in ligand if a.occ > 0 and not a.element.is_hydrogen])
        if q != r:
            xyz = transform(xyz, meta[q], meta[r])
    else:
        if group in ("blind", "explicit_graph"):
            folder = work / "runs" if group == "blind" else work
            name = "_".join(row[k] for k in ("query", "receiver", "seed", "arm"))
        else:
            folder = work / "witness"
            name = q + ("_fixed_local" if group == "witness" else "_fit")
        model = 1
        xyz = []
        for line in (folder / (name + ".pdbqt")).read_text().splitlines():
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            if (model == int(row["rank"]) and line.startswith(("ATOM", "HETATM"))
                    and line.split()[-1] not in ("HD", "H", "HS")):
                xyz.append([float(line[i:i+8]) for i in (30, 38, 46)])
        xyz = np.array(xyz)
    assert len(xyz) == len(meta[q]["serial"])
    if row["frame"] == "query":
        xyz = transform(xyz, meta[r], meta[q])
    return xyz


def environment(st, center, radius):
    # Every point in this sphere has a single nearest periodic image of a
    # specific space-group copy; the VHL cells are orthogonal with a,b >= 92 Å.
    assert 2 * radius < min(st.cell.a, st.cell.b, st.cell.c)
    assert [st.cell.alpha, st.cell.beta, st.cell.gamma] == [90, 90, 90]
    search = gemmi.NeighborSearch(st[0], st.cell, radius).populate(include_h=False)
    ref = gemmi.Position(*center)
    env = {"asu": [], "image": []}
    for mark in search.find_atoms(ref, radius=radius):
        cra = mark.to_cra(st[0])
        if cra.residue.entity_type != gemmi.EntityType.Polymer or cra.atom.occ <= 0:
            continue
        img = st.cell.find_nearest_pbc_image(ref, cra.atom.pos, mark.image_idx)
        if img.same_asu() and cra.chain.name in "ABC":
            continue
        part = "asu" if img.same_asu() else "image"
        xyz = st.cell.find_nearest_pbc_position(ref, cra.atom.pos, mark.image_idx)
        env[part].append(list(xyz))
    return {k: np.array(v) for k, v in env.items()}


def main(work, results):
    data = read_csv(results / "coordinates.csv")
    files = [(p, work / p / "metadata.json", work / "raw" / (p + ".pdb"))
             for p in ("4W9H", "5LLI")]
    meta = {p: json.loads(m.read_text()) for p, m, _ in files}
    structs = {p: gemmi.read_structure(str(s)) for p, _, s in files}
    for st in structs.values():
        st.setup_entities()
    errors = []
    for p, m in meta.items():
        rows = [row for row in data if row["lattice"] == p]
        xyz = [coords(row, work, meta, structs) for row in rows]
        center = np.array(m["center"])
        extent = np.linalg.norm(np.concatenate(xyz) - center, axis=1).max()
        env = environment(structs[p], center, extent + 18.0)
        for row, atoms in zip(rows, xyz):
            for part, cloud in env.items():
                separation = cdist(atoms, cloud)
                assert separation.min() < 18
                errors.append(abs(separation.min() - float(row[f"{part}_min_A"])))
                assert int((separation < 2.3).sum()) == int(row[f"{part}_pairs_lt_2p3"])
    assert max(errors) < 1e-8
    history = read_csv(Path(__file__).resolve().parent.parent / "results" / "modes.csv")
    baseline = {(r["query"], r["receptor"], r["seed"], r["arm"], r["rank"]): float(r["rmsd"])
                for r in history}
    check = [abs(float(r["rmsd_A"]) - baseline[r["query"], r["receiver"], r["seed"], r["arm"], r["rank"]])
             for r in data if r["group"] == "blind" and r["frame"] == "receiver"]
    assert len(check) == len(baseline) == 91 and max(check) < 1e-10
    summary = {"coordinates": len(data), "neighbor_distance_checks": len(errors),
               "neighbor_count_checks": len(errors), "max_distance_difference_A": max(errors),
               "original_blind_rmsd_checks": len(check), "max_blind_rmsd_difference_A": max(check)}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    main(args.work, args.results)
