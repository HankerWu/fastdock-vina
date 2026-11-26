"""Compare observed map support to a non-query donor map, donor geometry and an apo control."""

import argparse
import json
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from evaluate import CHAIN, ROOT, instances, models, read, save
from packing import atoms


def transform(source, target, pdb, ref):
    def positions(path, chain):
        by_id = {}
        for a in atoms(path):
            if (
                not a.polymer
                or a.name != "CA"
                or a.chain != chain
                or not a.number.isdigit()
            ):
                continue
            n = int(a.number)
            if not (pdb.startswith("6W") or pdb == "3KKX" or 44 <= n <= 166):
                continue
            if n not in by_id or a.occ > by_id[n][0]:
                by_id[n] = (a.occ, a.xyz)
        return by_id

    sa, tb = (
        positions(source, CHAIN.get(pdb, "A")),
        positions(target, CHAIN.get(ref, "A")),
    )
    ids = sorted(sa.keys() & tb.keys())
    a = np.array([sa[k][1] for k in ids])
    b = np.array([tb[k][1] for k in ids])
    if len(ids) < 75:
        raise ValueError((pdb, ref, len(ids)))
    u, _, vt = np.linalg.svd((a - a.mean(0)).T @ (b - b.mean(0)))
    return a.mean(0), u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt, b.mean(0)


def load_grid(folder, pdb):
    code = pdb.lower()
    name = f"{code}_validation_2fo-fc_map_coef.cif.gz"
    f = folder / name
    if not f.exists():
        url = f"https://files.rcsb.org/pub/pdb/validation_reports/{code[1:3]}/{code}/{name}"
        f.write_bytes(urlopen(url, timeout=60).read())
    doc = gemmi.cif.read(str(f))
    block = gemmi.as_refln_blocks(doc)[0]
    grid = block.transform_f_phi_to_map("pdbx_FWT", "pdbx_PHWT", sample_rate=3.0)
    grid.normalize()
    return grid


def evaluate(ca, brd, ca_runs, brd_runs, results, maps):
    cache = {}
    data = read(results / "candidates.csv")
    meta = json.loads((brd / "prepared/manifest.json").read_text())
    graph = {
        p: json.loads(
            (ROOT / f"benchmarks/mapping_ensemble/results/{p}/graphs.json").read_text()
        )
        for p in ("brd4", "hcaii_declared")
    }
    rows = []
    observed = []
    for target in ("hcaii", "brd4"):
        root = ca if target == "hcaii" else brd
        work = ca_runs if target == "hcaii" else brd_runs
        g = graph["hcaii_declared" if target == "hcaii" else "brd4"]
        apo = "3KS3" if target == "hcaii" else "2OSS"
        coordinates = lambda k: (
            maps / (k + ".pdb")
            if k in ("2OSS", "3KS3")
            else root / "raw" / (k + ".pdb")
        )
        if not coordinates(apo).exists():
            coordinates(apo).write_bytes(
                urlopen(f"https://files.rcsb.org/download/{apo}.pdb", timeout=60).read()
            )
        bins = defaultdict(list)
        for r in data:
            if r["target"] == target:
                bins[r["receptor"], r["ligand"], r["seed"], r["slot"]].append(r)
        for pdb in list(g) + [apo]:
            grid = load_grid(maps, pdb)
            structure = atoms(coordinates(pdb))
            host = CHAIN.get(pdb, "A")
            ligids = instances(structure, pdb) if pdb != apo else {}
            active = next(
                (
                    group
                    for (ch, num), group in ligids.items()
                    if ch == host and (target != "hcaii" or num == "303")
                ),
                [],
            )
            if active:
                lig_xyz = np.array([a.xyz for a in active])
                lig_density = grid.interpolate_position_array(lig_xyz)
                for a, den in zip(active, lig_density):
                    observed.append(
                        dict(
                            target=target,
                            pdb=pdb,
                            name=a.name,
                            element=a.element,
                            b=a.b,
                            normalized_2mFo_DFc_density=float(den),
                        )
                    )
                lig_tree = cKDTree(lig_xyz)
            else:
                lig_tree = None
            contaminants = cKDTree(
                [
                    a.xyz
                    for a in structure
                    if a.polymer or a.residue in ("ZN", "HOH", "DOD")
                ]
            )
            cache[pdb] = grid, lig_tree, contaminants
        for (rec, lig, seed, slot), group in bins.items():
            coords = models(
                work / f"{rec}_{lig}_{seed}_{slot}.pdbqt",
                meta[lig]["heavy_serials"] if target == "brd4" else None,
            )
            refs = dict(own=lig, donor=g[lig]["donor"], apo=apo)
            fits = {
                k: transform(coordinates(rec), coordinates(id), rec, id)
                for k, id in refs.items()
            }
            for r in group:
                xyz, _ = coords[int(r["model"])]
                record = {
                    k: r[k]
                    for k in (
                        "target",
                        "receptor",
                        "ligand",
                        "seed",
                        "pool",
                        "slot",
                        "model",
                        "rmsd",
                        "sign",
                        "energy",
                        "intra",
                        "native",
                        "factor",
                    )
                }
                for kind, id in refs.items():
                    o, R, d = fits[kind]
                    probe = np.asarray((xyz - o) @ R + d, order="C")
                    grid, lig_tree, contam = cache[id]
                    values = grid.interpolate_position_array(probe)
                    record[kind + "_density"] = float(values.mean())
                    if lig_tree is not None:
                        dlig = lig_tree.query(probe)[0]
                        # Mask protein/solvent peaks and distant background to isolate ligand-supported volume.
                        mask = (dlig < contam.query(probe)[0]) & (dlig <= 2)
                        record[kind + "_masked"] = float(
                            np.mean(np.where(mask, np.maximum(values, 0), 0))
                        )
                        record[kind + "_shape"] = float(np.exp(-dlig * dlig / 2).mean())
                rows.append(record)
    save(results / "observed_density.csv", observed)
    save(results / "density_candidates.csv", rows)
    groups = defaultdict(list)
    for r in rows:
        groups[r["target"], r["receptor"], r["ligand"], r["seed"], r["pool"]].append(r)
    cells = []
    methods = (
        "own_density",
        "donor_density",
        "apo_density",
        "own_masked",
        "donor_masked",
        "own_shape",
        "donor_shape",
    )
    for (target, rec, lig, seed, arm), rs in groups.items():
        native = min(rs, key=lambda r: float(r["native"]))
        for method in methods:
            best = min(
                rs,
                key=lambda r: (
                    -r[method],
                    float(r["native"]),
                    r["slot"] != "common",
                    int(r["model"]),
                ),
            )
            energy = (
                (
                    (float(best["native"]) - float(native["intra"]))
                    * float(best["factor"])
                )
                if target == "hcaii"
                else float(best["energy"])
            )
            cells.append(
                dict(
                    target=target,
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    pool=arm,
                    measure=method,
                    model=best["model"],
                    slot=best["slot"],
                    rmsd=best["rmsd"],
                    sign=best["sign"],
                    joint=int(float(best["rmsd"]) <= 2 and int(best["sign"])),
                    energy=energy,
                )
            )
    save(results / "density_cells.csv", cells)
    group = defaultdict(list)
    for r in cells:
        group[r["target"], r["pool"], r["measure"]].append(r)
    save(
        results / "density_summary.csv",
        [
            dict(
                target=t,
                pool=a,
                method=m,
                n=len(rs),
                rmsd_top=sum(float(r["rmsd"]) <= 2 for r in rs),
                joint_top=sum(r["joint"] for r in rs),
                median_rmsd=np.median([float(r["rmsd"]) for r in rs]),
            )
            for (t, a, m), rs in sorted(group.items())
        ],
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("ca", "brd", "ca_runs", "brd_runs", "results", "maps"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.maps.mkdir(exist_ok=True, parents=True)
    evaluate(a.ca, a.brd, a.ca_runs, a.brd_runs, a.results, a.maps)
