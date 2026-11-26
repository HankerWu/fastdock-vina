"""Test whether receptor-specific packing or choice of experimental copy explains pose errors."""

import argparse
import csv
import gzip
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET
from urllib.request import urlopen

import numpy as np
from scipy.spatial.transform import Rotation

from packing import atoms, crystal_neighbors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from select_zinc_pose import parse_atoms, read_models

COMPONENT = dict(
    zip(
        (
            "6WQ4",
            "6WQ5",
            "6WQ7",
            "6WQ8",
            "6WQ9",
            "4OGI",
            "4OGJ",
            "4O76",
            "4O77",
            "3MXF",
        ),
        ("U7A", "U7G", "U6V", "U77", "U7J", "R78", "2TA", "1M3", "2RE", "JQ1"),
        strict=True,
    )
)
CHAIN = {k: ("B" if k == "4OGJ" else "A") for k in COMPONENT}


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def instances(rows, code):
    by_id = defaultdict(list)
    for a in rows:
        if a.residue == COMPONENT[code]:
            by_id[a.chain, a.number].append(a)
    return by_id


def transformed_backbone(rows, chain, target_rows, target_chain, xyz):
    def ca(atoms, chain):
        return {
            int(a.number): np.array(a.xyz)
            for a in atoms
            if a.polymer
            and a.name == "CA"
            and a.chain == chain
            and a.alt in ("", "A")
            and a.number.isdigit()
            and 44 <= int(a.number) <= 166
        }

    left, right = ca(rows, chain), ca(target_rows, target_chain)
    keys = sorted(left.keys() & right.keys())
    if len(keys) < 75:
        raise ValueError("Too few common BRD4 protein CA atoms")
    a, b = np.array([left[i] for i in keys]), np.array([right[i] for i in keys])
    rot, _ = Rotation.align_vectors(b - b.mean(0), a - a.mean(0))
    return rot.apply(xyz - a.mean(0)) + b.mean(0)


def structure_summary(path, pdb, validation):
    seq = atoms(path)
    dest = validation / (pdb.lower() + "_validation.xml.gz")
    url = f"https://files.rcsb.org/pub/pdb/validation_reports/{pdb[1:3].lower()}/{pdb.lower()}/{dest.name}"
    if not dest.exists():
        dest.write_bytes(urlopen(url, timeout=60).read())
    info = ET.fromstring(gzip.decompress(dest.read_bytes()))
    resolution = float(info.find("Entry").get("PDB-resolution"))
    rows = []
    for (chain, num), group in sorted(instances(seq, pdb).items()):
        xyz = np.array([a.xyz for a in group])
        env = crystal_neighbors(path, chain, xyz)
        geom = env.contacts(xyz, [a.element for a in group])
        ndx = geom.pop("ligand_nearest")
        zn = np.array(
            [a.xyz for a in seq if a.residue == "ZN" and a.chain == chain]
        ).reshape(-1, 3)
        distance_to_zn = (
            float(np.min(np.linalg.norm(xyz[:, None] - zn, axis=-1))) if len(zn) else ""
        )
        (sg,) = [
            a
            for a in info.findall(".//ModelledSubgroup")
            if a.get("resname") == COMPONENT[pdb]
            and a.get("chain") == chain
            and a.get("resnum") == num
        ]
        rows.append(
            dict(
                pdb=pdb,
                component=COMPONENT[pdb],
                chain=chain,
                residue=num,
                resolution=resolution,
                atoms=len(group),
                occ_min=min(a.occ for a in group),
                altlocs=" ".join(sorted({a.alt for a in group})),
                b_mean=np.mean([a.b for a in group]),
                rscc=sg.get("rscc"),
                rsr=sg.get("rsr"),
                min_zinc=distance_to_zn,
                ligand_atom=group[ndx].name,
                **geom,
            )
        )
    return rows


def brd_refs(raw, meta, lig, rec):
    """Map deposit atom names into the reference's original order, then generate all graph-equivalent references."""
    seq = atoms(raw / (lig + ".pdb"))
    target = atoms(raw / (rec + ".pdb"))
    groups = instances(seq, lig)
    (chosen,) = [rs for (chain, _), rs in groups.items() if chain == CHAIN[lig]]
    src = np.array([a.xyz for a in chosen])
    names = np.array([a.name for a in chosen])
    src = transformed_backbone(seq, CHAIN[lig], target, CHAIN[rec], src)
    old = np.asarray(meta["alignments"][rec]["reference_by_match"])
    # Strict Cartesian match to the deposited source makes the audit independent of SDF atom numbering.
    distances = np.linalg.norm(old[0, :, None, :] - src[None, :, :], axis=-1)
    order = distances.argmin(axis=1)
    if len(set(order)) != len(order) or distances.min(axis=1).max() > 1e-4:
        raise ValueError("Wrong source correspondence")
    old_order = np.linalg.norm(old[:, :, None, :] - old[0, None, :, :], axis=-1).argmin(
        axis=2
    )
    out = []
    for (chain, num), rs in sorted(groups.items()):
        coords = {a.name: a.xyz for a in rs}
        xyz = transformed_backbone(
            seq, chain, target, CHAIN[rec], np.array([coords[k] for k in names[order]])
        )
        out.append((chain + "_" + num, xyz[old_order]))
    return out


def models(path, serials=None):
    converted = {}
    for num, lines in read_models(path.read_text().splitlines()).items():
        pts = parse_atoms(lines)
        ids = serials or sorted(
            s for s, a in pts.items() if a.atom_type not in ("H", "HD", "HS")
        )
        typemap = {"A": "C", "NA": "N", "OA": "O", "SA": "S"}
        converted[num] = (
            np.array([pts[s].xyz for s in ids]),
            [typemap.get(pts[s].atom_type, pts[s].atom_type).upper() for s in ids],
        )
    return converted


def evaluate(data, brd, ca_runs, brd_runs, out):
    out.mkdir(parents=True, exist_ok=True)
    validation = out / "validation"
    validation.mkdir(exist_ok=True)
    struct_rows = []
    for code in COMPONENT:
        raw = (data if code.startswith("6W") else brd) / "raw"
        struct_rows.extend(structure_summary(raw / (code + ".pdb"), code, validation))
    save(out / "crystals.csv", struct_rows)
    refs = json.loads((brd / "prepared/manifest.json").read_text())
    hmeta = json.loads((data / "prepared/metadata.json").read_text())
    ca = read(ca_runs / "candidates.csv")
    ca_prior = {
        (r["receptor"], r["ligand"], r["seed"], r["slot"], r["model"]): float(
            r["prior"]
        )
        for r in read(ca_runs / "rank_energies.csv")
        if r["arm"] == "both_states"
    }
    ca = [r for r in ca if r["slot"] in ("common", "second_reflected")]
    brows = read(brd_runs / "modes.csv")
    candidates = []
    cache = {}
    for target, rows, raw, work in (
        ("hcaii", ca, data / "raw", ca_runs),
        ("brd4", brows, brd / "raw", brd_runs),
    ):
        for rec in sorted({r["receptor"] for r in rows}):
            by_file = defaultdict(list)
            for r in rows:
                if r["receptor"] != rec:
                    continue
                slot = r["slot"] if target == "hcaii" else r["arm"]
                by_file[(r["ligand"], r["seed"], slot)].append(r)
            packed = {}
            for lig, seed, slot in by_file:
                packed[lig, seed, slot] = models(
                    work / f"{rec}_{lig}_{seed}_{slot}.pdbqt",
                    refs[lig]["heavy_serials"] if target == "brd4" else None,
                )
            xyz = np.concatenate(
                [xyz for group in packed.values() for xyz, _ in group.values()]
            )
            env = crystal_neighbors(raw / (rec + ".pdb"), CHAIN[rec], xyz)
            for (lig, seed, slot), group in by_file.items():
                if target == "brd4":
                    cache[rec, lig] = brd_refs(raw, refs[lig], lig, rec)
                for r in group:
                    n = int(r["model"])
                    xyz, els = packed[lig, seed, slot][n]
                    geom = env.contacts(xyz, els)
                    rms = float(r["rmsd"])
                    alternate = rms
                    altid = ""
                    if target == "brd4":
                        copies = cache[rec, lig]
                        fits = {
                            k: float(
                                np.sqrt(
                                    np.square(ref - xyz).sum(axis=2).mean(axis=1).min()
                                )
                            )
                            for k, ref in copies
                        }
                        (base,) = [v for k, v in fits.items() if k[0] == CHAIN[lig]]
                        if abs(rms - base) > 1e-4:
                            raise ValueError(f"Mismatch {lig} {rec} {rms} {base}")
                        altid = min(fits, key=fits.get)
                        alternate = fits[altid]
                    candidates.append(
                        dict(
                            target=target,
                            receptor=rec,
                            ligand=lig,
                            seed=seed,
                            pool="both_states" if target == "hcaii" else slot,
                            slot=slot,
                            model=n,
                            native=float(r["total"] if target == "hcaii" else r["raw"]),
                            energy=float(
                                r["selected_score"]
                                if target == "hcaii"
                                else r["native_affinity"]
                            ),
                            intra=float(r["intra"]) if target == "hcaii" else 0,
                            factor=float(r["flex_factor"]) if target == "hcaii" else 0,
                            prior=(
                                ca_prior[rec, lig, seed, slot, r["model"]]
                                if target == "hcaii"
                                else float(r["multi"])
                            ),
                            rmsd=rms,
                            alt_rmsd=alternate,
                            alt_id=altid,
                            sign=int(
                                r["same_sign_as_bound"]
                                if target == "hcaii"
                                else r["sign_ok"]
                            ),
                            **geom,
                        )
                    )
    save(out / "candidates.csv", candidates)
    blocks = defaultdict(list)
    for r in candidates:
        blocks[
            tuple(r[k] for k in ("target", "receptor", "ligand", "seed", "pool"))
        ].append(r)
    cells = []
    for (target, rec, lig, seed, pool), group in blocks.items():
        native = min(
            group, key=lambda r: (r["native"], r["slot"] != "common", r["model"])
        )
        if target == "brd4":
            # Native order uses affinity after explicit refinement.
            native = min(group, key=lambda r: r["model"])
        for alg in ("native", "native_pack", "prior", "prior_pack"):

            def raw(r):
                return (
                    r["native"]
                    + ("prior" in alg) * r["prior"]
                    + ("pack" in alg) * r["repulsion"]
                )

            top = (
                native
                if alg == "native"
                else min(
                    group, key=lambda r: (raw(r), r["slot"] != "common", r["model"])
                )
            )
            # Keep native binding energy at the selected pose for every arm. Packing is a crystal diagnostic.
            energy = (
                (top["native"] - native["intra"]) * top["factor"]
                if target == "hcaii"
                else top["energy"]
            )
            cells.append(
                dict(
                    target=target,
                    receptor=rec,
                    ligand=lig,
                    seed=seed,
                    pool=pool,
                    method=alg,
                    model=top["model"],
                    slot=top["slot"],
                    rmsd=top["rmsd"],
                    hit=int(top["rmsd"] <= 2 and top["sign"]),
                    alt_rmsd=top["alt_rmsd"],
                    alt_hit=int(top["alt_rmsd"] <= 2 and top["sign"]),
                    oracle=int(any(r["rmsd"] <= 2 and r["sign"] for r in group)),
                    alt_oracle=int(
                        any(r["alt_rmsd"] <= 2 and r["sign"] for r in group)
                    ),
                    energy=energy,
                    repulsion=top["repulsion"],
                )
            )
    save(out / "cells.csv", cells)
    # Retain pair-repeated / unique chemical comparisons. Solution labels are never assigned lattice energy.
    groups = defaultdict(list)
    for r in cells:
        groups[
            tuple(r[k] for k in ("target", "receptor", "seed", "pool", "method"))
        ].append(r)
    pairs = []
    for (target, rec, seed, pool, alg), group in groups.items():
        get = lambda r: (
            hmeta[r["ligand"]]["Ki_nM"]
            if target == "hcaii"
            else refs[r["ligand"]]["kd_uM"]
        )
        for a, b in itertools.combinations(sorted(group, key=get), 2):
            if get(b) < 3 * get(a):
                continue
            gap = b["energy"] - a["energy"]
            pairs.append(
                dict(
                    target=target,
                    receptor=rec,
                    seed=seed,
                    pool=pool,
                    method=alg,
                    strong=a["ligand"],
                    weak=b["ligand"],
                    gap=gap,
                    credit=1 if gap > 0 else 0.5 if gap == 0 else 0,
                )
            )
    save(out / "pairs.csv", pairs)
    summary = []
    for k in sorted({(r["target"], r["pool"], r["method"]) for r in cells}):
        g = [r for r in cells if (r["target"], r["pool"], r["method"]) == k]
        p = [r for r in pairs if (r["target"], r["pool"], r["method"]) == k]
        summary.append(
            dict(
                target=k[0],
                pool=k[1],
                method=k[2],
                cells=len(g),
                top=sum(r["hit"] for r in g),
                alt_top=sum(r["alt_hit"] for r in g),
                oracle=sum(r["oracle"] for r in g),
                alt_oracle=sum(r["alt_oracle"] for r in g),
                median=float(np.median([r["rmsd"] for r in g])),
                alt_median=float(np.median([r["alt_rmsd"] for r in g])),
                pair_credit=sum(r["credit"] for r in p),
                pairs=len(p),
            )
        )
    save(out / "summary.csv", summary)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for n in ("ca_data", "brd_data", "ca_runs", "brd_runs", "out"):
        p.add_argument(n, type=Path)
    a = p.parse_args()
    evaluate(a.ca_data, a.brd_data, a.ca_runs, a.brd_runs, a.out)
