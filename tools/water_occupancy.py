"""Prepare external water sites or marginalize their occupancies on a fixed pose pool."""

import argparse
import json
from pathlib import Path

import numpy as np
from extract_independent_receptor_evidence import (
    align_donor,
    aligned_waters,
    coordinate,
    dominant_atoms,
    structure_chain,
    water_sites,
)
from scipy.spatial.distance import cdist
from select_zinc_pose import parse_atoms, read_models

RT = 0.00198720425864 * 298.15
WATER_LIGAND_CONTACT = (
    -0.587439
)  # one native Vina H-bond depth, used as a single-contact proxy


def atom_arrays(lines):
    atoms = list(parse_atoms(lines).values())
    if not atoms:
        raise ValueError("PDBQT has no atoms")
    xyz = np.array([a.xyz for a in atoms])
    types = np.array([a.atom_type for a in atoms])
    heavy = ~np.isin(types, ["H", "HD", "HS", "W", "G0", "G1", "G2", "G3"])
    protic = np.zeros(len(types), dtype=bool)
    if np.any(types == "HD"):
        protic = np.isin(types, ["N", "NA", "O", "OA"]) & (
            cdist(xyz, xyz[types == "HD"]).min(axis=1) < 1.25
        )
    roles = np.isin(types, ["OA", "NA"]).astype(int) + 2 * protic.astype(int)
    return xyz[heavy], roles[heavy], types[heavy]


def build_sites(raw, rec, chain, receptor, center, size, donors):
    """Fit donor protein only. Inputs identify distinct experimental PDB entries and chains."""
    if len(donors) < 2:
        raise ValueError("At least two independent donor structures are required")
    _, q = structure_chain(raw, rec, chain)
    xyz, roles, types = atom_arrays(receptor.read_text().splitlines())
    metal = xyz[np.isin(types, ["Zn", "Mg", "Ca", "Fe", "Cu", "Mn", "Co", "Ni"])]
    reports, waters = [], []
    for name, donor_chain, path in sorted(donors):
        if name == rec or any(r["id"] == name for r in reports):
            raise ValueError("Water evidence must use distinct external PDB IDs")
        model, protein = structure_chain(path, name, donor_chain)
        report = {"id": name}
        _, rot, shift = align_donor(q, protein, center, report)
        obs = [
            w
            for w in aligned_waters(model, name, center, size, rot, shift)
            if w["chain"] == donor_chain
        ]
        waters.extend(obs)
        reports.append(
            {
                "id": name,
                "fit_rmsd": report["local_backbone_rmsd_angstrom"],
                "local_residues": report["local_residues"],
                "water_oxygens": len(obs),
            }
        )
    sites = []
    for cluster in water_sites(waters):
        count = cluster["independent_pdb_count"]
        point = np.array(cluster["xyz"])
        dist = np.linalg.norm(xyz - point, axis=1)
        nearby = (roles != 0) & (dist >= 2.4) & (dist <= 3.5)
        partners = xyz[nearby]
        if count < 2 or not len(partners) or dist.min() < 2.3:
            continue
        if np.any(np.linalg.norm(metal - point, axis=1) < 3.0):
            continue
        sites.append(
            {
                "xyz": point.tolist(),
                "partners": partners.tolist(),
                "partner_roles": roles[nearby].tolist(),
                "p": (count + 0.5) / (len(donors) + 1),
                "count": count,
                "donors": [w["pdb_id"] for w in cluster["hits"]],
                "source_residues": [w["source_resseq"] for w in cluster["hits"]],
            }
        )
    return {"receptor": rec, "chain": chain, "alignment": reports, "sites": sites}


def displacement_witnesses(raw, rec, chain, center, sites, sources):
    """Find sites occluded by specified ligands from other experimental structures."""
    _, reference = structure_chain(raw, rec, chain)
    support = [[] for _ in sites]
    site_points = np.array([site["xyz"] for site in sites], dtype=float).reshape(-1, 3)
    seen = set()
    for pdb, source_chain, component, path in sources:
        if pdb == rec or pdb in seen:
            raise ValueError("Displacement evidence must use distinct external PDB IDs")
        seen.add(pdb)
        model, protein = structure_chain(path, pdb, source_chain)
        _, rotation, translation = align_donor(reference, protein, center, {})
        matches = [
            residue
            for part in model
            if part.name == source_chain
            for residue in part
            if residue.name == component and residue.het_flag != "A"
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one {component} in {pdb} chain {source_chain}, found {len(matches)}"
            )
        atoms = [
            coordinate(atom)
            for atom in dominant_atoms(matches[0]).values()
            if atom.element.name not in ("D", "H")
        ]
        if not atoms:
            raise ValueError(
                f"No heavy atoms in {pdb} chain {source_chain} component {component}"
            )
        distances = cdist(site_points, np.asarray(atoms) @ rotation + translation)
        for index in np.where(distances.min(axis=1) < 2.3)[0]:
            support[index].append(pdb)
    return support


def ramp(x, lo, hi):
    t = np.clip((x - lo) / (hi - lo), 0, 1)
    return t * t * (3 - 2 * t)


def replacement(site, overlapping, roles):
    """One occluding atom may replace a fraction of a water's receptor polar contacts."""
    partners = np.array(site["partners"])
    reference = np.array(site["xyz"]) - partners
    reference /= np.linalg.norm(reference, axis=1, keepdims=True)
    delta = overlapping[:, None, :] - partners[None, :, :]
    dist = np.linalg.norm(delta, axis=2)
    cosine = np.sum(delta * reference[None, :, :], axis=2) / np.maximum(dist, 1e-12)
    angle = np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))
    geom = (
        ramp(dist, 2.4, 2.7)
        * (1 - ramp(dist, 3.0, 3.5))
        * (1 - ramp(angle, 20.0, 45.0))
    )
    other = np.asarray(site["partner_roles"])[None, :]
    ligand = roles[:, None]
    complementary = (((ligand & 1) != 0) & ((other & 2) != 0)) | (
        ((ligand & 2) != 0) & ((other & 1) != 0)
    )
    return geom.mean(axis=1).max(), (geom * complementary).mean(axis=1).max()


def score_components(sites, xyz, roles):
    """Return terms in kcal/mol after marginalizing water retention and displacement."""
    if not sites:
        return dict.fromkeys(
            ("exclusion", "contact", "joint", "replace_any", "replace_polar"), 0.0
        )
    loc = np.array([s["xyz"] for s in sites])
    p = np.array([s["p"] for s in sites])
    if not np.all((p > 0) & (p < 1)):
        raise ValueError("Site occupancy must be finite and strictly between 0 and 1")
    dist = cdist(loc, xyz)
    # One water can be excluded once and can earn at most one ligand interaction.
    accessible = ramp(dist.min(axis=1), 2.0, 2.6)
    strength = np.zeros(len(sites))
    any_contact = np.zeros(len(sites))
    polar_contact = np.zeros(len(sites))
    for i, s in enumerate(sites):
        overlap = dist[i] < 2.6
        if np.any(overlap):
            any_contact[i], polar_contact[i] = replacement(
                s, xyz[overlap], roles[overlap]
            )
        candidate = (roles != 0) & (dist[i] > 2.4) & (dist[i] < 3.5)
        if not np.any(candidate):
            continue
        u = (xyz[candidate] - loc[i]) / dist[i, candidate, None]
        v = np.array(s["partners"]) - loc[i]
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        # Tetrahedrally compatible receptor-O-ligand angle, no inferred water hydrogen positions.
        angle = np.rad2deg(np.arccos(np.clip(u @ v.T, -1, 1)))
        angular = (ramp(angle, 65.0, 85.0) * (1 - ramp(angle, 140.0, 160.0))).max(
            axis=1
        )
        r = dist[i, candidate]
        strength[i] = (ramp(r, 2.4, 2.7) * (1 - ramp(r, 3.0, 3.5)) * angular).max()
    interaction = np.exp(-WATER_LIGAND_CONTACT * strength / RT)

    def free_energy(weight):
        return float((-RT * np.log1p(p * (weight - 1))).sum())

    return {
        "exclusion": free_energy(accessible),
        "contact": free_energy(interaction),
        "joint": free_energy(accessible * interaction),
        "replace_any": free_energy(
            accessible * interaction + (1 - accessible) * any_contact
        ),
        "replace_polar": free_energy(
            accessible * interaction + (1 - accessible) * polar_contact
        ),
    }


def box(path):
    settings = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            name, value = line.split("=", 1)
            settings[name.strip()] = float(value)
    return tuple(
        np.asarray([settings[p + "_" + c] for c in "xyz"]) for p in ("center", "size")
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("prepare")
    s.add_argument("rec", help="receptor PDB ID")
    s.add_argument("chain")
    for key in ("pdb", "receptor", "box"):
        s.add_argument(key, type=Path)
    s.add_argument(
        "--donor",
        nargs=3,
        action="append",
        required=True,
        metavar=("PDB_ID", "CHAIN", "PDB_FILE"),
    )
    s.add_argument(
        "--displacer",
        nargs=4,
        action="append",
        default=[],
        metavar=("PDB_ID", "CHAIN", "COMPONENT", "PDB_FILE"),
        help="omit sites overlapped by a ligand in an external structure",
    )
    s = sub.add_parser("score")
    s.add_argument("sites", type=Path)
    s.add_argument("pool", type=Path)
    a = p.parse_args()
    if a.command == "prepare":
        result = build_sites(
            a.pdb,
            a.rec,
            a.chain,
            a.receptor,
            *box(a.box),
            [(i, c, Path(f)) for i, c, f in a.donor],
        )
        if a.displacer:
            sites = result["sites"]
            witnesses = displacement_witnesses(
                a.pdb,
                a.rec,
                a.chain,
                box(a.box)[0],
                sites,
                [(i, c, component, Path(f)) for i, c, component, f in a.displacer],
            )
            result["excluded_sites"] = [
                {**site, "witnesses": w}
                for site, w in zip(sites, witnesses, strict=True)
                if w
            ]
            result["sites"] = [
                site for site, w in zip(sites, witnesses, strict=True) if not w
            ]
    else:
        sites = json.loads(a.sites.read_text())["sites"]
        result = {}
        for model, lines in read_models(a.pool.read_text().splitlines()).items():
            xyz, polar, _ = atom_arrays(lines)
            result[model] = score_components(sites, xyz, polar)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
