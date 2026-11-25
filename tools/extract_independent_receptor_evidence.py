"""Align independent PDB structures and extract water and side-chain evidence."""

import argparse
from collections import defaultdict
from dataclasses import dataclass
import json
from math import atan2, cos, degrees, radians, sin
from pathlib import Path
import re

from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1
import gemmi
import numpy as np


BACKBONE = ("N", "CA", "C", "O")
MIN_SEQUENCE_IDENTITY = 0.90
MIN_SHORTER_CHAIN_COVERAGE = 0.75
MIN_MAPPED_RESIDUES = 25
LOCAL_RADIUS = 12.0
MIN_LOCAL_RESIDUES = 8
MAX_LOCAL_RMSD = 1.25
REGION_PADDING = 1.5
WATER_MIN_OCCUPANCY = 0.5
WATER_CLUSTER_DISTANCE = 1.0
SITE_RADIUS = 6.5
MAX_AUTO_SITES = 3
MIN_DISTINCT_CHI = 35.0
SAME_STATE_CHI = 30.0

# First path follows standard chi atom order, remaining paths cover branches.
# Ring atoms all rotate together around the exocyclic bond.
SIDECHAIN_PATHS = {
    "ARG": ("CA CB CG CD NE CZ NH1", "CZ NH2"),
    "ASN": ("CA CB CG OD1", "CG ND2"),
    "ASP": ("CA CB CG OD1", "CG OD2"),
    "CYS": ("CA CB SG",),
    "GLN": ("CA CB CG CD OE1", "CD NE2"),
    "GLU": ("CA CB CG CD OE1", "CD OE2"),
    "HIS": ("CA CB CG ND1 CE1 NE2 CD2",),
    "ILE": ("CA CB CG1 CD1", "CB CG2"),
    "LEU": ("CA CB CG CD1", "CG CD2"),
    "LYS": ("CA CB CG CD CE NZ",),
    "MET": ("CA CB CG SD CE",),
    "PHE": ("CA CB CG CD1 CE1 CZ CE2 CD2",),
    "SER": ("CA CB OG",),
    "THR": ("CA CB OG1", "CB CG2"),
    "TRP": ("CA CB CG CD1 NE1 CE2 CZ2 CH2 CZ3 CE3 CD2",),
    "TYR": ("CA CB CG CD1 CE1 CZ CE2 CD2", "CZ OH"),
    "VAL": ("CA CB CG1", "CB CG2"),
}
CHI_COUNTS = {
    "ARG": 4, "ASN": 2, "ASP": 2, "CYS": 1, "GLN": 3, "GLU": 3,
    "HIS": 2, "ILE": 2, "LEU": 2, "LYS": 4, "MET": 3, "PHE": 2,
    "SER": 1, "THR": 1, "TRP": 2, "TYR": 2, "VAL": 1,
}
SYMMETRIC_TERMINAL = {"ASP", "GLU", "PHE", "TYR"}


@dataclass(frozen=True)
class Residue:
    name: str
    atoms: dict


@dataclass(frozen=True)
class TargetSite:
    key: tuple
    name: str
    rows: dict
    baseline: tuple
    distance: float


def coordinate(atom):
    return np.array((atom.pos.x, atom.pos.y, atom.pos.z), dtype=float)


def dominant_atoms(residue):
    alternatives = defaultdict(list)
    for atom in residue:
        if atom.altloc not in {"\0", " "}:
            alternatives[atom.altloc].append(atom.occ)
    dominant = min(alternatives, key=lambda key: (-np.mean(alternatives[key]), key)) if alternatives else None
    atoms = {}
    for atom in residue:
        if atom.altloc in {"\0", " ", dominant}:
            prior = atoms.get(atom.name)
            if prior is None or atom.occ > prior.occ:
                atoms[atom.name] = atom
    return atoms


def structure_chain(path, expected_id, chain):
    structure = gemmi.read_structure(str(path))
    actual_id = structure.info["_entry.id"].strip().upper() if "_entry.id" in structure.info else ""
    actual_id = actual_id.removeprefix("PDB_0000")
    if actual_id and actual_id != expected_id:
        raise ValueError(f"structure entry {actual_id} is not manifest PDB {expected_id}")
    if not len(structure):
        raise ValueError("structure lacks a coordinate model")
    model = structure[0]
    residues = {}
    for part in model:
        if part.name != chain:
            continue
        for residue in part:
            if residue.het_flag != "A" or residue.name not in protein_letters_3to1:
                continue
            key = (residue.seqid.num, residue.seqid.icode.strip())
            atoms = {name: coordinate(atom) for name, atom in dominant_atoms(residue).items()
                     if atom.element.name not in {"H", "D"}}
            if "CA" in atoms and key not in residues:
                residues[key] = Residue(residue.name, atoms)
    if not residues:
        raise ValueError(f"no canonical protein chain {chain}")
    return model, residues


def prepared_receptor(path, chain, raw):
    lines = path.read_text().splitlines(keepends=True)
    groups = defaultdict(lambda: defaultdict(list))
    for index, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")) or line[21].strip() != chain:
            continue
        residue = line[17:20].strip()
        if residue not in protein_letters_3to1:
            continue
        key = (int(line[22:26]), line[26].strip(), residue)
        name = line[12:16].strip()
        xyz = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError(f"nonfinite prepared coordinate {chain}:{key}:{name}")
        groups[key][name].append((index, xyz, residue, line[76:].strip()))

    ca_by_name = defaultdict(list)
    for key, atoms in groups.items():
        ca_by_name[key[2]].extend((key, row) for row in atoms.get("CA", ()))
    by_coordinate = {}
    for key, residue in raw.items():
        candidates = sorted(
            (float(np.linalg.norm(ca[1] - residue.atoms["CA"])), ca[0], prepared_key, ca)
            for prepared_key, ca in ca_by_name[residue.name]
        )
        if candidates and candidates[0][0] <= 0.75 and (len(candidates) == 1 or candidates[1][0] > 0.75):
            by_coordinate[key] = candidates[0][2:]
    uses_coordinates = len(by_coordinate) >= MIN_LOCAL_RESIDUES
    if uses_coordinates:
        matches = by_coordinate
    else:
        # For a translated or rotated preparation, unambiguous same-identity
        # author labels can establish the backbone alignment independently.
        matches = {}
        for key, residue in raw.items():
            prepared_key = (*key, residue.name)
            ca = groups.get(prepared_key, {}).get("CA", ())
            if len(ca) == 1:
                matches[key] = (prepared_key, ca[0])

    prepared = {}
    mapping = []
    for key, (prepared_key, ca) in matches.items():
        source = raw[key]
        group = groups[prepared_key]
        rows = {"CA": ca}
        hydrogens = []
        for name, candidates in group.items():
            if name == "CA":
                continue
            if all(row[3] in {"H", "HD", "HS"} for row in candidates):
                hydrogens.extend(candidates)
                continue
            if len(candidates) == 1 and len(group["CA"]) == 1:
                rows[name] = candidates[0]
                continue
            if not uses_coordinates or name not in source.atoms:
                continue
            near = [row for row in candidates if np.linalg.norm(row[1] - source.atoms[name]) <= 0.75]
            if len(near) == 1:
                rows[name] = near[0]
        for row in hydrogens:
            # Repeated H names are normal. Bind each physical row to the
            # nearest heavy atom of its co-numbered prepared residue(s).
            _, host = min((np.linalg.norm(row[1] - heavy[1]), heavy[0])
                          for candidates in group.values() for heavy in candidates
                          if heavy[3] not in {"H", "HD", "HS"})
            if host in {heavy[0] for heavy in rows.values()}:
                rows[f"H@{row[0]}"] = row
        prepared[key] = rows
        mapping.append({"query": key, "prepared_chain": chain, "prepared_residue": prepared_key[:2],
                        "prepared_ca_row": ca[0] + 1})
    return lines, prepared, {"method": "same-residue C-alpha coordinate match" if uses_coordinates else
                             "unambiguous author labels followed by rigid registration", "residue_rows": mapping}


def fit_coordinates(source, destination):
    centered_src = source - source.mean(axis=0)
    centered_dst = destination - destination.mean(axis=0)
    u, singular, vt = np.linalg.svd(centered_src.T @ centered_dst)
    if singular[1] <= 1e-8:
        raise ValueError("degenerate backbone fit")
    u[:, -1] *= np.linalg.det(u @ vt)
    rotation = u @ vt
    translation = destination.mean(axis=0) - source.mean(axis=0) @ rotation
    error = np.linalg.norm(source @ rotation + translation - destination, axis=1)
    return rotation, translation, error


def verify_query_frame(raw, prepared):
    source, destination = [], []
    for key, residue in raw.items():
        for name in BACKBONE:
            row = prepared.get(key, {}).get(name)
            if row is not None and row[2] == residue.name and name in residue.atoms:
                source.append(residue.atoms[name])
                destination.append(row[1])
    if len(source) < 4 * MIN_LOCAL_RESIDUES:
        raise ValueError("too few matching query and prepared backbone atoms")
    source, destination = np.array(source), np.array(destination)
    direct_error = float(np.sqrt(np.mean(np.square(np.linalg.norm(source - destination, axis=1)))))
    if direct_error <= 0.01:
        rotation, translation, rmsd = np.eye(3), np.zeros(3), direct_error
        registered = raw
    else:
        rotation, translation, error = fit_coordinates(source, destination)
        rmsd = float(np.sqrt(np.mean(np.square(error))))
        if rmsd > 0.25:
            raise ValueError(f"query/prepared backbone did not match after frame registration (RMSD {rmsd:.3f} Angstrom)")
        registered = {key: Residue(residue.name, {name: xyz @ rotation + translation
                                                for name, xyz in residue.atoms.items()})
                      for key, residue in raw.items()}
    return registered, {"matched_atoms": len(source), "direct_rmsd_angstrom": direct_error,
                        "rmsd_angstrom": rmsd, "rotation": rotation.tolist(),
                        "translation": translation.tolist()}


def vector3(value, label):
    value = np.asarray(value, dtype=float)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise ValueError(f"{label} requires three finite values")
    return value


def extraction_region(query, protein, root):
    anchors = query.get("protein_anchors")
    box = query.get("box")
    if anchors and box is not None:
        raise ValueError("choose protein_anchors or an independent box")
    if anchors:
        points = []
        for anchor in anchors:
            key = (int(anchor["resseq"]), anchor.get("icode", "").strip())
            name = anchor.get("atom", "CA")
            if key not in protein or name not in protein[key].atoms:
                raise ValueError(f"missing protein anchor {key}:{name}")
            points.append(protein[key].atoms[name])
        center = np.mean(np.array(points), axis=0)
        size = vector3(query.get("box_size", [16, 16, 16]), "box_size")
        basis = {"kind": "query_protein_anchors", "anchors": anchors, "description": query.get("box_basis")}
    elif box is not None:
        description = query.get("box_basis")
        if not isinstance(description, str) or not description.strip():
            raise ValueError("box_basis must describe the independent provenance of the extraction box")
        if isinstance(box, str):
            config = {}
            for line in (root / box).read_text().splitlines():
                name, sep, value = line.partition("#")[0].partition("=")
                if sep:
                    config[name.strip()] = value.strip()
            box = {key: [config[f"{key}_{axis}"] for axis in "xyz"] for key in ("center", "size")}
        if isinstance(box, dict):
            center, size = vector3(box["center"], "center"), vector3(box["size"], "size")
        else:
            if len(box) != 6:
                raise ValueError("box array requires cx, cy, cz, sx, sy, sz")
            center, size = vector3(box[:3], "center"), vector3(box[3:], "size")
        basis = {"kind": "independent_box", "description": description, "input": query["box"]}
    else:
        raise ValueError("independent protein anchors or an independently defined box are required")
    if np.any(size <= 0):
        raise ValueError("extraction box has a nonpositive axis")
    return center, size, basis


def align_donor(query, donor, center, report):
    query_keys, donor_keys = list(query), list(donor)
    sequence_query, sequence_donor = (
        "".join(protein_letters_3to1[record.name] for record in residues.values())
        for residues in (query, donor)
    )
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                              open_gap_score=-8, extend_gap_score=-0.5, end_gap_score=0)
    pair = aligner.align(sequence_query, sequence_donor)[0]
    mapping = {}
    count = 0
    for (a, b), (c, d) in zip(*pair.aligned):
        for qi, di in zip(range(a, b), range(c, d)):
            count += 1
            if sequence_query[qi] == sequence_donor[di]:
                mapping[query_keys[qi]] = donor_keys[di]
    identity = len(mapping) / count if count else 0
    coverage = count / min(len(query), len(donor))
    report.update(sequence_identity=identity, shorter_chain_coverage=coverage, mapped_identical_residues=len(mapping))
    if identity < MIN_SEQUENCE_IDENTITY or coverage < MIN_SHORTER_CHAIN_COVERAGE or len(mapping) < MIN_MAPPED_RESIDUES:
        raise ValueError("donor sequence identity, chain coverage, or residue mapping below acceptance threshold")

    local = [key for key, mapped in mapping.items()
             if np.linalg.norm(query[key].atoms["CA"] - center) <= LOCAL_RADIUS
             and all(name in query[key].atoms and name in donor[mapped].atoms for name in BACKBONE)]
    report["local_residues"] = len(local)
    if len(local) < MIN_LOCAL_RESIDUES:
        raise ValueError("too few matching local backbone residues")
    q = np.array([query[key].atoms[name] for key in local for name in BACKBONE])
    d = np.array([donor[mapping[key]].atoms[name] for key in local for name in BACKBONE])
    rotation, translation, error = fit_coordinates(d, q)
    rmsd = float(np.sqrt(np.mean(np.square(error))))
    report.update(local_backbone_rmsd_angstrom=rmsd, max_fit_atom_error_angstrom=float(max(error)),
                  rotation=rotation.tolist(), translation=translation.tolist(),
                  local_residue_mapping=[{"query": key, "donor": mapping[key]} for key in local])
    if rmsd > MAX_LOCAL_RMSD:
        raise ValueError("local backbone RMSD above acceptance threshold")
    return mapping, rotation, translation


def aligned_waters(model, pdb_id, center, size, rotation, translation):
    for chain in model:
        for residue in chain:
            if not residue.is_water():
                continue
            oxygens = [a for a in dominant_atoms(residue).values()
                       if a.element.name == "O" and a.occ >= WATER_MIN_OCCUPANCY]
            if not oxygens:
                continue
            oxygen = max(oxygens, key=lambda atom: atom.occ)
            xyz = coordinate(oxygen) @ rotation + translation
            if np.all(np.abs(xyz - center) <= size / 2 + REGION_PADDING):
                yield {"pdb_id": pdb_id, "chain": chain.name, "source_resseq": residue.seqid.num,
                       "source_icode": residue.seqid.icode.strip(), "atom": oxygen.name,
                       "occupancy": oxygen.occ, "source_xyz": coordinate(oxygen).tolist(), "xyz": xyz.tolist()}


def water_sites(waters):
    groups = []
    for water in sorted(waters, key=lambda w: (w["pdb_id"], -w["occupancy"], w["chain"],
                                               w["source_resseq"], w["source_icode"])):
        point = np.asarray(water["xyz"])
        choices = [(float(np.linalg.norm(point - group["seed"])), i) for i, group in enumerate(groups)]
        nearest = min(choices, default=(float("inf"), 0))
        if nearest[0] > WATER_CLUSTER_DISTANCE:
            groups.append({"seed": point, "hits": {}, "duplicates": []})
            group = groups[-1]
        else:
            group = groups[nearest[1]]
        if water["pdb_id"] in group["hits"]:
            group["duplicates"].append(water)
        else:
            group["hits"][water["pdb_id"]] = water
    result = []
    for index, group in enumerate(groups, 1):
        hits = list(group["hits"].values())
        xyz = np.mean([hit["xyz"] for hit in hits], axis=0)
        result.append({"id": index, "xyz": xyz.tolist(), "independent_pdb_count": len(hits),
                       "hits": hits, "same_pdb_repeats": group["duplicates"]})
    return result


def chis(name, atoms):
    path = ["N", *SIDECHAIN_PATHS[name][0].split()]
    windows = [path[i:i + 4] for i in range(CHI_COUNTS[name])]
    angles = []
    for window in windows:
        if not all(atom in atoms for atom in window):
            return None
        a, b, c, d = (atoms[atom] for atom in window)
        norm = np.linalg.norm(c - b)
        if norm < 1e-8:
            return None
        axis = (c - b) / norm
        v, w = (offset - np.dot(offset, axis) * axis for offset in (a - b, d - c))
        if min(np.linalg.norm(v), np.linalg.norm(w)) < 1e-8:
            return None
        angles.append(degrees(atan2(np.dot(np.cross(axis, v), w), np.dot(v, w))))
    return tuple(angles)


def angle_difference(name, index, target, current):
    symmetric = name in SYMMETRIC_TERMINAL and index == CHI_COUNTS[name] - 1
    period = 180 if symmetric else 360
    return (target - current + period / 2) % period - period / 2


def prepared_heavy_atoms(name):
    return set(BACKBONE) | {atom for path in SIDECHAIN_PATHS[name] for atom in path.split()}


def receptor_sites(target, query, prepared, center):
    requested = target.get("sidechain_sites")
    if requested is None:
        keys = list(query)
    else:
        keys = [(int(site["resseq"]), site.get("icode", "").strip()) for site in requested]
        if len(keys) != len(set(keys)):
            raise ValueError("duplicate configured side-chain residues")

    sites, rejected = [], []
    for key in keys:
        residue = query.get(key)
        rows = prepared.get(key, {})
        name = residue.name if residue else None
        if name not in SIDECHAIN_PATHS:
            if requested is not None:
                rejected.append({"key": key, "reason": "unsupported or absent reference residue"})
            continue
        heavy = prepared_heavy_atoms(name)
        extra = set(rows) - heavy - {"OXT"}
        if not heavy.issubset(rows) or any(rows[a][3] not in {"H", "HD", "HS"} for a in extra):
            if requested is not None:
                rejected.append({"key": key, "reason": "incomplete or unexpected prepared atom topology"})
            continue
        if any(row[2] != name for row in rows.values()):
            if requested is not None:
                rejected.append({"key": key, "reason": "query and prepared residue identity differ"})
            continue
        atoms = {atom: row[1] for atom, row in rows.items()}
        baseline = chis(name, atoms)
        if baseline is None:
            if requested is not None:
                rejected.append({"key": key, "reason": "query chi is undefined"})
            continue
        distance = min(float(np.linalg.norm(atoms[atom] - center)) for atom in heavy - set(BACKBONE))
        if requested is not None or distance <= SITE_RADIUS:
            sites.append(TargetSite(key, name, rows, baseline, distance))
    if requested is None:
        sites = sorted(sites, key=lambda site: (site.distance, site.key))[:MAX_AUTO_SITES]
    return sites, rejected


def max_state_difference(sites, first, second):
    return max((abs(angle_difference(site.name, i, angle_a, angle_b))
                for site, a, b in zip(sites, first, second)
                for i, (angle_a, angle_b) in enumerate(zip(a, b))), default=0.0)


def sidechain_states(sites, relationships):
    baseline = tuple(site.baseline for site in sites)
    groups = [{"chi": baseline, "hits": []}]
    missing = []
    for source, mapping, donor in relationships:
        angles, provenance = [], []
        for site in sites:
            mapped = mapping.get(site.key)
            sidechain = donor.get(mapped)
            values = chis(site.name, sidechain.atoms) if sidechain and sidechain.name == site.name else None
            if values is None:
                missing.append({"pdb_id": source, "missing_query_residue": site.key})
                break
            angles.append(values)
            provenance.append({"query": site.key, "donor": mapped, "chi_degrees": values})
        if len(angles) != len(sites):
            continue
        hit = {"pdb_id": source, "residues": provenance}
        if max_state_difference(sites, angles, baseline) <= MIN_DISTINCT_CHI:
            groups[0]["hits"].append(hit)
            continue
        distances = [(max_state_difference(sites, angles, group["chi"]), index)
                     for index, group in enumerate(groups) if index > 0]
        distance, index = min(distances, default=(float("inf"), 0))
        if distance > SAME_STATE_CHI:
            groups.append({"chi": tuple(angles), "hits": []})
            index = len(groups) - 1
        groups[index]["hits"].append(hit)
    return groups, missing


def move_sidechain(site, target):
    paths = [path.split() for path in SIDECHAIN_PATHS[site.name]]
    chain = ["N", *paths[0]]
    adjacency = defaultdict(set)
    for path in paths:
        for atom, other in zip(path, path[1:]):
            adjacency[atom].add(other)
            adjacency[other].add(atom)
    atoms = {name: xyz.copy() for name, (_, xyz, _, _) in site.rows.items()}
    heavy = prepared_heavy_atoms(site.name)
    for name, (_, xyz, _, kind) in site.rows.items():
        if kind in {"H", "HD", "HS"}:
            _, host = min((np.linalg.norm(xyz - atoms[atom]), atom) for atom in heavy)
            adjacency[host].add(name)
    moved = set()
    for index, desired in enumerate(target):
        current = chis(site.name, atoms)[index]
        change = angle_difference(site.name, index, desired, current)
        if abs(change) < 1e-10:
            continue
        fixed, pivot = chain[index + 1:index + 3]
        stack, component = [pivot], {fixed}
        while stack:
            name = stack.pop()
            if name not in component:
                component.add(name)
                stack.extend(adjacency[name])
        component.difference_update({fixed, pivot})
        origin = atoms[fixed]
        axis = atoms[pivot] - origin
        axis /= np.linalg.norm(axis)
        angle = radians(change)
        for name in component:
            offset = atoms[name] - origin
            atoms[name] = origin + offset * cos(angle) + np.cross(axis, offset) * sin(angle) + axis * np.dot(axis, offset) * (1 - cos(angle))
        moved.update(component)
    return {site.rows[name][0]: atoms[name] for name in moved}


def write_sidechains(directory, lines, sites, groups, missing, rejected):
    total = sum(1 + len(group["hits"]) for group in groups)
    alternate = []
    provenance = {"selected_query_sites": [{"key": site.key, "name": site.name,
                                             "native_chi_degrees": site.baseline, "pocket_center_distance": site.distance}
                                            for site in sites],
                  "unavailable_query_sites": rejected, "incomplete_donor_observations": missing,
                  "baseline_observations": groups[0]["hits"], "population_prior_per_state": 1}
    for index, group in enumerate(groups[1:], 1):
        changed = {}
        included = []
        for site, angles in zip(sites, group["chi"]):
            if max_state_difference((site,), (angles,), (site.baseline,)) > MIN_DISTINCT_CHI:
                changed.update(move_sidechain(site, angles))
                included.append(site.key)
        if not changed:
            raise ValueError("distinct side-chain state produced no alternate coordinates")
        changed_lines = lines.copy()
        for row, xyz in changed.items():
            line = lines[row]
            fields = "".join(f"{value:8.3f}" for value in xyz)
            if len(fields) != 24 or not np.isfinite(xyz).all():
                raise ValueError("alternate coordinate outside PDBQT representable range")
            changed_lines[row] = line[:30] + fields + line[54:]
        identity = f"state-{index:02d}"
        filename = identity + ".pdbqt"
        (directory / filename).write_text("".join(changed_lines))
        alternate.append({"id": identity, "receptor": filename, "population": (1 + len(group["hits"])) / total,
                          "source_pdbs": [hit["pdb_id"] for hit in group["hits"]],
                          "moved_query_sites": included, "representative_chi_degrees": group["chi"],
                          "observations": group["hits"], "moved_atom_count": len(changed)})
    result = {"baseline_population": (1 + len(groups[0]["hits"])) / total,
              "alternates": alternate, "provenance": provenance}
    (directory / "sidechains.json").write_text(json.dumps(result, indent=2) + "\n")
    return len(alternate)


def extract_target(target, root, directory):
    directory.mkdir()
    query = target["query"]
    pdb_id = query["pdb_id"].upper()
    query_path = (root / query["structure"]).resolve(strict=True)
    receptor_path = (root / query["receptor"]).resolve(strict=True)
    _, residues = structure_chain(query_path, pdb_id, query["chain"])
    lines, prepared, prepared_mapping = prepared_receptor(receptor_path, query.get("receptor_chain", query["chain"]), residues)
    residues, frame = verify_query_frame(residues, prepared)
    frame["prepared_mapping"] = prepared_mapping
    center, size, basis = extraction_region(query, residues, root)
    sites, rejected = receptor_sites(target, residues, prepared, center)
    recorded, water_hits, observations, seen = [], [], [], set()

    for source in sorted(target["donors"], key=lambda donor: donor["pdb_id"].upper()):
        donor_id = source["pdb_id"].upper()
        report = {**source, "pdb_id": donor_id, "status": "rejected"}
        recorded.append(report)
        try:
            if donor_id == pdb_id:
                raise ValueError("query PDB cannot be its own donor")
            if donor_id in seen:
                raise ValueError("duplicate donor PDB")
            if source.get("accession") != target["accession"] or not source.get("condition"):
                raise ValueError("accession mismatch or missing independent donor condition")
            donor_path = (root / source["structure"]).resolve(strict=True)
            if donor_path.samefile(query_path) or donor_path.samefile(receptor_path):
                raise ValueError("query structure or prepared receptor cannot be its own donor")
            donor_model, donor = structure_chain(donor_path, donor_id, source["chain"])
            mapping, rotation, translation = align_donor(residues, donor, center, report)
            hits = list(aligned_waters(donor_model, donor_id, center, size, rotation, translation))
            water_hits.extend(hits)
            observations.append((donor_id, mapping, donor))
            report.update(status="accepted", candidate_water_oxygens=len(hits), structure=str(donor_path))
            seen.add(donor_id)
        except (OSError, RuntimeError, ValueError) as exc:
            report["reason"] = str(exc)

    waters = water_sites(water_hits)
    groups, missing = sidechain_states(sites, observations)
    alternates = write_sidechains(directory, lines, sites, groups, missing, rejected)
    (directory / "waters.xyz").write_text("".join(
        " ".join(f"{value:.6f}" for value in water["xyz"]) + "\n" for water in waters
    ))
    provenance = {"id": target["id"], "accession": target["accession"], "query": query,
                  "query_frame": frame, "box": {"center": center.tolist(), "size": size.tolist(), "basis": basis},
                  "donors": recorded, "water_sites": waters, "alternate_receptors": alternates}
    (directory / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return {"id": target["id"], "accepted_donors": len(seen), "water_sites": len(waters),
            "alternate_receptors": alternates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="independent query protein region and donor sources")
    parser.add_argument("--output-dir", type=Path, required=True, help="new result directory")
    args = parser.parse_args()
    manifest = args.manifest.resolve(strict=True)
    targets = json.loads(manifest.read_text())["targets"]
    ids = [target["id"] for target in targets]
    if not ids or len(set(ids)) != len(ids) or not all(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", key) for key in ids):
        parser.error("target IDs must be nonempty, unique, and path-safe")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    summaries = []
    for target in targets:
        directory = args.output_dir / target["id"]
        try:
            summaries.append(extract_target(target, manifest.parent, directory))
        except (OSError, RuntimeError, ValueError) as exc:
            # Preserve the unavailable query as an explicit baseline-only row.
            if any(directory.iterdir()):
                raise
            failure = {"id": target["id"], "status": "unavailable_query", "reason": str(exc)}
            (directory / "waters.xyz").write_text("")
            (directory / "sidechains.json").write_text(json.dumps(
                {"baseline_population": 1.0, "alternates": [], "provenance": failure}, indent=2) + "\n")
            (directory / "provenance.json").write_text(json.dumps(failure, indent=2) + "\n")
            summaries.append({**failure, "accepted_donors": 0, "water_sites": 0, "alternate_receptors": 0})
    text = json.dumps({"manifest": str(manifest), "targets": summaries}, indent=2, allow_nan=False) + "\n"
    (args.output_dir / "summary.json").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
