#!/usr/bin/env python3
"""Select a pose for a declared monodentate donor at a three-coordinate protein zinc site."""

import argparse
from dataclasses import dataclass
from itertools import combinations
import json
import math
from pathlib import Path
from typing import Tuple


COORDINATION_DISTANCES = {"N": (2.00, 2.50), "O": (2.10, 2.50), "S": (2.25, 2.80)}
# AutoDock types are case-sensitive. Se is an alias for S in Vina's scoring,
# but this geometry is defined for actual nitrogen, oxygen and sulfur atoms.
COORDINATOR_TYPES = {"N": "N", "NA": "N", "O": "O", "OA": "O", "S": "S", "SA": "S"}
HALOGEN_TYPES = {"F", "Cl", "Br", "I", "At"}
PROTEIN_RESIDUES = set("ALA ARG ASN ASP ASH CYS CYM CYX GLN GLU GLH GLY HIS HID HIE HIP HSD HSE HSP ILE LEU LYS LYN MET PHE PRO SER THR TRP TYR VAL".split())
MIN_PLANE_OFFSET = 0.20
TETRAHEDRAL_ANGLE = math.degrees(math.acos(-1 / 3))
SITE_ANGLE_TOLERANCE = 25.0
DONOR_ANGLE_LIMIT = 30.0
SHORT_DISTANCE_TOLERANCE = 0.35

NATIVE_RESULT = "REMARK VINA RESULT:"
PRO_RESULT = "REMARK VINA PRO RESULT:"
PRO_NATIVE = "REMARK VINA NATIVE:"
PRO_SCALE = "REMARK VINA PRO GAUSS2 SCALE:"
SCORE_HEADERS = (NATIVE_RESULT, PRO_RESULT, PRO_NATIVE, PRO_SCALE)

Point = Tuple[float, float, float]


@dataclass(frozen=True)
class Atom:
    serial: int
    residue: str
    atom_type: str
    xyz: Point


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    x, y, z = a
    u, v, w = b
    return (y * w - z * v, z * u - x * w, x * v - y * u)


def angle_deg(a, b):
    cosine = dot(a, b) / (math.hypot(*a) * math.hypot(*b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def parse_atoms(lines):
    atoms = {}
    for line in lines:
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        serial = int(line[6:11])
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        atom_type = line[77:].strip()
        if serial in atoms or not atom_type or not all(map(math.isfinite, xyz)):
            raise ValueError(f"duplicate serial or invalid PDBQT atom {serial}")
        atoms[serial] = Atom(serial, line[17:20].strip().upper(), atom_type, xyz)
    return atoms


def nearby_coordinators(atoms, zinc_xyz):
    for atom in atoms.values():
        element = COORDINATOR_TYPES.get(atom.atom_type)
        if element is not None:
            distance = math.dist(atom.xyz, zinc_xyz)
            if distance <= COORDINATION_DISTANCES[element][1]:
                yield atom, element, distance


def assess_site(receptor, zinc_serial):
    zinc = receptor.get(zinc_serial)
    site = {"zinc_serial": zinc_serial, "zinc_xyz": None, "coordinators": [], "applicable": False}
    if zinc is None or zinc.atom_type != "Zn":
        return dict(site, reason="requested_zinc_not_found")
    shell = list(nearby_coordinators(receptor, zinc.xyz))
    site["zinc_xyz"] = zinc.xyz
    site["coordinators"] = [
        {"serial": atom.serial, "residue": atom.residue, "element": element, "distance_angstrom": distance}
        for atom, element, distance in shell
    ]
    largest_shell = max(upper for _, upper in COORDINATION_DISTANCES.values())
    nearby_halogens = [atom.serial for atom in receptor.values() if atom.atom_type in HALOGEN_TYPES
                       and math.dist(atom.xyz, zinc.xyz) <= largest_shell]
    if nearby_halogens:
        return dict(site, reason="nearby_receptor_halogen", nearby_halogen_serials=nearby_halogens)
    if len(shell) != 3:
        return dict(site, reason="expected_three_receptor_coordinators")
    if any(atom.residue not in PROTEIN_RESIDUES for atom, _, _ in shell):
        return dict(site, reason="non_protein_coordinator")
    if any(distance < COORDINATION_DISTANCES[element][0] - SHORT_DISTANCE_TOLERANCE
           for _, element, distance in shell):
        return dict(site, reason="receptor_coordinator_too_close")

    a, b, c = [atom.xyz for atom, _, _ in shell]
    normal = cross(subtract(b, a), subtract(c, a))
    magnitude = math.hypot(*normal)
    if magnitude < 1e-8:
        return dict(site, reason="degenerate_coordinator_plane")
    direction = tuple(component / magnitude for component in normal)
    plane_offset = dot(subtract(zinc.xyz, a), direction)
    site["plane_offset_angstrom"] = abs(plane_offset)
    if abs(plane_offset) < MIN_PLANE_OFFSET:
        return dict(site, reason="zinc_too_close_to_coordinator_plane")
    if plane_offset < 0:
        direction = tuple(-component for component in direction)

    occupied = [subtract(atom.xyz, zinc.xyz) for atom, _, _ in shell]
    occupied_angles = [angle_deg(first, second) for first, second in combinations(occupied, 2)]
    vacancy_angles = [angle_deg(direction, ray) for ray in occupied]
    site.update(vacant_direction=direction, occupied_angles_deg=occupied_angles, vacancy_angles_deg=vacancy_angles)
    if any(abs(angle - TETRAHEDRAL_ANGLE) > SITE_ANGLE_TOLERANCE for angle in occupied_angles + vacancy_angles):
        return dict(site, reason="outside_tetrahedral_site_tolerance")
    return dict(site, applicable=True, reason=None)


def read_models(lines):
    models = {}
    active = None
    for line in lines:
        tag = line[:6].strip()
        if tag == "MODEL":
            number = int(line.split()[1])
            if active is not None or number in models or number <= 0:
                raise ValueError(f"misplaced or repeated MODEL {number}")
            active = number
            models[number] = []
        elif tag == "ENDMDL" and active is not None:
            active = None
        elif active is not None:
            if line.startswith("BEGIN_RES"):
                raise ValueError("flexible-receptor MODEL output is unsupported")
            models[active].append(line)
        elif line.strip():
            raise ValueError("expected scored PDBQT MODEL blocks")
    if active is not None:
        raise ValueError(f"MODEL {active} is missing ENDMDL")
    return models


def read_score_header(model, lines):
    def fields(index, prefix, length):
        if index >= len(lines) or not lines[index].startswith(prefix):
            raise ValueError(f"MODEL {model} missing {prefix} in its opening score header")
        tokens = lines[index][len(prefix):].split()
        if len(tokens) != length or not all(math.isfinite(float(token)) for token in tokens):
            raise ValueError(f"MODEL {model} has an invalid {prefix} score header")
        return tokens

    if lines and lines[0].startswith(PRO_RESULT):
        pro = float(fields(0, PRO_RESULT, 3)[0])
        native, rank = fields(1, PRO_NATIVE, 2)
        native, rank = float(native), int(rank)
        scale = float(fields(2, PRO_SCALE, 1)[0])
        if rank <= 0 or not 0 <= scale < 1:
            raise ValueError(f"MODEL {model} has an invalid Pro rank or scale")
        end = 3
    else:
        native = float(fields(0, NATIVE_RESULT, 3)[0])
        rank, pro, scale, end = model, None, None, 1
    return {"model": model, "native_rank": rank, "native_score": native, "pro_score": pro, "gauss2_scale": scale,
            "ignored_trailing_score_headers": sum(line.startswith(SCORE_HEADERS) for line in lines[end:])}


def assess_pose(atoms, donor_serial, site):
    donor = atoms.get(donor_serial)
    element = COORDINATOR_TYPES.get(donor.atom_type) if donor else None
    zinc = site["zinc_xyz"]
    distance = math.dist(donor.xyz, zinc) if donor and zinc is not None else None
    result = {"donor_element": element, "zinc_distance_angstrom": distance, "direction_error_deg": None,
              "virtual_point_error_angstrom": None, "other_nearby_coordinator_serials": [], "passes_geometry": False}
    if not site["applicable"]:
        return dict(result, reason="site_inapplicable")
    if element is None:
        return dict(result, reason="donor_missing_or_not_n_o_s")
    ideal, upper = COORDINATION_DISTANCES[element]
    vacancy = site["vacant_direction"]
    virtual_point = tuple(x + ideal * d for x, d in zip(zinc, vacancy))
    deviation = angle_deg(subtract(donor.xyz, zinc), vacancy) if distance > 0 else None
    competitors = [atom.serial for atom, _, _ in nearby_coordinators(atoms, zinc) if atom.serial != donor_serial]
    result.update(direction_error_deg=deviation, virtual_point_error_angstrom=math.dist(donor.xyz, virtual_point),
                  other_nearby_coordinator_serials=competitors)
    if not ideal - SHORT_DISTANCE_TOLERANCE <= distance <= upper:
        return dict(result, reason="donor_distance_outside_window")
    if deviation > DONOR_ANGLE_LIMIT:
        return dict(result, reason="donor_outside_open_direction")
    if competitors:
        return dict(result, reason="other_ligand_coordinators_in_shell")
    return dict(result, passes_geometry=True, reason=None)


def score_sulfonamide_pose(pose, site, mode):
    score = {"scoring_eligible": False, "scoring_delta_angstrom": None,
             "geometry_reward": None, "ranking_score": None}
    if not site["applicable"]:
        return dict(score, scoring_reason="site_inapplicable")
    if pose["donor_element"] != "N":
        return dict(score, scoring_reason="donor_missing_or_not_nitrogen")
    nitrogen_ideal, _ = COORDINATION_DISTANCES["N"]
    distance = pose["zinc_distance_angstrom"]
    if distance < nitrogen_ideal - SHORT_DISTANCE_TOLERANCE:
        return dict(score, scoring_reason="donor_distance_below_minimum")
    if pose["other_nearby_coordinator_serials"]:
        return dict(score, scoring_reason="other_ligand_coordinators_in_shell")

    delta = (pose["virtual_point_error_angstrom"] if mode == "continuous"
             else abs(distance - nitrogen_ideal))
    reward = 3.0 * math.exp(-0.5 * delta * delta)
    return dict(score, scoring_eligible=True, scoring_reason=None, scoring_delta_angstrom=delta,
                geometry_reward=reward, ranking_score=pose["native_score"] - reward)


def select_pose(receptor, models, zinc_serial, donor_serial, mode="hard"):
    site = assess_site(receptor, zinc_serial)
    poses = []
    for model, lines in models.items():
        score = read_score_header(model, lines)
        poses.append(dict(score, **assess_pose(parse_atoms(lines), donor_serial, site)))
    if len({pose["gauss2_scale"] for pose in poses}) > 1 or len({pose["native_rank"] for pose in poses}) != len(poses):
        raise ValueError("pose pool mixes score modes or repeats native ranks")
    if mode == "hard":
        eligible = [pose for pose in poses if pose["passes_geometry"]]
        score_field = "native_score"
    else:
        for pose in poses:
            pose.update(score_sulfonamide_pose(pose, site, mode))
        eligible = [pose for pose in poses if pose["scoring_eligible"]]
        score_field = "ranking_score"
    best = min(eligible, key=lambda p: (p[score_field], p["native_rank"]), default=None)
    reason = None
    if not poses:
        reason = "no_candidates"
    elif not site["applicable"]:
        reason = "no_applicable_site"
    elif best is None:
        reason = "no_geometry_pass" if mode == "hard" else "no_scoring_eligible_pose"
    result = {"status": "selected" if best else "abstain", "reason": reason, "site": site,
              "donor_serial": donor_serial, "selected_model": best["model"] if best else None, "poses": poses}
    if mode != "hard":
        result["mode"] = mode
    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="Use a receptor site with no coordinated species beyond the three protein N/O/S atoms. "
               "AutoDock atom types are case-sensitive (Na is sodium, NA is nitrogen). "
               "A receptor F, Cl, Br, I, or At within 2.8 angstrom of zinc triggers abstention.",
    )
    parser.add_argument("--receptor", required=True, type=Path, help="rigid receptor PDBQT used for docking")
    parser.add_argument("--poses", required=True, type=Path, help="one ligand's scored docking PDBQT")
    parser.add_argument("--zinc-serial", required=True, type=int, help="receptor zinc atom serial")
    parser.add_argument("--donor-serial", required=True, type=int, help="declared monodentate ligand atom serial in the poses")
    parser.add_argument("--mode", choices=("hard", "continuous", "radial"), default="hard",
                        help="pose ranking rule (default: hard); choose continuous or radial only with prior "
                             "evidence for direct monodentate sulfonamide nitrogen coordination")
    args = parser.parse_args()
    if args.zinc_serial <= 0 or args.donor_serial <= 0:
        parser.error("atom serials must be positive")
    try:
        with args.receptor.open() as receptor_file, args.poses.open() as pose_file:
            receptor = parse_atoms(receptor_file)
            models = read_models(line.rstrip("\r\n") for line in pose_file)
        result = select_pose(receptor, models, args.zinc_serial, args.donor_serial, args.mode)
    except (OSError, ValueError, IndexError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
