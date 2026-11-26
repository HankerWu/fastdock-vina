"""Identify movable carboxamides or generate a distinguishable alternate planar rotamer."""

import argparse
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from make_scaffold_anchors import query_graph
from rdkit import Chem
from score_physical_pool import mapped_ligand_bonds


def _carboxamides(graph):
    """Yield non-ring carboxamide C–N bonds with at least one carbon on the nitrogen."""
    for carbon, oxygen, nitrogen in graph.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=[OX1])-[NX3+0]")
    ):
        if graph.GetBondBetweenAtoms(carbon, nitrogen).IsInRing():
            continue
        acyl = [
            a
            for a in graph.GetAtomWithIdx(carbon).GetNeighbors()
            if a.GetIdx() not in (nitrogen, oxygen)
        ]
        groups = [
            a
            for a in graph.GetAtomWithIdx(nitrogen).GetNeighbors()
            if a.GetIdx() != carbon
        ]
        if any(a.GetAtomicNum() != 6 for a in acyl + groups) or not groups:
            continue
        yield oxygen, carbon, nitrogen, groups


def amides(text):
    """Return the original distinguishable rotamer sites in PDBQT serial order."""
    graph, serials = query_graph(text)
    ranks = Chem.CanonicalRankAtoms(graph, breakTies=False)
    sites = []
    for oxygen, carbon, nitrogen, groups in _carboxamides(graph):
        if len({ranks[a.GetIdx()] for a in groups}) != len(groups):
            continue
        tip = max(groups, key=lambda a: ranks[a.GetIdx()]).GetIdx()
        sites.append(tuple(serials[i] for i in (oxygen, carbon, nitrogen, tip)))
    return sites


def active_amide_sites(text, geometry):
    """Return each movable O, C, N and its measured N-carbon substituent(s)."""
    active = {
        frozenset(map(int, line.split()[1:]))
        for line in text.splitlines()
        if line.startswith("BRANCH ")
    }
    if geometry == "distinguished":
        return [s for s in amides(text) if frozenset(s[1:3]) in active]
    if geometry != "all-arms":
        raise ValueError(f"Unknown amide geometry: {geometry}")

    graph, serials = query_graph(text)
    sites = []
    for oxygen, carbon, nitrogen, groups in _carboxamides(graph):
        if frozenset((serials[carbon], serials[nitrogen])) not in active:
            continue
        # An imide nitrogen shares its conjugation between two acyl groups.
        if any(
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(arm).GetAtomicNum() == 8
            for arm in groups
            for bond in arm.GetBonds()
        ):
            continue
        arms = sorted(serials[a.GetIdx()] for a in groups)
        sites.append((serials[oxygen], serials[carbon], serials[nitrogen], *arms))
    return sorted(sites)


def check_barrier(barrier):
    if not math.isfinite(barrier) or barrier < 0:
        raise ValueError("Barrier must be a finite, nonnegative energy in kcal/mol")


def search_amide_records(text, barrier):
    """Build all-arms site rows for the search objective using prepared atom serials."""
    check_barrier(barrier)
    sites = active_amide_sites(text, "all-arms")
    if not sites:
        raise ValueError("No eligible movable carboxamide in the prepared ligand")
    return "".join(" ".join(map(str, site)) + f" {barrier!r}\n" for site in sites)


def dihedral_degrees(site, xyz):
    """Signed O--C--N--substituent angle in the prepared serial convention."""
    o, c, n, x = [np.asarray(xyz[i]) for i in site]
    axis = n - c
    length = np.linalg.norm(axis)
    if length < 1e-8:
        raise ValueError("Zero-length amide axis")
    axis = axis / length
    u, v = o - c, x - n
    u = u - np.dot(axis, u) * axis
    v = v - np.dot(axis, v) * axis
    if min(np.linalg.norm(u), np.linalg.norm(v)) < 1e-8:
        raise ValueError("Undefined amide torsion")
    return float(np.degrees(np.arctan2(np.dot(axis, np.cross(u, v)), np.dot(u, v))))


def coordinates(text):
    return {
        int(line[6:11]): np.asarray([float(line[i : i + 8]) for i in (30, 38, 46)])
        for line in text.splitlines()
        if line.startswith(("ATOM ", "HETATM "))
    }


def rotate(text, site, *, alternate=True, planar=False):
    """Rotate the nitrogen side, including its explicit hydrogens, about the C--N bond."""
    if site not in amides(text):
        raise ValueError(
            "The requested site is not a fixed, distinguishable acyclic carboxamide"
        )
    _, carbon, nitrogen, _ = site
    adjacency = defaultdict(set)
    for a, b in mapped_ligand_bonds(text):
        if {a, b} != {carbon, nitrogen}:
            adjacency[a].add(b)
            adjacency[b].add(a)
    side, todo = {nitrogen}, [nitrogen]
    while todo:
        atom = todo.pop()
        unseen = adjacency[atom] - side
        side.update(unseen)
        todo.extend(unseen)
    if carbon in side:
        raise ValueError("Cannot divide the ligand at a cyclic C--N bond")
    xyz = coordinates(text)
    center = xyz[nitrogen]
    axis = center - xyz[carbon]
    axis = axis / np.linalg.norm(axis)
    angle = np.pi if alternate else 0.0
    if planar:
        phi = np.deg2rad(dihedral_degrees(site, xyz))
        angle += np.rint(phi / np.pi) * np.pi - phi
    sn, cs = np.sin(angle), np.cos(angle)
    out = []
    for line in text.splitlines(keepends=True):
        if line.startswith(("ATOM ", "HETATM ")):
            idx = int(line[6:11])
            if idx in side and idx != nitrogen:
                offset = xyz[idx] - center
                moved = (
                    center
                    + cs * offset
                    + sn * np.cross(axis, offset)
                    + (1 - cs) * np.dot(offset, axis) * axis
                )
                field = "".join(f"{a:8.3f}" for a in moved)
                if len(field) != 24:
                    raise ValueError(
                        "Rotated coordinates exceed the fixed PDBQT field width"
                    )
                line = line[:30] + field + line[54:]
        out.append(line)
    return "".join(out)


def lock_amide(text, site):
    """Merge the C--N child into its parent rigid body; retain all other branches."""
    if site not in amides(text):
        raise ValueError("Selected C--N pair is not an eligible carboxamide")
    target = set(site[1:3])
    if not any(
        set(map(int, l.split()[1:])) == target
        for l in text.splitlines()
        if l.startswith("BRANCH ")
    ):
        return text
    header, footer = [], []
    root = {"atoms": [], "children": []}
    stack = [root]
    root_seen = False
    joined = 0
    for line in text.splitlines():
        tag = line.split(maxsplit=1)[0] if line.strip() else ""
        if tag == "ROOT":
            if root_seen:
                raise ValueError("Expected one ligand root")
            root_seen = True
        elif tag == "ENDROOT":
            continue
        elif tag == "BRANCH":
            pair = tuple(map(int, line.split()[1:]))
            node = {"atoms": [], "children": [], "pair": pair}
            stack.append(node)
        elif tag == "ENDBRANCH":
            if len(stack) < 2:
                raise ValueError("Unmatched ENDBRANCH")
            node = stack.pop()
            if tuple(map(int, line.split()[1:])) != node["pair"]:
                raise ValueError("Inconsistent branch delimiters")
            if set(node["pair"]) == target:
                joined += 1
                stack[-1]["atoms"].extend(node["atoms"])
                stack[-1]["children"].extend(node["children"])
            else:
                stack[-1]["children"].append(node)
        elif tag in ("ATOM", "HETATM"):
            stack[-1]["atoms"].append(line)
        elif tag == "TORSDOF":
            count = int(line.split()[1])
            if count < 1 or footer:
                raise ValueError("Bad TORSDOF")
            footer.append(f"TORSDOF {count - 1}")
        elif not root_seen:
            header.append(line)
        elif tag:
            raise ValueError(f"Unsupported record in ligand torsion tree: {tag}")
    if joined != 1 or not footer or len(stack) != 1:
        raise ValueError("Incomplete ligand torsion tree")

    def branch(node):
        i, j = node["pair"]
        yield f"BRANCH {i:3d} {j:3d}"
        yield from node["atoms"]
        for child in node["children"]:
            yield from branch(child)
        yield f"ENDBRANCH {i:3d} {j:3d}"

    lines = header + ["ROOT", *root["atoms"], "ENDROOT"]
    for child in root["children"]:
        lines.extend(branch(child))
    return "\n".join(lines + footer) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ligand", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        help="write one alternate ligand; omit to list eligible serials",
    )
    parser.add_argument(
        "--search-amide",
        type=Path,
        metavar="SITES",
        help="write movable O C N X [X2] B sites for all-arms amide search",
    )
    parser.add_argument(
        "--barrier",
        type=float,
        metavar="B",
        help="energy at a 90-degree twist in kcal/mol, required with --search-amide",
    )
    parser.add_argument(
        "--bond",
        type=int,
        nargs=2,
        metavar=("C", "N"),
        help="select carbonyl-C and amide-N by serial",
    )
    parser.add_argument(
        "--lock",
        action="store_true",
        help="freeze the selected C--N rotor in the output torsion tree",
    )
    parser.add_argument(
        "--keep-orientation",
        action="store_true",
        help="keep the input orientation when freezing the bond",
    )
    parser.add_argument(
        "--planar",
        action="store_true",
        help="project the O–C–N–substituent torsion to its chosen planar well",
    )
    args = parser.parse_args()
    if args.keep_orientation and not args.lock:
        parser.error("--keep-orientation requires --lock")
    text = args.ligand.read_text()
    if args.search_amide:
        if (
            args.output
            or args.bond
            or args.lock
            or args.keep_orientation
            or args.planar
        ):
            parser.error("--search-amide cannot be combined with rotamer generation")
        if args.barrier is None:
            parser.error("--search-amide requires --barrier")
        try:
            records = search_amide_records(text, args.barrier)
        except ValueError as exc:
            parser.error(str(exc))
        with args.search_amide.open("x") as out:
            out.write(records)
        return
    if args.barrier is not None:
        parser.error("--barrier requires --search-amide")
    sites = amides(text)
    if args.bond:
        sites = [s for s in sites if s[1:3] == tuple(args.bond)]
        if not sites:
            parser.error("Selected bond is not an eligible acyclic carboxamide")
    if args.output:
        if len(sites) != 1:
            parser.error(
                "Writing requires one eligible site; select --bond when there are several"
            )
        result = rotate(
            text, sites[0], alternate=not args.keep_orientation, planar=args.planar
        )
        if args.lock:
            result = lock_amide(result, sites[0])
        with args.output.open("x") as f:
            f.write(result)
    else:
        print("oxygen,carbon,nitrogen,substituent")
        for s in sites:
            print(",".join(map(str, s)))


if __name__ == "__main__":
    main()
