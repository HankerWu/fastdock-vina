"""Generate the alternate planar rotamer of a acyclic carboxamide in a prepared ligand."""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem

from make_scaffold_anchors import query_graph
from score_physical_pool import mapped_ligand_bonds


def amides(text):
    """Return (oxygen, carbonyl carbon, nitrogen, distinguished substituent) PDBQT serials."""
    graph, serials = query_graph(text)
    ranks = Chem.CanonicalRankAtoms(graph, breakTies=False)
    sites = []
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
        if len({ranks[a.GetIdx()] for a in groups}) != len(groups):
            continue
        tip = max(groups, key=lambda a: ranks[a.GetIdx()]).GetIdx()
        sites.append(tuple(serials[i] for i in (oxygen, carbon, nitrogen, tip)))
    return sites


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
