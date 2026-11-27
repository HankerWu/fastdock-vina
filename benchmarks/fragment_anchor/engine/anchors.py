"""Select ligand-side rigid anchors and reuse the original pocket triangle generator."""
from pathlib import Path
import sys

import numpy as np
from scipy.spatial.distance import pdist

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools"))
from pocket_triangle import _ligand_triangle


def rigid_components(ligand, bonds, serials):
    """Return heavy-atom connected components after cutting the native PDBQT torsion axes."""
    atoms, cuts = {}, set()
    for line in Path(ligand).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            serial = int(line[6:11])
            if serial in atoms:
                raise ValueError(f"duplicate PDBQT atom serial {serial}")
            atoms[serial] = line.split()[-1]
        elif line.startswith("BRANCH "):
            cuts.add(tuple(sorted(map(int, line.split()[1:]))))
    heavy = {s for s, t in atoms.items() if t not in ("H", "HD", "HS")}
    if set(serials) != heavy or len(serials) != len(heavy):
        raise ValueError("source heavy atom serials disagree with ligand")
    neighbors = {s: set() for s in heavy}
    edges = set()
    for line in Path(bonds).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        a, b = sorted(map(int, line.split()))
        if a not in atoms or b not in atoms or a == b:
            raise ValueError(f"invalid graph edge {a} {b}")
        if a in heavy and b in heavy:
            edges.add((a, b))
            if (a, b) not in cuts:
                neighbors[a].add(b)
                neighbors[b].add(a)
    if any(a in heavy and b in heavy and (a, b) not in edges for a, b in cuts):
        raise ValueError("missing native torsion axis in graph")
    remaining, components = set(heavy), []
    while remaining:
        stack = [min(remaining)]
        component = set()
        while stack:
            atom = stack.pop()
            if atom not in component:
                component.add(atom)
                stack.extend(neighbors[atom] - component)
        remaining -= component
        components.append(sorted(component))
    return sorted(components, key=lambda c: (-len(c), c))


def select_anchor(xyz, serials, ligand, bonds):
    """Return three source row indices, selected component serials, and fallback flag."""
    rows = {int(s): i for i, s in enumerate(serials)}
    for component in rigid_components(ligand, bonds, serials):
        if len(component) < 3:
            continue
        indices = np.array([rows[s] for s in component])
        try:
            anchor = _ligand_triangle(xyz[indices])
        except ValueError:  # collinear or coincident; try next rigid component
            continue
        return indices[anchor], np.asarray(component), False
    indices = np.argsort(serials)
    return indices[_ligand_triangle(xyz[indices])], np.asarray(serials)[indices], True


def all_transforms(pocket, xyz, ligand_triangle):
    """Evaluation only, all site triples in original enumeration order before any edge gate."""
    source = xyz[ligand_triangle]
    centers, rotations = pocket.fit_triangles(xyz, source, pocket.triangles)
    return dict(center=centers, matched_rotation=rotations, site_triangles=pocket.triangles,
                edge_residual_A2=((pocket.edges - pdist(source)) ** 2).sum(axis=1))


def place(pocket, xyz, serials, ligand, bonds, seed, arm):
    """Generate a fixed-size bank for a single arm; matched_* is the actual chosen geometry."""
    x, serials = np.asarray(xyz, float), np.asarray(serials, int)
    if arm == "whole":
        order = np.argsort(serials)
        anchor = order[_ligand_triangle(x[order])]
        component, fallback = serials[order], False
    elif arm in ("fragment", "sham"):
        anchor, component, fallback = select_anchor(x, serials, ligand, bonds)
    else:
        raise ValueError(f"unknown geometry {arm}")
    kind = "sham" if arm == "sham" else "matched"
    result = pocket.place(x, seed, arms=(kind,), ligand_triangle=anchor)
    if kind == "sham":
        for key in ("rotation", "score", "best"):
            result["matched_" + key] = result.pop("sham_" + key)
    if len(result["center"]) != 1024:
        raise ValueError("this experiment requires exactly 1024 cavity correspondences")
    # Identify each prefilter triangle without retaining a second, redundant transform bank.
    n = len(pocket.sites)
    weights = np.array([n*n, n, 1])
    result.update(anchor_serials=serials[anchor], component_serials=component, anchor_fallback=fallback,
                  triangle_id=np.searchsorted(pocket.triangles @ weights, result["site_triangles"] @ weights))
    return result
