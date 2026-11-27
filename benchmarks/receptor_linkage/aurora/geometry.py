"""Cross-docked AURKA ligands in independent crystal frames, without ligand fitting."""

import argparse
import csv
import json
import tarfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation

REFERENCES = {
    "danusertib": (("2J50", "627", ("A", "B")),),
    "at9283": (("2W1G", "L0G", ("A",)), ("6CPG", "35R", ("A", "D"))),
}
RECEPTORS = {"2J4Z": "626", "2X6E": "YM4"}
METHODS = ("official", "pro", "double")
SEEDS = (60913, 60943, 61211, 61231)


def xyz(line):
    return [float(line[k : k + 8]) for k in (30, 38, 46)]


def write_csv(out, name, rows):
    with (out / name).open("w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)


def ca_atoms(lines, chain):
    # Mirror the source structure preparation: blank/A alternate site.
    atoms = {}
    for line in lines:
        if (
            line.startswith("ATOM  ")
            and line[21] == chain
            and line[16] in " A"
            and line[12:16].strip() == "CA"
        ):
            key = line[22:27].strip(), line[17:20].strip()
            if key in atoms:
                raise ValueError(("repeated backbone position", key))
            atoms[key] = np.array(xyz(line))
    return atoms


def ligand_atoms(lines, component, chain):
    atoms = {}
    for line in lines:
        if (
            line.startswith("HETATM")
            and line[21] == chain
            and line[17:20] == component
            and line[76:78].strip() not in {"H", "D"}
            and line[16] in " A"
        ):
            name = line[12:16].strip()
            if name in atoms:
                raise ValueError(
                    ("multiple copies/sites with same atom", component, chain, name)
                )
            atoms[name] = (line[76:78].strip(), int(line[6:11]), xyz(line))
    return atoms


def receiver_frames(archive):
    frames = {}
    with tarfile.open(archive) as f:
        for pdb, component in RECEPTORS.items():
            raw = f.extractfile("raw/" + pdb + ".pdb").read().decode().splitlines()
            prepared = (
                f.extractfile("prepared/" + pdb + "-protein.pdb")
                .read()
                .decode()
                .splitlines()
            )
            lig = np.array([a[2] for a in ligand_atoms(raw, component, "A").values()])
            raw_ca = ca_atoms(raw, "A")
            processed_ca = ca_atoms(prepared, "A")
            selected = {
                key: pos
                for key, pos in raw_ca.items()
                if np.min(np.linalg.norm(lig - pos, axis=1)) <= 12
            }
            for key, pos in selected.items():
                if (
                    key not in processed_ca
                    or np.linalg.norm(pos - processed_ca[key]) > 0.002
                ):
                    raise ValueError(
                        ("prepared pocket disagrees with raw PDB", pdb, key)
                    )
            frames[pdb] = selected
    return frames


def transform(target_ca, ref_ca):
    keys = sorted(target_ca.keys() & ref_ca.keys())
    if len(keys) < 8:
        raise ValueError("Insufficient CA matches; sequence alignment is required")
    dest = np.array([target_ca[k] for k in keys])
    source = np.array([ref_ca[k] for k in keys])
    sc, dc = source.mean(axis=0), dest.mean(axis=0)
    x, y = source - sc, dest - dc
    u, _, vh = np.linalg.svd(x.T @ y)
    sign = np.linalg.det(u @ vh)
    matrix = u @ np.diag([1, 1, sign]) @ vh
    rotational, _ = Rotation.align_vectors(y, x)
    rotation_gap = np.max(np.abs(rotational.as_matrix().T - matrix))
    if rotation_gap > 1e-10 or np.linalg.det(matrix) <= 0:
        raise ArithmeticError("Independent protein rotations disagree")
    shift = dc - sc @ matrix
    delta = source @ matrix + shift - dest
    return (
        matrix,
        shift,
        {
            "ca_residues": keys,
            "aligned_ca": len(keys),
            "rmsd_A": float(np.sqrt(np.mean(np.sum(delta**2, axis=1)))),
            "max_rotation_matrix_error": float(rotation_gap),
        },
    )


def single_bond_graph(elements, pairs):
    graph = Chem.RWMol()
    for elem in elements:
        a = Chem.Atom(elem)
        a.SetNoImplicit(True)
        graph.AddAtom(a)
    for i, j in pairs:
        graph.AddBond(i, j, Chem.BondType.SINGLE)
    return graph.GetMol()


def ccd_data(directory, component):
    block = gemmi.cif.read_file(str(directory / (component + ".cif"))).sole_block()
    atoms = [
        (a, e)
        for a, e in block.find("_chem_comp_atom.", ["atom_id", "type_symbol"])
        if e not in {"H", "D"}
    ]
    names = [a for a, _ in atoms]
    bonds = [
        (a, b)
        for a, b in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2"])
        if a in names and b in names
    ]
    inchikey = next(
        d
        for t, d in block.find("_pdbx_chem_comp_descriptor.", ["type", "descriptor"])
        if t == "InChIKey"
    )
    return (names, [e for _, e in atoms], bonds, inchikey)


def reference_mappings(molecule, atoms, ccd, lines):
    names, elements, bonds, crystal_key = ccd
    pred_key = Chem.MolToInchiKey(molecule)
    # The last character distinguishes neutral and protonated parents.
    if pred_key[:-1] != crystal_key[:-1]:
        raise ValueError(
            ("chemical identity or stereo mismatch", pred_key, crystal_key)
        )
    if set(names) != set(atoms) or any(
        atoms[n][0] != e for n, e in zip(names, elements)
    ):
        raise ValueError("Incomplete crystal heavy-atom identity")
    indices = {atoms[n][1]: n for n in names}
    pdb_bonds = set()
    for line in lines:
        if line.startswith("CONECT"):
            serials = [
                int(line[k : k + 5])
                for k in range(6, len(line.rstrip()), 5)
                if line[k : k + 5].strip()
            ]
            for other in serials[1:]:
                if serials[0] in indices and other in indices:
                    pdb_bonds.add(frozenset((indices[serials[0]], indices[other])))
    if pdb_bonds != {frozenset(edge) for edge in bonds}:
        raise ValueError("PDB CONECT and CCD heavy graphs differ")
    reference = single_bond_graph(
        elements, ((names.index(a), names.index(b)) for a, b in bonds)
    )
    graph = single_bond_graph(
        [a.GetSymbol() for a in molecule.GetAtoms()],
        [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in molecule.GetBonds()],
    )
    matches = reference.GetSubstructMatches(graph, uniquify=False, maxMatches=10000)
    if not matches or graph.GetNumAtoms() != reference.GetNumAtoms():
        raise ValueError("No full graph isomorphism")
    chemical_autos = molecule.GetSubstructMatches(
        molecule, uniquify=False, useChirality=True, maxMatches=10000
    )
    accepted = {tuple(matches[0][j] for j in auto) for auto in chemical_autos}
    if accepted != set(matches) or len(matches) == 10000:
        raise ValueError(
            "Element-only graph includes non-equivalent maps; require typed adjudication"
        )
    crystal_coords = np.array([atoms[n][2] for n in names])
    # Signed volumes compare the same three ordered neighbours at each specified
    # tetrahedral centre, independent of protonation, reference rotation and origin.
    ref_volumes = []
    input_coords = molecule.GetConformer().GetPositions()
    for k, tag in Chem.FindMolChiralCenters(molecule, includeUnassigned=True):
        if tag == "?":
            raise ValueError("Unspecified prepared ligand chirality")
        neighbors = [n.GetIdx() for n in molecule.GetAtomWithIdx(k).GetNeighbors()][:3]
        v_in = np.linalg.det(input_coords[neighbors] - input_coords[k])
        for match in matches:
            c = crystal_coords[list(match)]
            v_ref = np.linalg.det(c[neighbors] - c[k])
            if v_in * v_ref <= 0:
                raise ValueError(("Inverted tetrahedral geometry", k))
        ref_volumes.append((k, tag, float(v_in), float(v_ref)))
    return crystal_coords[np.array(matches)], {
        "query_inchikey": pred_key,
        "ccd_inchikey": crystal_key,
        "heavy_atoms": len(names),
        "graph_mappings": len(matches),
        "checked_stereo": ref_volumes,
        "ccd_names_in_query_order": [[names[j] for j in match] for match in matches],
    }


def serial_map(ligand_path, mol):
    lines = ligand_path.read_text().splitlines()
    mapping = {}
    for l in lines:
        if l.startswith("REMARK INDEX MAP "):
            vals = list(map(int, l.split()[3:]))
            mapping.update(zip(vals[::2], vals[1::2]))
    atoms = {int(l[6:11]): xyz(l) for l in lines if l.startswith(("ATOM  ", "HETATM"))}
    # SDF heavy atoms precede hydrogens in the frozen prepared inputs. Verify each
    # mapped heavy coordinate against the actual SDF before applying to output modes.
    mapping = [mapping[i] for i in range(1, mol.GetNumAtoms() + 1)]
    deviation = np.abs(
        np.array([atoms[i] for i in mapping]) - mol.GetConformer().GetPositions()
    ).max()
    if deviation > 0.00061 or len(set(mapping)) != len(mapping):
        raise ValueError(("PDBQT/SDF atom correspondence", ligand_path.name, deviation))
    return mapping, float(deviation)


def mode_coords(path, serials):
    modes = []
    for section in path.read_text().split("ENDMDL"):
        lines = section.splitlines()
        coords = {
            int(l[6:11]): xyz(l) for l in lines if l.startswith(("ATOM  ", "HETATM"))
        }
        if coords:
            modes.append(np.array([coords[i] for i in serials]))
    if not modes:
        raise ValueError(("empty pool", path))
    return np.array(modes)


def summary_rows(cells):
    strata = defaultdict(list)
    for r in cells:
        keys = [
            (r["reference"], r["chain"], r["receptor"], r["method"]),
            (r["reference"], r["chain"], "both", r["method"]),
        ]
        if r["primary"]:
            keys.extend(
                [
                    ("primary", "A", r["receptor"], r["method"]),
                    ("primary", "A", "both", r["method"]),
                ]
            )
        for key in keys:
            strata[key].append(r)
    for (ref, chain, rec, method), rows in sorted(strata.items()):
        for col in ("top1", "best5", "oracle", "regret"):
            x = [r[col] for r in rows]
            q25, med, q75, p90, p95 = np.percentile(x, [25, 50, 75, 90, 95])
            yield {
                "reference": ref,
                "chain": chain,
                "receptor": rec,
                "method": method,
                "metric": col,
                "n": len(x),
                "mean": np.mean(x),
                "q25": q25,
                "median": med,
                "q75": q75,
                "p90": p90,
                "p95": p95,
                "maximum": max(x),
                "above10": sum(v > 10 for v in x),
            }


def evaluate(run, ref_dir, output, kinase_archive):
    frames = receiver_frames(kinase_archive)
    modes, cells, paired, checks = [], [], [], {}
    for drug, entries in REFERENCES.items():
        ligand = run / "prepared" / drug
        molecule = Chem.MolFromMolFile(str(ligand.with_suffix(".sdf")), removeHs=True)
        serials, err = serial_map(ligand.with_suffix(".pdbqt"), molecule)
        for ref, component, chains in entries:
            lines = (ref_dir / (ref + ".pdb")).read_text().splitlines()
            for chain in chains:
                ccd = ccd_data(ref_dir, component)
                mapped, check = reference_mappings(
                    molecule, ligand_atoms(lines, component, chain), ccd, lines
                )
                check.update(receptors={}, sdf_pdbqt_max_abs_A=err)
                checks[f"{ref}_{chain}"] = check
                for rec, pocket in frames.items():
                    matrix, shift, note = transform(pocket, ca_atoms(lines, chain))
                    check["receptors"][rec] = note
                    transformed = mapped @ matrix + shift
                    for seed in SEEDS:
                        case = {
                            "compound": drug,
                            "reference": ref,
                            "chain": chain,
                            "primary": ref != "6CPG" and chain == "A",
                            "receptor": rec,
                            "seed": seed,
                        }
                        arms = {}
                        for method in METHODS:
                            path = run / "pools" / f"{method}_{rec}_{drug}_{seed}.pdbqt"
                            poses = mode_coords(path, serials)
                            squared = (
                                poses[:, None, :, :] - transformed[None, :, :, :]
                            ) ** 2
                            rmsds = np.sqrt(
                                squared.sum(axis=-1).mean(axis=-1).min(axis=-1)
                            )
                            for i, val in enumerate(rmsds, 1):
                                modes.append(
                                    dict(**case, method=method, rank=i, rmsd=val)
                                )
                            scores = {
                                "top1": rmsds[0],
                                "best5": min(rmsds[:5]),
                                "oracle": min(rmsds),
                                "regret": rmsds[0] - min(rmsds),
                            }
                            arms[method] = scores
                            cells.append(
                                dict(**case, method=method, count=len(rmsds), **scores)
                            )
                        for before, after in combinations(METHODS, 2):
                            paired.append(
                                dict(
                                    **case,
                                    before=before,
                                    after=after,
                                    **{
                                        k: arms[after][k] - v
                                        for k, v in arms[before].items()
                                    },
                                )
                            )
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output, "geometry_modes.csv", modes)
    write_csv(output, "geometry_cells.csv", cells)
    write_csv(output, "geometry_pairs.csv", paired)
    write_csv(output, "geometry_distributions.csv", list(summary_rows(cells)))
    (output / "geometry_checks.json").write_text(json.dumps(checks, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    p.add_argument("references", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--kinase-archive", type=Path, required=True)
    a = p.parse_args()
    evaluate(a.run, a.references, a.output, a.kinase_archive)
