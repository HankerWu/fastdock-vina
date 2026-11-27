"""Reference-informed placement and relaxation on the fixed Aurora A panel."""

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from geometry import (
    REFERENCES,
    SEEDS,
    ca_atoms,
    ccd_data,
    ligand_atoms,
    mode_coords,
    receiver_frames,
    reference_mappings,
    serial_map,
    transform,
    write_csv,
    xyz,
)
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


def get_tree(path):
    lines = path.read_text().splitlines()
    ids, sdf_map, branches, stack = [], {}, [], []
    for l in lines:
        if l.startswith("REMARK INDEX MAP "):
            v = list(map(int, l.split()[3:]))
            sdf_map.update(zip(v[1::2], [i - 1 for i in v[::2]]))
        if l.startswith("BRANCH "):
            a, b = map(int, l.split()[1:])
            node = (a, b, [])
            branches.append(node)
            stack.append(node)
        if l.startswith("ENDBRANCH "):
            node = stack.pop()
            assert tuple(map(int, l.split()[1:])) == node[:2]
        if l.startswith(("ATOM  ", "HETATM")):
            ids.append(int(l[6:11]))
            for _, _, nodes in stack:
                nodes.append(len(ids) - 1)
    assert not stack and len(ids) == len(sdf_map)
    pts = np.array([xyz(l) for l in lines if l.startswith(("ATOM  ", "HETATM"))])
    arms = [(ids.index(a), ids.index(b), np.array(m)) for a, b, m in branches]
    return lines, ids, np.array([sdf_map[i] for i in ids]), pts, arms


def register(points, heavy, target):
    x, y = points[heavy], target
    u, _, vt = np.linalg.svd((x - x.mean(axis=0)).T @ (y - y.mean(axis=0)))
    R = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
    return (points - x.mean(axis=0)) @ R + y.mean(axis=0)


def fit_torsions(points, heavy, arms, templates, *, accept_finite_limit=False):
    def oriented(angles, target):
        x = points.copy()
        for theta, (a, b, nodes) in zip(angles, arms):
            axis = x[b] - x[a]
            x[nodes] = (
                Rotation.from_rotvec(axis * (theta / np.linalg.norm(axis))).apply(
                    x[nodes] - x[a]
                )
                + x[a]
            )
        return register(x, heavy, target)

    # Construct an informed start from each observable branch dihedral.
    # A rigid rotation about axis a->b aligns the projected outgoing and
    # incoming neighbour directions, preserving the prepared ring geometry.
    inv = {j: i for i, j in enumerate(heavy)}

    def twist(x, a, b, c, d):
        axis = x[b] - x[a]
        axis /= np.linalg.norm(axis)
        u, v = x[c] - x[a], x[d] - x[b]
        u -= axis * (u @ axis)
        v -= axis * (v @ axis)
        return np.arctan2(np.cross(u, v) @ axis, u @ v)

    def residual(angles, target):
        return (oriented(angles, target)[heavy] - target).ravel()

    best = None
    for t in templates:
        init = []
        for a, b, nodes in arms:
            inner = set(nodes)
            outer_heavy = [i for i in heavy if i != a and i not in inner]
            inner_heavy = [i for i in heavy if i != b and i in inner]
            c = min(outer_heavy, key=lambda i: np.linalg.norm(points[i] - points[a]))
            d = min(inner_heavy, key=lambda i: np.linalg.norm(points[i] - points[b]))
            init.append(
                twist(t, inv[a], inv[b], inv[c], inv[d]) - twist(points, a, b, c, d)
            )
        for first in (init, np.zeros(len(arms))):
            res = least_squares(
                residual,
                first,
                args=(t,),
                max_nfev=200,
                ftol=1e-10,
                xtol=1e-10,
                gtol=1e-10,
            )
            finite_limit = accept_finite_limit and res.status == 0 and np.all(np.isfinite(res.fun))
            if not res.success and not finite_limit:
                raise ArithmeticError(res.message)
            rmsd = np.linalg.norm(res.fun) / np.sqrt(len(heavy))
            if best is None or rmsd < best[0]:
                best = rmsd, oriented(res.x, t), res.x
    return best


def crystal_geometry(sdf, ref, sdf_order):
    # Replace the prepared heavy geometry, then place hydrogens using this
    # SAME neutral graph. Heavy atoms remain fixed throughout MMFF minimization.
    heavy = Chem.RemoveHs(sdf)
    for i, pos in enumerate(ref):
        heavy.GetConformer().SetAtomPosition(i, pos)
    mol = Chem.AddHs(heavy, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    for i in range(len(ref)):
        ff.AddFixedPoint(i)
    ff.Initialize()
    status = ff.Minimize(maxIts=1500)
    assert status == 0
    xyz_new = mol.GetConformer().GetPositions()
    assert np.abs(xyz_new[: len(ref)] - ref).max() < 1e-8
    coords = []
    for i in sdf_order:
        atom = sdf.GetAtomWithIdx(int(i))
        if atom.GetAtomicNum() == 1:
            a = atom.GetNeighbors()[0].GetIdx()
            old = [
                b.GetIdx()
                for b in sdf.GetAtomWithIdx(a).GetNeighbors()
                if b.GetAtomicNum() == 1
            ]
            new = [
                b.GetIdx()
                for b in mol.GetAtomWithIdx(a).GetNeighbors()
                if b.GetAtomicNum() == 1
            ]
            i = new[old.index(i)]
        coords.append(xyz_new[i])
    return np.array(coords)


def store_pose(path, lines, coords):
    output = lines.copy()
    atoms = [i for i, l in enumerate(lines) if l.startswith(("ATOM  ", "HETATM"))]
    assert len(atoms) == len(coords)
    for i, p in zip(atoms, coords):
        output[i] = lines[i][:30] + "".join(f"{a:8.3f}" for a in p) + lines[i][54:]
    path.write_text("\n".join(output) + "\n")


def run_logged(command, log):
    result = subprocess.run(
        list(map(str, command)),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{log}: {result.stdout[-1900:]}")
    return result.stdout


def score_pose(score_command, path):
    output = run_logged(
        [
            *score_command,
            "--pose_coordinates",
            path,
            "--score_only",
            "--unbound_energy",
            "0",
        ],
        path.with_name(path.stem + "_score.log"),
    )
    match = re.search(r"Estimated Free Energy of Binding\s*:\s*([-+.\d]+)", output)
    if not match:
        raise ValueError(output)
    return float(match[1])


def local_refinement(command, path, bonds, output, log, options=()):
    return run_logged(
        [
            *command,
            "--ligand",
            path,
            "--ligand_bonds",
            bonds,
            "--local_only",
            *options,
            "--out",
            output,
        ],
        log,
    )


def initial_objective(command, path, bonds, serials):
    """Read the optimizer objective before any coordinate update."""
    static = path.with_name(path.stem + "_one_eval.pdbqt")
    message = local_refinement(
        command,
        path,
        bonds,
        static,
        static.with_suffix(".log"),
        ["--refinement_optimizer", "armijo", "--max_evals", "1"],
    )
    report = re.search(r"gradient_evals=1 objective_kcal_mol=([^\s]+)", message)
    if not report or not np.array_equal(
        mode_coords(path, serials), mode_coords(static, serials)
    ):
        raise ValueError(
            f"Expected one gradient evaluation and unchanged coordinates: {path}"
        )
    return float(report[1])


def mapped_rmsd(pose, references):
    return np.sqrt(
        np.min(np.mean(np.sum((pose[None] - references) ** 2, axis=-1), axis=-1))
    )


def max_distance_change(reference, coords, pairs):
    return max(
        abs(
            np.linalg.norm(coords[i] - coords[j])
            - np.linalg.norm(reference[i] - reference[j])
        )
        for i, j in pairs
    )


def run(input_dir, refs, out, archive, exe, parser):

    out.mkdir(exist_ok=True, parents=True)
    raw = out / "raw"
    raw.mkdir(exist_ok=True)
    frames = receiver_frames(archive)
    stage_rows, pair_rows, audits = [], [], {}

    for drug, entries in REFERENCES.items():
        reference, component, _ = entries[0]
        lines = (refs / (reference + ".pdb")).read_text().splitlines()
        source = input_dir / "prepared" / (drug + ".pdbqt")
        sdf = Chem.MolFromMolFile(str(source.with_suffix(".sdf")), removeHs=False)
        mol = Chem.RemoveHs(sdf)
        templates, _ = reference_mappings(
            mol, ligand_atoms(lines, component, "A"), ccd_data(refs, component), lines
        )
        template, ids, sdf_order, points, arms = get_tree(source)
        by_sdf = {i: j for j, i in enumerate(sdf_order)}
        heavy = np.array([by_sdf[i] for i in range(mol.GetNumAtoms())])
        parsed = subprocess.check_output([parser, source])
        (raw / (drug + "_prepared_topology.txt")).write_bytes(parsed)
        native_bonds = raw / (drug + "_native.bonds")
        native_bonds.write_text(
            "".join(
                l[5:] + "\n"
                for l in parsed.decode().splitlines()
                if l.startswith("bond ")
            )
        )
        if subprocess.check_output([parser, source, native_bonds]) != parsed:
            raise ValueError("Explicit native bonds changed preparing graph")
        best, fitted, angles = fit_torsions(points, heavy, arms, templates)
        rigid = min(
            np.linalg.norm(register(points, heavy, t)[heavy] - t) / np.sqrt(len(heavy))
            for t in templates
        )
        exact = crystal_geometry(sdf, templates[0], sdf_order)
        heavy_serials, _ = serial_map(source, mol)
        # Independent Cartesian distance checks on all covalent edges and all
        # angles in the input graph, including retained polar hydrogens.
        pairs = set()
        for b in sdf.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i in by_sdf and j in by_sdf:
                pairs.add(tuple(sorted([by_sdf[i], by_sdf[j]])))
        angular = set()
        for a in sdf.GetAtoms():
            if a.GetIdx() not in by_sdf:
                continue
            ns = [by_sdf[x.GetIdx()] for x in a.GetNeighbors() if x.GetIdx() in by_sdf]
            for i in ns:
                for j in ns:
                    if i < j:
                        angular.add((i, j))

        audits[drug] = {
            "rigid_only_RMSD": rigid,
            "legal_torsion_fit_RMSD": best,
            "torsion_radians": angles.tolist(),
            "branch_count": len(arms),
            "prepared_atom_serials": ids,
            "bonds_checked": len(pairs),
            "one_three_checked": len(angular),
            "legal_max_bond_change_A": max_distance_change(points, fitted, pairs),
            "legal_max_one_three_change_A": max_distance_change(
                points, fitted, angular
            ),
            "crystal_max_bond_change_A": max_distance_change(points, exact, pairs),
            "crystal_max_one_three_change_A": max_distance_change(
                points, exact, angular
            ),
        }
        assert (
            max_distance_change(points, fitted, pairs) < 1e-10
            and max_distance_change(points, fitted, angular) < 1e-10
        )
        # Same prepared native graph for all arms. If its bond inference
        # differs from the SDF, repeat the same stored inputs on the full graph.
        native_edges = {
            tuple(map(int, line.split()))
            for line in native_bonds.read_text().splitlines()
        }
        chemical_edges = {tuple(sorted((ids[i], ids[j]))) for i, j in pairs}
        audits[drug].update(
            missing_edges=sorted(chemical_edges - native_edges),
            extra_edges=sorted(native_edges - chemical_edges),
        )
        graphs = {"native": native_bonds}
        if chemical_edges != native_edges:
            chem = raw / (drug + "_chemical.bonds")
            chem.write_text("".join(f"{a} {b}\n" for a, b in sorted(chemical_edges)))
            graphs["sdf"] = chem
        for rec, ca in frames.items():
            R, shift, _ = transform(ca, ca_atoms(lines, "A"))
            refs_in_frame = templates @ R + shift
            center = json.loads((input_dir / "prepared" / (rec + ".json")).read_text())[
                "center"
            ]
            common = [
                exe,
                "--receptor",
                input_dir / "prepared" / (rec + ".pdbqt"),
                "--cpu",
                "1",
                "--seed",
                "60913",
            ]
            for axis, v in zip("xyz", center):
                common += ["--center_" + axis, str(v), "--size_" + axis, "26"]
            rec_lines = [
                l
                for l in (input_dir / "prepared" / (rec + ".pdbqt"))
                .read_text()
                .splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
                and l.split()[-1] not in {"H", "HD", "HS"}
            ]
            rec_coords = np.array([xyz(l) for l in rec_lines])
            for graph, bond_file in graphs.items():
                signature_ref = subprocess.check_output([parser, source, bond_file])

                score_command = [*common, "--ligand", source]
                if graph != "native":
                    score_command += ["--ligand_bonds", bond_file]
                winners = {}
                for seed in SEEDS:
                    blind = (
                        input_dir / "pools" / f"pro_{rec}_{drug}_{seed}.pdbqt"
                    ).read_text()
                    path = raw / f"{rec}_{drug}_{graph}_{seed}_winner.pdbqt"
                    path.write_text(blind.split("ENDMDL", 1)[0] + "ENDMDL\n")
                    winners[seed] = score_pose(score_command, path)
                for label, arr in [("crystal", exact), ("legal", fitted)]:
                    stem = f"{rec}_{drug}_{graph}_{label}"
                    original = raw / (stem + "_placed.pdbqt")
                    store_pose(original, template, arr @ R + shift)
                    signature = subprocess.check_output([parser, original])
                    (raw / (stem + "_topology.txt")).write_bytes(signature)
                    locked = subprocess.check_output([parser, original, bond_file])
                    if locked != signature_ref:
                        raise ValueError(
                            f"Locked graph differs from prepared input for {original}"
                        )
                    if graph == "native":
                        audits[drug].setdefault("placed_native_inference_equal", {})[
                            stem
                        ] = signature == parsed
                    relaxed = raw / (stem + "_relaxed.pdbqt")
                    local_refinement(
                        common,
                        original,
                        bond_file,
                        relaxed,
                        raw / (stem + "_local.log"),
                    )
                    for stage, path in (("placed", original), ("relaxed", relaxed)):
                        pose = mode_coords(path, heavy_serials)[0]
                        rmsd = mapped_rmsd(pose, refs_in_frame)
                        energy = score_pose(score_command, path)
                        distances = np.linalg.norm(
                            pose[:, None] - rec_coords[None, :], axis=-1
                        )
                        ai, bi = np.unravel_index(np.argmin(distances), distances.shape)
                        edge = float(np.min(13 - np.abs(pose - np.array(center))))
                        stage_rows.append(
                            {
                                "receptor": rec,
                                "compound": drug,
                                "graph": graph,
                                "witness": label,
                                "stage": stage,
                                "rmsd_A": rmsd,
                                "zero_reference_score": energy,
                                "grid_edge_margin_A": edge,
                                "ligand_closest_serial": heavy_serials[ai],
                                "receiver_closest_atom": rec_lines[bi][12:26].strip(),
                                "closest_heavy_A": distances[ai, bi],
                                "inter_heavy_pairs_below_2A": int(
                                    np.sum(distances < 2)
                                ),
                            }
                        )
                        for seed, s in winners.items():
                            pair_rows.append(
                                {
                                    "receptor": rec,
                                    "compound": drug,
                                    "graph": graph,
                                    "witness": label,
                                    "stage": stage,
                                    "seed": seed,
                                    "winner_score": s,
                                    "witness_score": energy,
                                    "winner_minus_witness": s - energy,
                                }
                            )
                    # Test whether a completely unmoved native local result is
                    # stationary using rigid translations and the existing Armijo search.
                    start_xyz = mode_coords(original, ids)[0]
                    final_xyz = mode_coords(relaxed, ids)[0]
                    if np.array_equal(start_xyz, final_xyz):
                        base_obj = initial_objective(common, original, bond_file, ids)
                        derivatives = []
                        for axis in range(3):
                            scores = []
                            objectives = []
                            for delta in (-0.02, 0.02):
                                xyz_shifted = start_xyz.copy()
                                xyz_shifted[:, axis] += delta
                                probe = raw / f"{stem}_shift{axis}_{delta}.pdbqt"
                                store_pose(probe, template, xyz_shifted)
                                scores.append(score_pose(score_command, probe))
                                objectives.append(
                                    initial_objective(common, probe, bond_file, ids)
                                )
                            derivatives.append(
                                {
                                    "axis": "xyz"[axis],
                                    "minus_score": scores[0],
                                    "plus_score": scores[1],
                                    "central_slope_per_A": (scores[1] - scores[0])
                                    / 0.04,
                                    "minus_objective": objectives[0],
                                    "plus_objective": objectives[1],
                                    "objective_secant_per_A": (
                                        objectives[1] - objectives[0]
                                    )
                                    / 0.04,
                                }
                            )
                        alternative = raw / (stem + "_armijo.pdbqt")
                        local_refinement(
                            common,
                            original,
                            bond_file,
                            alternative,
                            raw / (stem + "_armijo.log"),
                            ["--refinement_optimizer", "armijo"],
                        )
                        alt_pose = mode_coords(alternative, heavy_serials)[0]
                        alt_rmsd = mapped_rmsd(alt_pose, refs_in_frame)
                        alt_score = score_pose(score_command, alternative)
                        audits[drug].setdefault("unchanged_native_local", {})[stem] = {
                            "initial_objective": base_obj,
                            "central_translation": derivatives,
                            "armijo_rmsd_A": alt_rmsd,
                            "armijo_score": alt_score,
                            "winner_minus_armijo": {
                                str(seed): s - alt_score for seed, s in winners.items()
                            },
                        }
    write_csv(out, "probe_stages.csv", stage_rows)
    write_csv(out, "probe_pairs.csv", pair_rows)
    audits["binary_sha256"] = hashlib.sha256(exe.read_bytes()).hexdigest()
    (out / "probe_checks.json").write_text(json.dumps(audits, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "references", "output", "archive", "binary", "parser"):
        p.add_argument(name, type=lambda x: Path(x).resolve())
    a = p.parse_args()
    run(a.run, a.references, a.output, a.archive, a.binary, a.parser)
