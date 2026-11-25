"""Freeze shape-qualified edges and recompute RBFE from OpenFE per-replica leg records."""
import argparse
import bz2
import csv
from datetime import datetime, timezone
import hashlib
import json
from math import hypot, isfinite, log, sqrt
from pathlib import Path
import os
import re
import shutil
import subprocess
from statistics import mean, stdev
from time import perf_counter

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign
import yaml


SCORES = ("native_score", "burial", "geometry", "compensated")
R_KCAL = 0.00198720425864083


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def put_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def put_csv(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)


def alpha_carbons(path, chain):
    return {
        (line[17:20], line[22:27]): [float(line[i:i+8]) for i in (30, 38, 46)]
        for line in path.read_text().splitlines()
        if line.startswith("ATOM ") and line[12:16] == " CA " and line[21] == chain and line[16] in " A"
    }


def receptor_alignment(sources, config):
    old = alpha_carbons(sources["historical_protein"], config["historical_chain"])
    new = alpha_carbons(sources["openfe_protein"], config["openfe_chain"])
    keys = sorted(old.keys() & new.keys())
    if len(keys) < config["minimum_protein_ca_matches"]:
        raise ValueError("Insufficient same-residue C-alpha matches")
    left, right = (np.asarray([atoms[key] for key in keys]) for atoms in (old, new))
    old_center, new_center = left.mean(axis=0), right.mean(axis=0)
    cov = (right-new_center).T @ (left-old_center)
    u, _, vt = np.linalg.svd(cov)
    handedness = np.eye(3)
    handedness[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ handedness @ vt
    translate = old_center - new_center @ rotation
    rmsd = float(np.sqrt(np.mean(np.sum((right @ rotation + translate - left) ** 2, axis=1))))
    if rmsd > config["maximum_protein_alignment_rmsd_angstrom"]:
        raise ValueError(f"Receptor alignment RMSD {rmsd:.3f} exceeds threshold")
    return rotation, translate, dict(ca_matches=len(keys), receptor_ca_rmsd=rmsd,
                                     rotation=rotation.tolist(), translation=translate.tolist())


def ligand_molecules(path):
    items = list(Chem.SDMolSupplier(str(path), removeHs=False))
    if any(item is None for item in items):
        raise ValueError(f"Unparseable 3D ligand in {path}")
    result = {item.GetProp("_Name"): Chem.RemoveHs(item) for item in items}
    if len(result) != len(items):
        raise ValueError(f"Repeated ligand name in {path}")
    return result


def native_pose(pdbqt):
    import meeko

    model_pattern = r"^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$"
    blocks = re.findall(model_pattern, pdbqt.read_text(), re.M | re.S)
    rank_1 = [(model, block) for model, block in blocks
              if re.search(r"^REMARK VINA NATIVE: \S+ 1\s*$", block, re.M)]
    (model, block), = rank_1
    recovered, = meeko.RDKitMolCreate.from_pdbqt_mol(meeko.PDBQTMolecule(block))
    return int(model), Chem.RemoveHs(recovered)


def active_atoms(molecule, query):
    return molecule.GetSubstructMatches(query, useChirality=True, uniquify=False, maxMatches=2000)


def fixed_pose_rmsd(pose, target, rotation, translation):
    atoms = active_atoms(pose, target)
    if not atoms:
        raise ValueError("Prepared pose and analog identity differ")
    mapping = [[(idx, i) for i, idx in enumerate(atom_map)] for atom_map in atoms]
    transformed = Chem.Mol(target)
    for i, coord in enumerate(target.GetConformer().GetPositions() @ rotation + translation):
        transformed.GetConformer().SetAtomPosition(i, coord)
    return float(rdMolAlign.CalcRMS(pose, transformed, map=mapping, symmetrizeConjugatedTerminalGroups=True))


def eligible_ligand(name, openfe, old, measurements, sources, alignment):
    old_name = "lig_" + name
    if name not in openfe or old_name not in old or old_name not in measurements:
        return dict(name=name, reason="unmatched_source_name")
    a, b = openfe[name], old[old_name]
    record = measurements[old_name]
    assay = record["measurement"]
    if (assay["doi"] != sources["config"]["measurement_doi"] or assay["type"] != "ki" or assay["unit"] != "uM"
            or not assay["comment"].startswith(sources["config"]["measurement_locator"])):
        return dict(name=name, reason="unmatched_Ki_conditions")
    smiles = Chem.MolToSmiles
    if smiles(a) != smiles(b) or smiles(a) != smiles(Chem.MolFromSmiles(record["smiles"])):
        return dict(name=name, reason="unmatched_stereochemical_molecular_graph")
    dock = sources["docked"] / old_name / "complete.json"
    if not dock.is_file():
        return dict(name=name, reason="missing_pose")
    source_pool = Path(json.loads(dock.read_text())["poses"])
    try:
        model, pose = native_pose(source_pool)
        if smiles(pose) != smiles(a):
            raise ValueError("native pose chemical graph differs")
        rmsd = fixed_pose_rmsd(pose, a, alignment[0], alignment[1])
    except (ValueError, RuntimeError, OSError) as error:
        return dict(name=name, reason="native_pose_invalid", detail=str(error))
    return dict(name=name, old_name=old_name, native_model=model, rmsd_angstrom=rmsd, pool=str(source_pool),
                charge=Chem.GetFormalCharge(a), smiles=smiles(a), id_key=Chem.MolToInchiKey(a),
                reason="eligible" if rmsd < sources["config"]["maximum_fixed_pose_rmsd_angstrom"] else "pose_RMSD")


def analyse_mcs(a, b):
    result = rdFMCS.FindMCS([a, b], timeout=5, ringMatchesRingOnly=True, completeRingsOnly=True,
                           matchValences=False, bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    if result.canceled or not result.numAtoms:
        raise ValueError("No resolved heavy-atom MCS")
    counts = [m.GetNumHeavyAtoms() for m in (a, b)]
    return (sum(counts) - 2*result.numAtoms, result.numAtoms/max(counts), result.smartsString)


def load_rows(sources):
    with bz2.open(sources["archive"], "rt") as stream:
        raw = json.load(stream)
    conf = sources["config"]
    return [row for row in raw["ddg"] if row["system_name"] == conf["target"]
            and row["system_group"] == conf["system_group"]]


def freeze(sources, output):
    config = sources["config"]
    started = perf_counter()
    alignment = receptor_alignment(sources, config)
    openfe = ligand_molecules(sources["openfe_ligands"])
    old = ligand_molecules(sources["historical_ligands"])
    measurements = yaml.safe_load(sources["measurements"].read_text())
    edges = {(row["ligand_a"], row["ligand_b"]) for row in load_rows(sources)}
    ligands = {name: eligible_ligand(name, openfe, old, measurements, sources, alignment)
               for name in sorted({x for edge in edges for x in edge})}
    records = []
    for a, b in sorted(edges):
        first, second = ligands[a], ligands[b]
        record = dict(a=a, b=b, status="excluded")
        if first["reason"] != "eligible" or second["reason"] != "eligible":
            record["reason"] = f"A {first['reason']}; B {second['reason']}"
        elif first["charge"] != second["charge"]:
            record["reason"] = "formal_charge_change"
        else:
            try:
                edits, fraction, smarts = analyse_mcs(openfe[a], openfe[b])
                record.update(heavy_edits=edits, common_heavy_fraction=fraction, common_heavy_smarts=smarts,
                              max_pose_rmsd=max(first["rmsd_angstrom"], second["rmsd_angstrom"]))
                if not 0 < edits <= config["maximum_heavy_atom_edits"] or fraction < config["minimum_common_heavy_fraction"]:
                    record["reason"] = "heavy_edit_or_core_threshold"
                else:
                    record.update(status="candidate", reason="")
            except ValueError as error:
                record["reason"] = str(error)
        records.append(record)
    used = set()
    selected = []
    for row in sorted((r for r in records if r["status"] == "candidate"),
                      key=lambda r: (r["heavy_edits"], r["max_pose_rmsd"], r["a"], r["b"])):
        if used.intersection((row["a"], row["b"])):
            row["reason"] = "shares_selected_endpoint"
        elif len(selected) >= config["edge_count"]:
            row["reason"] = "after_predefined_budget"
        else:
            row["status"] = "selected"
            selected.append(row)
            used.update((row["a"], row["b"]))
    with output.open("x") as stream:
        json.dump(dict(frozen_utc=now(), wall_seconds=perf_counter()-started,
                       configuration_sha256=digest(sources["config_path"]), archive_sha256=digest(sources["archive"]),
                       alignment=alignment[2], selected=selected, eligibility=records, ligands=ligands),
                  stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(f"{len(records)} named OpenFE edges screened; {len(selected)} frozen. "
          f"Pairs: {[(r['a'], r['b']) for r in selected]}", flush=True)


def kcal(q):
    if q["unit"] not in ("kilocalorie / mole", "kilocalories_per_mole", "kilocalorie_per_mole"):
        raise ValueError(f"Unexpected archived energy units: {q['unit']}")
    x = float(q["magnitude"])
    if not isfinite(x):
        raise ValueError("Nonfinite archived energy")
    return x


def compute_legs(row, config):
    legs = {}
    replicates = []
    for leg in ("complex", "solvent"):
        values = list(map(kcal, row[f"dgs_{leg}"]))
        mbars = list(map(kcal, row[f"{leg}_mbar_errors"]))
        overlaps = row[f"{leg}_smallest_mbar_overlaps"]
        mixings = row[f"{leg}_smallest_replica_mixing"]
        equil = row[f"{leg}_equilibration_iterations"]
        prod = row[f"{leg}_production_iterations"]
        for replica, (value, sigma, overlap, mixing, equil_iterations, prod_iterations) in enumerate(
                zip(values, mbars, overlaps, mixings, equil, prod, strict=True), 1):
            if sigma < 0 or not 0 <= overlap <= 1 or not 0 <= mixing <= 1 or prod_iterations <= 0:
                raise ValueError(f"Invalid {leg}/{replica} diagnostics")
            replicates.append(dict(leg=leg, replica=replica, delta_kcal_mol=value, mbar_sd=sigma,
                                   minimum_mbar_overlap=overlap, minimum_replica_mixing=mixing,
                                   equilibration_iterations=equil_iterations, production_iterations=prod_iterations))
        n = len(values)
        if n < 2:
            raise ValueError(f"Fewer than two archived {leg} repeats")
        legs[leg] = dict(n=n, mean=mean(values), replica_sd=stdev(values), replica_sem=stdev(values)/sqrt(n),
                         within_replicate_mbar_sem=sqrt(sum(v*v for v in mbars))/n,
                         min_overlap=min(overlaps), min_mixing=min(mixings), mean_mbar_sd=mean(mbars))
    c, s = legs["complex"], legs["solvent"]
    delta = c["mean"] - s["mean"]
    if abs(delta - kcal(row["ddg"])) > 1e-6:
        raise ValueError("Recomputed complex–solvent difference does not match parent result")
    return dict(delta=delta, published_uncertainty=kcal(row["ddg_uncertainty"]),
                replica_sem=hypot(c["replica_sem"], s["replica_sem"]),
                within_replicate_mbar_sem=hypot(c["within_replicate_mbar_sem"], s["within_replicate_mbar_sem"]),
                min_overlap=min(c["min_overlap"], s["min_overlap"]),
                min_mixing=min(c["min_mixing"], s["min_mixing"]),
                quality_flags=[flag for flag, value, limit in (
                    ("low_overlap", min(c["min_overlap"], s["min_overlap"]), config["minimum_overlap_for_flag"]),
                    ("low_mixing", min(c["min_mixing"], s["min_mixing"]), config["minimum_mixing_for_flag"]))
                               if value < limit], legs=legs), replicates


def compare(sources, selection_path, output):
    started = perf_counter()
    config = sources["config"]
    frozen = json.loads(selection_path.read_text())
    if frozen["configuration_sha256"] != digest(sources["config_path"]) or frozen["archive_sha256"] != digest(sources["archive"]):
        raise ValueError("Frozen protocol or archive changed before analysis")
    raw_rows = {(row["ligand_a"], row["ligand_b"]): row for row in load_rows(sources)}
    measurements = yaml.safe_load(sources["measurements"].read_text())
    with sources["fast_scores"].open() as csvfile:
        fast_scores = {r["name"]: r for r in csv.DictReader(csvfile) if r["target"] == config["target"]}
    comparisons = []
    repeats = []
    archived_inputs = []
    for edge in frozen["selected"]:
        a, b = edge["a"], edge["b"]
        row = dict(edge)
        try:
            names = [frozen["ligands"][name]["old_name"] for name in (a, b)]
            ki = [float(measurements[name]["measurement"]["value"]) * 1000 for name in names]
            if not all(isfinite(x) and x > 0 for x in ki):
                raise ValueError("Missing valid exact positive measurement")
            row["A_Ki_nM"], row["B_Ki_nM"] = ki
            row["experimental_delta_kcal_mol"] = R_KCAL * config["temperature_kelvin"] * log(ki[1]/ki[0])
            for arm in SCORES:
                row[f"{arm}_B_minus_A"] = float(fast_scores[names[1]][arm]) - float(fast_scores[names[0]][arm])
            source_row = raw_rows[a, b]
            archived_inputs.append(source_row)
            computed, details = compute_legs(source_row, config)
            row.update(computed, status="complete")
            row["RBFE_error"] = computed["delta"] - row["experimental_delta_kcal_mol"]
            row["direction_correct"] = computed["delta"] * row["experimental_delta_kcal_mol"] > 0
            repeats.extend(dict(a=a, b=b, **r) for r in details)
        except (KeyError, ValueError, OSError) as error:
            row.update(status="failed", reason=str(error))
        comparisons.append(row)
    output.mkdir()
    put_json(output / "archived_edges.json", archived_inputs)
    put_csv(output / "by_repeat.csv", repeats)
    put_csv(output / "by_edge.csv", [{k: v for k, v in row.items() if k != "legs"} for row in comparisons])
    put_json(output / "summary.json", dict(frozen_utc=frozen["frozen_utc"], analysis_utc=now(),
        local_analysis_wall_seconds=perf_counter()-started, planned=config["edge_count"],
        selected=len(comparisons), completed=sum(r["status"] == "complete" for r in comparisons),
        edges=comparisons))
    print([(r["a"], r["b"], r["status"]) for r in comparisons], flush=True)


def prepare_openfe_inputs(sources, selection_path, output):
    """Preserve source SDF charges and poses for ordinary OpenFE planning."""
    frozen = json.loads(selection_path.read_text())
    if frozen["configuration_sha256"] != digest(sources["config_path"]):
        raise ValueError("Frozen configuration changed")
    source_sdf = sources["openfe_ligands"]
    payload = source_sdf.read_text()
    blocks = payload.split("$$$$\n")
    sdf_records = {block.splitlines()[0]: block for block in blocks if block.strip()}
    molecules = ligand_molecules(source_sdf)
    pairs = frozen["selected"]
    if not pairs:
        raise ValueError("There are no frozen edges")
    output.mkdir()
    for edge in pairs:
        a, b = edge["a"], edge["b"]
        owner = output / f"{a}_to_{b}"
        owner.mkdir()
        for name in (a, b):
            if Chem.MolToInchiKey(molecules[name]) != frozen["ligands"][name]["id_key"]:
                raise ValueError(f"{name}: OpenFE input no longer matches frozen identity")
        with (owner / "ligands.sdf").open("x") as stream:
            for name in (a, b):
                stream.write(sdf_records[name] + "$$$$\n")
        if set(ligand_molecules(owner / "ligands.sdf")) != {a, b}:
            raise ValueError("Extracted SDF names differ from frozen pair")
        with (owner / "protein.pdb").open("xb") as stream:
            stream.write(sources["openfe_protein"].read_bytes())
    print(f"Prepared {len(pairs)} exact coordinate/charge bundles in {output}", flush=True)


def plan_openfe(executable, pair_dir, settings, repeats):
    """Build OpenFE's two Transformation JSONs for an isolated two-molecule edge."""
    executable, env = openfe_cli_environment(executable)
    pair_dir = pair_dir.resolve(strict=True)
    for name in ("ligands.sdf", "protein.pdb"):
        if not (pair_dir / name).is_file():
            raise FileNotFoundError(pair_dir / name)
    plan = pair_dir / "network"
    if plan.exists():
        raise FileExistsError(plan)
    command = [str(executable), "plan-rbfe-network", "-M", str(pair_dir / "ligands.sdf"),
               "-p", str(pair_dir / "protein.pdb"), "-o", str(plan),
               "--n-protocol-repeats", str(repeats)]
    if settings is not None:
        command.extend(("-s", str(settings.resolve(strict=True))))
    subprocess.run(command, check=True, env=env)
    transformations = list((plan / "transformations").glob("*.json"))
    if len(transformations) != 2 or any(sum(f"_{leg}_" in p.stem for p in transformations) != 1
                                       for leg in ("complex", "solvent")):
        raise ValueError("Expected one complex and one solvent OpenFE Transformation")
    for path in transformations:
        print(path, flush=True)


def plan_openfe_hybrid(pair_dir, config_path, output, platform):
    """Build directed hybrid-topology inputs with the official OpenFE Python API."""
    import openfe
    from kartograf import KartografAtomMapper
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from openff.units import unit

    pair_dir = pair_dir.resolve(strict=True)
    profile = json.loads(config_path.read_text())
    source = Chem.SDMolSupplier(str(pair_dir / "ligands.sdf"), removeHs=False)
    a_mol, b_mol = list(source)
    if a_mol is None or b_mol is None or Chem.GetFormalCharge(a_mol) != Chem.GetFormalCharge(b_mol):
        raise ValueError("Prepare two valid ligands of the same formal charge")
    a = openfe.SmallMoleculeComponent.from_rdkit(a_mol)
    b = openfe.SmallMoleculeComponent.from_rdkit(b_mol)
    if a.name == b.name:
        raise ValueError("Two distinguishable inputs are required")
    for component in (a, b):
        molecule = component.to_rdkit()
        charges = [atom.GetDoubleProp("PartialCharge") for atom in molecule.GetAtoms()
                   if atom.HasProp("PartialCharge")]
        if len(charges) != molecule.GetNumAtoms() or not all(map(isfinite, charges)) or (
            abs(sum(charges) - Chem.GetFormalCharge(molecule)) > 1e-4
        ):
            raise ValueError(f"{component.name}: preassigned atomic charges do not match formal charge")
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=True)
    mapping = next(mapper.suggest_mappings(a, b))
    atom_map = mapping.componentA_to_componentB
    if any((a_mol.GetAtomWithIdx(i).GetAtomicNum() == 1) != (b_mol.GetAtomWithIdx(j).GetAtomicNum() == 1)
           for i, j in atom_map.items()):
        raise ValueError("OpenFE mapping contains a hydrogen-to-heavy-atom transition")
    unmatched = {
        "A": {i: atom.GetSymbol() for i, atom in enumerate(a_mol.GetAtoms()) if i not in atom_map},
        "B": {j: atom.GetSymbol() for j, atom in enumerate(b_mol.GetAtoms()) if j not in atom_map.values()},
    }
    protein = openfe.ProteinComponent.from_pdb_file(str(pair_dir / "protein.pdb"))
    solvent = openfe.SolventComponent(ion_concentration=profile["salt_molar"] * unit.molar)
    output.mkdir()
    metadata = {}
    for leg, components in (("complex", {"protein": protein}), ("solvent", {})):
        settings = RelativeHybridTopologyProtocol.default_settings()
        settings.protocol_repeats = profile["repeats"]
        settings.simulation_settings.equilibration_length = profile["equilibration_ps"] * unit.picosecond
        settings.simulation_settings.production_length = profile["production_ps"] * unit.picosecond
        settings.simulation_settings.time_per_iteration = profile["iteration_ps"] * unit.picosecond
        settings.simulation_settings.minimization_steps = profile["minimization_steps"]
        if "timestep_fs" in profile:
            settings.integrator_settings.timestep = profile["timestep_fs"] * unit.femtosecond
        settings.output_settings.checkpoint_interval = profile["iteration_ps"] * unit.picosecond
        settings.output_settings.positions_write_frequency = profile["iteration_ps"] * unit.picosecond
        settings.engine_settings.compute_platform = platform
        protocol = RelativeHybridTopologyProtocol(settings)
        state_a, state_b = (openfe.ChemicalSystem({**components, "solvent": solvent, "ligand": ligand},
                                                  name=f"{ligand.name}_{leg}") for ligand in (a, b))
        transformation = openfe.Transformation(stateA=state_a, stateB=state_b,
            mapping=mapping, protocol=protocol, name=f"{a.name}_to_{b.name}_{leg}")
        path = output / f"{leg}.json"
        transformation.to_json(file=path)
        metadata[leg] = str(path)
    put_json(output / "protocol.json", dict(openfe_version=openfe.__version__, profile=profile, platform=platform,
        pair=[a.name, b.name], atom_map=atom_map, unmatched_atoms=unmatched, transformations=metadata))


def openfe_cli_environment(executable):
    path = shutil.which(str(executable))
    if path is None:
        raise FileNotFoundError(executable)
    executable = Path(path).absolute()
    env = os.environ.copy()
    env["PATH"] = str(executable.parent) + os.pathsep + env.get("PATH", "")
    if shutil.which("openfe_analysis", path=str(executable.parent)) is None:
        raise FileNotFoundError(f"openfe_analysis is not installed beside {executable}")
    return executable, env


def analyse_openfe_leg(attempt, output):
    """Rerun the official energy and structure analyses on preserved simulation files."""
    from gufe.tokenization import JSON_HANDLER
    from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocolUnit
    from openfe.protocols.openmm_utils.multistate_analysis import MultistateEquilFEAnalysis
    from openff.units import unit
    from openmmtools.multistate import MultiStateReporter

    import sys
    attempt = attempt.resolve(strict=True)
    executable, env = openfe_cli_environment(Path(sys.executable).with_name("openfe"))
    os.environ["PATH"] = env["PATH"]
    files = [attempt / name for name in ("simulation.nc", "checkpoint.chk", "hybrid_system.pdb", "db.json")]
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    output.mkdir()
    output = output.resolve()
    for path in files:
        (output / path.name).symlink_to(path)
    reporter = MultiStateReporter(str(output / "simulation.nc"), open_mode="r", checkpoint_storage="checkpoint.chk")
    start = perf_counter()
    try:
        mbar = MultistateEquilFEAnalysis(reporter, sampling_method="repex", result_units=unit.kilocalorie_per_mole)
        with (output / "mbar.json").open("x") as stream:
            json.dump(mbar.unit_results_dict, stream, cls=JSON_HANDLER.encoder, indent=2)
        mbar.plot(filepath=output, filename_prefix="")
        mbar.close()
    finally:
        reporter.close()
    structural = RelativeHybridTopologyProtocolUnit.structural_analysis(output, output)
    with (output / "analysis_execution.json").open("x") as stream:
        json.dump(dict(source_attempt=str(attempt), openfe=str(executable),
                       structural_results={key: value.decode(errors="replace") if isinstance(value, bytes) else str(value)
                                           for key, value in structural.items()},
                       analysis_wall_seconds=perf_counter()-start),
                  stream, indent=2)
    if "structural_analysis_error" in structural:
        raise RuntimeError(structural["structural_analysis_error"].decode(errors="replace"))


def run_openfe(executable, transformation, output):
    """Run an individually scheduled leg with its own checkpoint and result directory."""
    executable, env = openfe_cli_environment(executable)
    transformation = transformation.resolve(strict=True)
    output = output.resolve()
    output.mkdir()
    command = [str(executable), "quickrun", str(transformation), "-d", str(output),
               "-o", str(output / "openfe_result.json")]
    put_json(output / "command.json", command)
    start = perf_counter()
    with (output / "quickrun.log").open("x") as log_file:
        process = subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT, env=env)
    put_json(output / "execution.json", dict(returncode=process.returncode, wall_seconds=perf_counter()-start))
    process.check_returncode()


def configuration(path):
    root = path.resolve().parent
    config = json.loads(path.read_text())
    prior = root / config["historical_run"]
    diagnostic = root / config["diagnostic_run"]
    target = config["target"]
    return dict(config=config, config_path=path, archive=root/config["archive"],
                openfe_ligands=root/config["openfe_ligands"], openfe_protein=root/config["openfe_protein"],
                historical_ligands=prior/f"data/source/{target}/ligands.sdf",
                historical_protein=prior/f"data/source/{target}/protein.pdb",
                measurements=prior/f"data/source/{target}/ligands.yml",
                docked=prior/f"data/docked/polar/{target}", fast_scores=diagnostic/"replay/ki_by_ligand.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    children = parser.add_subparsers(dest="stage", required=True)
    for stage in ("freeze", "compare", "prepare-openfe"):
        child = children.add_parser(stage)
        child.add_argument("--config", type=Path, required=True)
        child.add_argument("--selection", type=Path, required=True)
        if stage != "freeze":
            child.add_argument("--output", type=Path, required=True)
    plan = children.add_parser("plan-openfe")
    plan.add_argument("--openfe", default="openfe")
    plan.add_argument("--pair-dir", type=Path, required=True)
    plan.add_argument("--settings", type=Path)
    plan.add_argument("--repeats", type=int, default=3, choices=range(1, 9))
    hybrid = children.add_parser("plan-openfe-hybrid")
    hybrid.add_argument("--pair-dir", type=Path, required=True)
    hybrid.add_argument("--config", type=Path, required=True)
    hybrid.add_argument("--output", type=Path, required=True)
    hybrid.add_argument("--platform", choices=("cuda", "cpu"), required=True)
    run = children.add_parser("run-openfe")
    run.add_argument("--openfe", default="openfe")
    run.add_argument("--transformation", type=Path, required=True)
    run.add_argument("--output", type=Path, required=True)
    analysis = children.add_parser("analyse-openfe-leg")
    analysis.add_argument("--attempt", type=Path, required=True)
    analysis.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "analyse-openfe-leg":
        analyse_openfe_leg(args.attempt, args.output)
    elif args.stage == "plan-openfe":
        plan_openfe(args.openfe, args.pair_dir, args.settings, args.repeats)
    elif args.stage == "plan-openfe-hybrid":
        plan_openfe_hybrid(args.pair_dir, args.config, args.output, args.platform)
    elif args.stage == "run-openfe":
        run_openfe(args.openfe, args.transformation, args.output)
    else:
        sources = configuration(args.config)
        if args.stage == "freeze":
            freeze(sources, args.selection)
        elif args.stage == "compare":
            compare(sources, args.selection, args.output)
        else:
            prepare_openfe_inputs(sources, args.selection, args.output)


if __name__ == "__main__":
    main()
