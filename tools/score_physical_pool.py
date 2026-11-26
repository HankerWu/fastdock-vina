"""Score fixed Vina poses with optional structure-supported waters and side chains."""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import json
from math import isfinite, log
from pathlib import Path
import re
from shutil import copyfileobj
import subprocess
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem

ARMS = ("native", "direction", "burial", "combined", "water", "strain", "sidechain", "joint")
SATISFACTION_ARMS = ("satisfied", "satisfied_water")
SATISFACTION_FIELDS = ("direct_credit_area", "additional_water_credit_area",
                       "direct_credit", "additional_water_credit",
                       "satisfied_delta", "satisfied_water_delta")
STRAIN_WEIGHT = 0.1
MMFF_VARIANT = "MMFF94s"
GAS_CONFORMERS = 8
GAS_EMBED_SEED = 777
RT_KCAL = 0.592
CNN_CONSENSUS_WEIGHT = 0.592
CNN_ENDPOINT_EPSILON = 1e-6
CLASH_TOLERANCE_ANGSTROM = 0.4
CLASH_WEIGHT = 2.0
ATOM_RADII = {
    "C": 1.70, "A": 1.70, "N": 1.55, "NA": 1.55, "OA": 1.52, "O": 1.52,
    "S": 1.80, "SA": 1.80, "P": 1.80,
}
RECEPTOR_HYDROGEN = {"H", "HD", "HS"}
PDBQT_ELEMENT_ALIASES = {"A": "C", "NA": "N", "OA": "O", "SA": "S",
                         **{f"CG{i}": "C" for i in range(4)}}
PDBQT_DUMMY_TYPES = {"W", "G0", "G1", "G2", "G3"}
AROMATIC_SCORE_WINDOW = 0.5
STACK_DISTANCE_MIN = 0.5
STACK_DISTANCE_MAX = 5.5
STACK_ANGLE_DEVIATION = 30.0
STACK_OFFSET_MAX = 2.0
CCD_MAPPING_LIMIT = 4096
# Broad covalent aromatic C-C/C-N interval and planarity guard used to
# reject residues assembled from corrupted PDBQT identities.
RECEPTOR_RING_BOND_RANGE = (1.2, 1.6)
AROMATIC_RING_PLANE_RMS_MAX = 0.2
PHENYL_RING = ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")
HISTIDINE_RING = ("CG", "ND1", "CE1", "NE2", "CD2")
PROTEIN_RINGS = {
    "PHE": (PHENYL_RING,), "TYR": (PHENYL_RING,),
    "HIS": (HISTIDINE_RING,), "HID": (HISTIDINE_RING,), "HIE": (HISTIDINE_RING,),
    "HIP": (HISTIDINE_RING,), "HSD": (HISTIDINE_RING,), "HSE": (HISTIDINE_RING,),
    "HSP": (HISTIDINE_RING,),
    "TRP": (("CG", "CD1", "NE1", "CE2", "CD2"), ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")),
}


@dataclass(frozen=True)
class Pose:
    model: int
    native_score: float
    native_rank: int
    text: str
    native_intra: float | None
    native_unbound: float | None


def read_poses(path):
    """Read original affinity/rank and one unwrapped PDBQT ligand per MODEL."""
    lines = path.read_text().splitlines()
    poses = []
    model = None
    block = []
    for line in lines:
        if line.startswith("MODEL "):
            if model is not None:
                raise ValueError(f"Nested MODEL in {path}")
            model, block = int(line.split()[1]), []
        elif line.strip() == "ENDMDL":
            if model is None:
                raise ValueError(f"Orphan ENDMDL in {path}")
            native = [l.split()[3:] for l in block if l.startswith("REMARK VINA NATIVE:")]
            if native:
                [[score, rank]] = native
                rank = int(rank)
            else:
                [score] = [l.split()[3] for l in block if l.startswith("REMARK VINA RESULT:")]
                rank = len(poses) + 1
            if not isfinite(float(score)) or not any(l.startswith(("ATOM ", "HETATM ")) for l in block):
                raise ValueError(f"Invalid native score or atom records in {path}, MODEL {model}")
            body = [l for l in block if not l.startswith(("REMARK VINA ", "REMARK INTER",
                                                           "REMARK INTRA", "REMARK UNBOUND"))]
            def saved_energy(label):
                record = [float(line.split()[-1]) for line in block if line.startswith(label)]
                return record[0] if len(record) == 1 and isfinite(record[0]) else None

            poses.append(Pose(
                model, float(score), rank, "\n".join(body) + "\n",
                saved_energy("REMARK INTRA:"), saved_energy("REMARK UNBOUND:"),
            ))
            model = None
        elif model is not None:
            block.append(line)
        elif line.strip():
            raise ValueError(f"Unexpected line outside MODEL in {path}")
    if model is not None or not poses:
        raise ValueError(f"Missing complete MODEL records in {path}")
    if (
        [p.model for p in poses] != list(range(1, len(poses) + 1))
        or sorted(p.native_rank for p in poses) != list(range(1, len(poses) + 1))
    ):
        raise ValueError(f"Invalid model/rank identities in {path}")
    return poses


def pdbqt_chemistry_record(text):
    """Return the ordered ligand representation, excluding atom coordinates."""
    return tuple(line[:30] + line[54:] if line.startswith(("ATOM ", "HETATM ")) else line
                 for line in text.splitlines()
                 if line.startswith(("ATOM ", "HETATM ", "ROOT", "ENDROOT", "BRANCH ", "ENDBRANCH ",
                                     "TORSDOF ", "REMARK SMILES ", "REMARK H PARENT")))


def serial_index_pairs(lines, label):
    tokens = [token for line in lines if line.startswith(label) for token in line[len(label):].split()]
    if any(not token.isascii() or not token.isdecimal() or int(token) < 1 for token in tokens):
        raise ValueError(f"{label.strip()} requires positive decimal indices")
    indices = map(int, tokens)
    return list(zip(indices, indices, strict=True))


def mapped_ligand_bonds(text):
    """Recover the full atom adjacency represented by one prepared Meeko PDBQT."""
    lines = text.splitlines()
    if any(line.startswith(("MODEL ", "ENDMDL")) for line in lines):
        raise ValueError("Ligand bond mapping requires one unwrapped prepared PDBQT")
    [smiles] = [line.removeprefix("REMARK SMILES ") for line in lines
                if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX ")]
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError("Prepared Meeko SMILES must encode one connected ligand")
    atoms = {}
    for line in lines:
        if not line.startswith(("ATOM ", "HETATM ")):
            continue
        serial, atom_type = int(line[6:11]), line.split()[-1]
        if serial < 1 or serial in atoms:
            raise ValueError("Prepared PDBQT requires unique positive atom serials")
        if atom_type in PDBQT_DUMMY_TYPES or atom_type.startswith("CG"):
            raise ValueError("Ligand bond mapping does not support virtual or macrocycle atoms")
        symbol = ("H" if atom_type in RECEPTOR_HYDROGEN
                  else PDBQT_ELEMENT_ALIASES.get(atom_type, atom_type.upper()))
        position = np.array([float(line[i:i + 8]) for i in (30, 38, 46)])
        if not np.isfinite(position).all():
            raise ValueError("Ligand bond mapping requires finite prepared coordinates")
        atoms[serial] = symbol, position
    encoded = serial_index_pairs(lines, "REMARK SMILES IDX ")
    index_to_serial = {index: serial for index, serial in encoded}
    if (len(index_to_serial) != len(encoded) or len({serial for _, serial in encoded}) != len(encoded)
            or set(index_to_serial) != set(range(1, molecule.GetNumAtoms() + 1))):
        raise ValueError("REMARK SMILES IDX must map each encoded atom one-to-one")
    atom_numbers = {}
    for atom in molecule.GetAtoms():
        serial = index_to_serial[atom.GetIdx() + 1]
        if serial not in atoms or atoms[serial][0] != atom.GetSymbol().upper():
            raise ValueError(f"REMARK SMILES IDX element mismatch for PDBQT serial {serial}")
        atom_numbers[serial] = atom.GetAtomicNum()
    edges = {tuple(sorted((index_to_serial[bond.GetBeginAtomIdx()+1], index_to_serial[bond.GetEndAtomIdx()+1])))
             for bond in molecule.GetBonds()}
    marked = Counter()
    for parent, serial in serial_index_pairs(lines, "REMARK H PARENT"):
        if parent not in index_to_serial or serial in atom_numbers or atoms.get(serial, (None,))[0] != "H":
            raise ValueError(f"Invalid explicit hydrogen mapping for PDBQT serial {serial}")
        attached_atom = molecule.GetAtomWithIdx(parent - 1)
        marked[parent] += 1
        if attached_atom.GetAtomicNum() == 1 or marked[parent] > attached_atom.GetTotalNumHs():
            raise ValueError(f"SMILES cannot justify explicit hydrogen parent at index {parent}")
        atom_numbers[serial] = 1
        edges.add(tuple(sorted((index_to_serial[parent], serial))))
    if set(atom_numbers) != atoms.keys():
        raise ValueError("Meeko remarks do not account for every prepared PDBQT atom")
    if any(tuple(sorted((left, right))) not in edges for left, right in serial_index_pairs(lines, "BRANCH ")):
        raise ValueError("PDBQT BRANCH edge is absent from the reconstructed Meeko ligand graph")
    periodic_table = Chem.GetPeriodicTable()
    for left, right in sorted(edges):
        expected = sum(periodic_table.GetRcovalent(atom_numbers[serial]) for serial in (left, right))
        distance = np.linalg.norm(atoms[left][1] - atoms[right][1])
        if not expected - 0.4 <= distance <= expected + 0.35:
            raise ValueError(f"Prepared bond {left}–{right} has implausible length {distance:.3f} Å")
    return sorted(edges)


def validate_ligand_bonds(path, expected):
    edges = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        if len(tokens) != 2 or any(not token.isascii() or not token.isdecimal() for token in tokens):
            raise ValueError("Ligand bonds file requires two positive serials per line")
        edges.append(tuple(sorted(map(int, tokens))))
    if sorted(edges) != expected:
        raise ValueError("Ligand bonds file differs from the prepared Meeko graph")


def fixed_coordinate_inputs(poses, prepared_ligand, ligand_bonds):
    coordinate_inputs = {}
    if ligand_bonds and not prepared_ligand:
        raise ValueError("Ligand bonds require a prepared ligand")
    if prepared_ligand:
        prepared_text = prepared_ligand.read_text()
        if any(line.startswith(("MODEL ", "ENDMDL")) for line in prepared_text.splitlines()):
            raise ValueError("Prepared ligand must be one unwrapped PDBQT")
        chemistry = pdbqt_chemistry_record(prepared_text)
        if not any(line.startswith(("ATOM ", "HETATM ")) for line in chemistry):
            raise ValueError("Prepared ligand has no PDBQT atoms")
        for pose in poses:
            if pdbqt_chemistry_record(pose.text) != chemistry:
                raise ValueError(f"MODEL {pose.model} differs from the prepared PDBQT outside coordinates")
        coordinate_inputs["prepared_ligand"] = prepared_ligand
        if ligand_bonds:
            validate_ligand_bonds(ligand_bonds, mapped_ligand_bonds(prepared_text))
            coordinate_inputs["ligand_bonds"] = ligand_bonds
    return coordinate_inputs


def fields(stdout, label, count):
    [values] = [line[len(label):].split() for line in stdout.splitlines() if line.startswith(label)]
    numbers = tuple(map(float, values))
    if len(numbers) != count or not all(map(isfinite, numbers)):
        raise ValueError(f"Malformed engine record: {label}")
    return numbers


def evaluate_coordinates(binary, ligand, receptor, waters, timeout, *, satisfaction=False,
                         prepared_ligand=None, ligand_bonds=None):
    mode = "combined"
    if satisfaction:
        mode = "satisfied_water" if waters else "satisfied"
    command = [str(binary), "--receptor", str(receptor), "--ligand", str(ligand),
               "--score_only", "--net_polar", mode, "--cpu", "1", "--seed", "42",
               "--verbosity", "0"]
    if prepared_ligand:
        command[command.index("--ligand") + 1] = str(prepared_ligand)
        command.extend(("--pose_coordinates", str(ligand)))
    if ligand_bonds:
        command.extend(("--ligand_bonds", str(ligand_bonds)))
    if waters:
        command += ["--net_polar_waters", str(waters)]
    started = perf_counter()
    proc = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    if proc.returncode:
        raise RuntimeError(f"Scoring failed for {receptor}: {proc.stderr or proc.stdout}")
    term_names = ("native_hbond", "directed_hbond", "direction_raw", "buried_area", "water_area_credit")
    values = fields(proc.stdout, "Net polar terms:", 5)
    direction, burial, combined, water = fields(proc.stdout, "Net polar deltas:", 4)
    inter, intra, zero_reference = fields(proc.stdout, "Net polar native:", 3)
    (scale,) = fields(proc.stdout, "Net polar normalization:", 1)
    (water_count,) = fields(proc.stdout, "Net polar water points:", 1)
    if not 0 < scale <= 1 or water_count != int(water_count) or bool(water_count) != bool(waters):
        raise ValueError("Unexpected flexibility scale or independent water point count")
    result = {**dict(zip(term_names, values, strict=True)), "inter": inter, "intra": intra,
              "zero_reference": zero_reference, "scale": scale, "direction_delta": direction,
              "burial_delta": burial, "combined_delta": combined, "water_delta": water,
              "water_count": int(water_count), "wall_seconds": perf_counter() - started}
    if satisfaction:
        values = fields(proc.stdout, "Net polar satisfaction:", len(SATISFACTION_FIELDS))
        result.update(zip(SATISFACTION_FIELDS, values, strict=True))
    return result


def decode_ligand(pose):
    """Use Meeko's PDBQT atom/SMILES map and preserve the prepared microstate."""
    [ligand] = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(pose.text, skip_typing=True))
    if ligand is None:
        raise ValueError(f"Meeko atom mapping missing for MODEL {pose.model}")
    [expected] = [line.split()[2] for line in pose.text.splitlines()
                  if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX ")]
    reference = Chem.MolFromSmiles(expected)
    if reference is None or Chem.MolToSmiles(Chem.RemoveHs(ligand)) != Chem.MolToSmiles(reference):
        raise ValueError(f"Decoded identity differs from prepared SMILES for MODEL {pose.model}")
    return ligand, Chem.MolToSmiles(reference)

def export_gnina_sdf(pool, destination):
    """Write one molecule per native MODEL with Meeko-mapped heavy coordinates."""
    poses = read_poses(pool)
    with TemporaryDirectory(prefix=".gnina-export-", dir=destination.parent) as work:
        staged = Path(work) / "poses.sdf"
        with Chem.SDWriter(str(staged)) as writer:
            for pose in poses:
                ligand, _ = decode_ligand(pose)
                ligand.SetProp("_Name", f"native_model_{pose.model}")
                ligand.SetIntProp("native_model", pose.model)
                writer.write(ligand)
        records = Chem.SDMolSupplier(str(staged), sanitize=False, removeHs=False)
        for pose, molecule in zip(poses, records, strict=True):
            check_heavy_atoms(molecule, pdbqt_heavy_atoms(pose), pose.model, "Meeko export")
        with staged.open("rb") as source:
            target = destination.open("xb")
            try:
                with target:
                    copyfileobj(source, target)
            except BaseException:
                destination.unlink()
                raise


def pdbqt_heavy_atoms(pose):
    atoms = Counter()
    for line in pose.text.splitlines():
        if line.startswith(("ATOM ", "HETATM ")):
            atom_type = line.split()[-1]
            if atom_type not in RECEPTOR_HYDROGEN and atom_type not in PDBQT_DUMMY_TYPES:
                xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
                element = PDBQT_ELEMENT_ALIASES.get(atom_type, atom_type.upper())
                atoms[element, xyz] += 1
    return atoms


def check_heavy_atoms(molecule, expected, model, source):
    if molecule is None:
        raise ValueError(f"{source} record could not be read for MODEL {model}")
    xyz = molecule.GetConformer().GetPositions()
    observed = Counter((atom.GetSymbol().upper(), tuple(xyz[atom.GetIdx()]))
                       for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1)
    if observed != expected:
        raise ValueError(f"{source} heavy-atom elements or coordinates differ for MODEL {model}")


def read_gnina_scores(poses, scores_sdf):
    """Read raw GNINA scores only after certifying one unchanged pose per MODEL."""
    expected = {pose.model: pdbqt_heavy_atoms(pose) for pose in poses}
    scores = {}
    for record, molecule in enumerate(Chem.SDMolSupplier(str(scores_sdf), sanitize=False, removeHs=False), start=1):
        if molecule is None or not molecule.HasProp("native_model") or not molecule.HasProp("CNNscore"):
            raise ValueError(f"GNINA SDF record {record} needs native_model and CNNscore")
        model = int(molecule.GetProp("native_model"))
        if model in scores or model not in expected:
            raise ValueError(f"Unexpected or duplicate GNINA native_model {model}")
        cnn = float(molecule.GetProp("CNNscore"))
        if not isfinite(cnn) or not 0 <= cnn <= 1:
            raise ValueError(f"GNINA CNNscore outside [0, 1] for MODEL {model}")
        check_heavy_atoms(molecule, expected[model], model, "GNINA")
        scores[model] = cnn
    if scores.keys() != expected.keys():
        raise ValueError(f"GNINA SDF missing native_model values {sorted(expected.keys() - scores.keys())}")
    return scores


def gnina_rank_terms(pose, cnn, *, native_score=None):
    c = min(1 - CNN_ENDPOINT_EPSILON, max(CNN_ENDPOINT_EPSILON, cnn))
    score = pose.native_score if native_score is None else native_score
    return {"gnina": cnn, "cnn_consensus": score - CNN_CONSENSUS_WEIGHT * log(c / (1 - c))}


def score_gnina_pool(pool, scores_sdf):
    poses = read_poses(pool)
    scores = read_gnina_scores(poses, scores_sdf)
    rows = [{"model": pose.model, "native": pose.native_score, "native_rank": pose.native_rank,
             **gnina_rank_terms(pose, scores[pose.model])} for pose in poses]
    return {"poses": rows, "orders": rank_orders(rows, ("native", "gnina", "cnn_consensus"))}


VINARDO_TERMS = {
    "inter": re.compile(r"^\(1\) Final Intermolecular Energy\s+:\s+([-\d.]+)", re.M),
    "intra": re.compile(r"^\(2\) Final Total Internal Energy\s+:\s+([-\d.]+)", re.M),
}


def add_vinardo_order(result, poses, binary, receptor, box, timeout=120, *,
                       prepared_ligand=None, ligand_bonds=None):
    """Rank identical ligand coordinates with explicit Vinardo inter plus intra."""
    fixed_coordinate_inputs(poses, prepared_ligand, ligand_bonds)
    started = perf_counter()
    rows = {row["model"]: row for row in result["poses"]}
    for row in rows.values():
        row.update(vinardo=None, vinardo_terms=None)
    failures = []
    with TemporaryDirectory(prefix="vina-fixed-vinardo-") as temporary:
        ligand = Path(temporary) / "pose.pdbqt"
        command = [
            str(binary.resolve()), "--receptor", str(receptor), "--ligand", str(ligand),
            "--config", str(box), "--cpu", "1", "--score_only", "--unbound_energy", "0",
            "--scoring", "vinardo",
        ]
        if prepared_ligand:
            command[command.index("--ligand") + 1] = str(prepared_ligand)
            command.extend(("--pose_coordinates", str(ligand)))
        if ligand_bonds:
            command.extend(("--ligand_bonds", str(ligand_bonds)))
        for pose in poses:
            timer = perf_counter()
            try:
                ligand.write_text(pose.text)
                stdout = subprocess.run(command, capture_output=True, text=True, check=True,
                                        timeout=timeout).stdout
                terms = {}
                for name, pattern in VINARDO_TERMS.items():
                    [value] = pattern.findall(stdout)
                    terms[name] = Decimal(value)
                total = float(sum(terms.values()))
                if not isfinite(total):
                    raise ValueError("Non-finite Vinardo energy")
                rows[pose.model].update(
                    vinardo=total, vinardo_terms={**{key: float(value) for key, value in terms.items()},
                                                 "wall_seconds": perf_counter() - timer})
            except (OSError, ValueError, InvalidOperation, subprocess.SubprocessError) as error:
                reason = str(error)
                if isinstance(error, subprocess.CalledProcessError) and error.stderr:
                    reason = f"exit {error.returncode}: {error.stderr.strip()}"
                failures.append({"model": pose.model, "reason": reason,
                                 "wall_seconds": perf_counter() - timer})
    if failures:
        result.setdefault("abstentions", {})["vinardo"] = (
            f"Fixed-coordinate scoring failed for {len(failures)} of {len(poses)} models"
        )
        result["orders"]["vinardo"] = result["orders"]["native"].copy()
    else:
        result["orders"].update(rank_orders(result["poses"], ("vinardo",)))
    result["vinardo_scoring"] = {
        "binary": str(binary.resolve()), "receptor": str(receptor), "box": str(box),
        "reference": "Explicit intermolecular plus intramolecular energy; unbound energy fixed at zero",
        "model_count": len(poses), "scored_count": len(poses) - len(failures),
        "failures": failures, "wall_seconds": perf_counter() - started,
    }
    if prepared_ligand:
        result["vinardo_scoring"]["prepared_ligand"] = str(prepared_ligand)
        result["vinardo_scoring"]["ligand_bonds"] = str(ligand_bonds) if ligand_bonds else None
    return result


def read_ccd_graph(path):
    """Build an independent heavy-atom graph from CCD atom and bond tables."""
    from gemmi import cif

    block = cif.read_file(str(path)).sole_block()
    graph = Chem.RWMol()
    names = {}
    for name, symbol, charge in block.find("_chem_comp_atom.", ["atom_id", "type_symbol", "charge"]):
        name = cif.as_string(name)
        if symbol in {"H", "D", "T"}:
            continue
        if name in names:
            raise ValueError(f"Duplicate CCD atom name {name}")
        atom = Chem.Atom(symbol.capitalize())
        atom.SetFormalCharge(int(charge))
        atom.SetProp("ccd_atom_id", name)
        names[name] = graph.AddAtom(atom)
    if not names:
        raise ValueError("CCD has no heavy atoms")
    bond_types = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                  "TRIP": Chem.BondType.TRIPLE}
    neighbors = {name: 0 for name in names}
    for left, right, order in block.find("_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        left, right = map(cif.as_string, (left, right))
        if left in names and right in names:
            if order not in bond_types:
                raise ValueError(f"Unsupported CCD heavy bond order {order}")
            graph.AddBond(names[left], names[right], bond_types[order])
        elif left in names:
            neighbors[left] += 1
        elif right in names:
            neighbors[right] += 1
    for name, index in names.items():
        atom = graph.GetAtomWithIdx(index)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(neighbors[name])
    molecule = graph.GetMol()
    Chem.SanitizeMol(molecule)
    return molecule


def aromatic_cycles(molecule):
    return {tuple(sorted(ring)): tuple(ring) for ring in Chem.GetSymmSSSR(molecule)
            if all(molecule.GetAtomWithIdx(index).GetIsAromatic() for index in ring)
            and all(molecule.GetBondBetweenAtoms(first, second).GetIsAromatic()
                    for first, second in zip(ring, (*ring[1:], ring[0])))}


def ligand_ring_mapping(ligand, ccd):
    """Accept symmetry only when the entire aromatic-ring set is invariant."""
    # The prepared SMILES graph supplies the coordinate indexing. CCD matching
    # verifies connectivity, bond orders, charge and H count independently.
    if Chem.MolToSmiles(ligand, isomericSmiles=False) != Chem.MolToSmiles(ccd, isomericSmiles=False):
        raise ValueError("Prepared heavy graph or microstate differs from CCD")
    candidates = ligand.GetSubstructMatches(ccd, uniquify=False, maxMatches=CCD_MAPPING_LIMIT)
    if not candidates or len(candidates) == CCD_MAPPING_LIMIT:
        raise ValueError("CCD graph map absent or too symmetric to certify")
    rings = aromatic_cycles(ligand)
    ccd_cycles = aromatic_cycles(ccd)
    ccd_labels = defaultdict(set)
    for match in candidates:
        mapped = {tuple(sorted(match[index] for index in cycle)): cycle for cycle in ccd_cycles}
        if mapped.keys() != rings.keys():
            raise ValueError("CCD mapping changes aromatic ring membership")
        for ring, source in mapped.items():
            ccd_labels[ring].add(tuple(sorted(ccd.GetAtomWithIdx(index).GetProp("ccd_atom_id")
                                              for index in source)))
    return [(rings[key], sorted(ccd_labels[key])) for key in sorted(rings)], len(candidates)


def ring_plane(positions):
    positions = np.asarray(positions, dtype=float)
    if not np.isfinite(positions).all():
        raise ValueError("Nonfinite aromatic ring coordinates")
    center = positions.mean(axis=0)
    _, widths, axes = np.linalg.svd(positions - center, full_matrices=False)
    if widths[1] <= 1e-6:
        raise ValueError("Degenerate aromatic ring coordinates")
    return center, axes[-1], float(widths[-1] / np.sqrt(len(positions)))


def checked_aromatic_plane(coordinates, names, lower, upper):
    xyz = np.asarray(coordinates)
    # The input order is the chemical cycle, including the closing bond.
    bond_lengths = np.linalg.norm(xyz - np.roll(xyz, -1, axis=0), axis=1)
    low, high, _ = np.broadcast_arrays(lower, upper, bond_lengths)
    invalid = np.flatnonzero(~((low <= bond_lengths) & (bond_lengths <= high)))
    if invalid.size:
        i = invalid[0]
        following = (i + 1) % len(names)
        raise ValueError(f"Implausible aromatic bond {names[i]}–{names[following]} "
                         f"({bond_lengths[i]:.3f} Å, expected {low[i]:.2f}–{high[i]:.2f} Å)")
    plane = ring_plane(xyz)
    if plane[2] > AROMATIC_RING_PLANE_RMS_MAX:
        raise ValueError(f"Nonplanar aromatic ring ({plane[2]:.3f} Å plane RMS)")
    return plane


def receptor_ring_planes(path):
    residues = defaultdict(lambda: defaultdict(list))
    for identity, atom_type, xyz in receptor_atoms(path):
        atom_name = identity[0][6:].strip()
        residue = identity[1]
        if residue[:3] in PROTEIN_RINGS:
            residues[residue][atom_name].append((atom_type, xyz))
    rings, excluded = [], []
    for residue, atoms in sorted(residues.items()):
        for names in PROTEIN_RINGS[residue[:3]]:
            label = {"residue": residue.strip(), "atoms": names}
            try:
                coordinates = []
                for name in names:
                    entries = atoms.get(name, ())
                    if len(entries) != 1:
                        raise ValueError(f"Ring atom {name} is missing or duplicated")
                    atom_type, xyz = entries[0]
                    if PDBQT_ELEMENT_ALIASES.get(atom_type, atom_type.upper()) != name[0]:
                        raise ValueError(f"Ring atom {name} has inconsistent element")
                    coordinates.append(xyz)
                # Template order supplies the cyclic bonds. PDBQT has no
                # receptor bond table.
                plane = checked_aromatic_plane(coordinates, names, *RECEPTOR_RING_BOND_RANGE)
                rings.append((label, plane))
            except ValueError as error:
                excluded.append({**label, "reason": str(error)})
    return rings, excluded


def aromatic_geometry(ligand, receptor):
    ligand_center, ligand_normal, ligand_plane_rms = ligand
    receptor_center, receptor_normal, receptor_plane_rms = receptor
    delta = ligand_center - receptor_center
    distance = float(np.linalg.norm(delta))
    if distance >= STACK_DISTANCE_MAX:
        return None
    angle = float(np.degrees(np.arccos(np.clip(abs(ligand_normal @ receptor_normal), 0.0, 1.0))))
    ligand_offset = float(np.linalg.norm(delta - ligand_normal * (delta @ ligand_normal)))
    receptor_offset = float(np.linalg.norm(delta - receptor_normal * (delta @ receptor_normal)))
    stack = None
    if distance > STACK_DISTANCE_MIN and min(ligand_offset, receptor_offset) < STACK_OFFSET_MAX:
        if angle < STACK_ANGLE_DEVIATION:
            stack = "parallel"
        elif angle > 90.0 - STACK_ANGLE_DEVIATION:
            stack = "t_shaped"
    return {"kind": stack, "center_distance_angstrom": distance, "normal_angle_degrees": angle,
            "ligand_plane_offset_angstrom": ligand_offset,
            "receptor_plane_offset_angstrom": receptor_offset,
            "ligand_plane_rms_angstrom": ligand_plane_rms,
            "receptor_plane_rms_angstrom": receptor_plane_rms}


def add_aromatic_order(result, poses, receptor, ccd):
    """Attach a lexicographic geometry order with native fallback for the whole pool."""
    started = perf_counter()
    rows = {row["model"]: row for row in result["poses"]}
    native_order = result["orders"]["native"]
    scores = {pose.model: Decimal(str(pose.native_score)) for pose in poses}
    cutoff = min(scores.values()) + Decimal(str(AROMATIC_SCORE_WINDOW))
    eligible = [model for model in native_order if scores[model] <= cutoff]
    info = {"ccd": str(ccd), "near_native_max_delta_kcal_mol": AROMATIC_SCORE_WINDOW,
            "eligible_models": eligible, "receptor_rings": [], "excluded_receptor_rings": [],
            "ligand_rings": [], "ccd_isomorphisms": None}
    reports = {}
    remove_hydrogens = Chem.RemoveHsParameters()
    # Stereo-defining explicit H is part of the prepared identity check in
    # decode_ligand, but has no bearing on the heavy-ring correspondence.
    remove_hydrogens.removeDefiningBondStereo = True
    try:
        graph = read_ccd_graph(ccd)
        receptor_rings, excluded = receptor_ring_planes(receptor)
        info["receptor_rings"] = [label for label, _ in receptor_rings]
        info["excluded_receptor_rings"] = excluded
        if not receptor_rings:
            raise ValueError("No fully identified standard aromatic receptor ring")
        for pose in poses:
            try:
                ligand, _ = decode_ligand(pose)
                # Meeko reads PDBQT coordinates as float32. Match that
                # representation exactly, then use PDBQT's 0.001 Å precision.
                expected = Counter({(element, tuple(np.float32(xyz))): count
                                    for (element, xyz), count in pdbqt_heavy_atoms(pose).items()})
                check_heavy_atoms(ligand, expected, pose.model, "Meeko aromatic mapping")
                molecule = Chem.RemoveHs(ligand, remove_hydrogens)
                if pose is poses[0]:
                    rings, n_mappings = ligand_ring_mapping(molecule, graph)
                    if not rings:
                        raise ValueError("CCD/Meeko ligand has no aromatic cycles")
                    info["ligand_rings"] = [{"meeko_heavy_indices": sorted(ring), "ccd_ring_names": labels}
                                            for ring, labels in rings]
                    info["ccd_isomorphisms"] = n_mappings
                    periodic_table = Chem.GetPeriodicTable()
                    ring_bond_limits = []
                    for ring, _ in rings:
                        radii = np.array([periodic_table.GetRcovalent(molecule.GetAtomWithIdx(i).GetAtomicNum())
                                          for i in ring])
                        expected = radii + np.roll(radii, -1)
                        # Single-bond covalent radii bound the size range of
                        # heterogeneous aromatic bonds without a C/N-only gate.
                        ring_bond_limits.append((expected - 0.30, expected + 0.20))
                    indexed_graph = Chem.MolToSmiles(molecule, canonical=False)
                    # Canonical=False SMILES can still traverse equivalent
                    # atom orders. Match each prepared atom index explicitly.
                    atom_bonds = indexed_bonds(molecule)
                elif (Chem.MolToSmiles(molecule, canonical=False) != indexed_graph
                      or indexed_bonds(molecule) != atom_bonds):
                    raise ValueError("Prepared ligand indexing or chemistry differs across MODELs")
                coordinates = np.round(molecule.GetConformer().GetPositions(), 3)
                nearby = []
                for ring_index, ((ring, _), limits) in enumerate(zip(rings, ring_bond_limits, strict=True)):
                    ligand_plane = checked_aromatic_plane(coordinates[list(ring)], ring, *limits)
                    for receptor_index, (_, receptor_plane) in enumerate(receptor_rings):
                        geometry = aromatic_geometry(ligand_plane, receptor_plane)
                        if geometry is not None:
                            nearby.append({"ligand_ring": ring_index, "receptor_ring": receptor_index,
                                           **geometry})
                reports[pose.model] = nearby
            except (RuntimeError, ValueError) as error:
                raise ValueError(f"MODEL {pose.model}: {error}") from error
    except (OSError, RuntimeError, ValueError) as error:
        result.setdefault("abstentions", {})["aromatic"] = str(error)

    for model, row in rows.items():
        nearby = reports.get(model)
        row["aromatic_contacts"] = (
            None if nearby is None else sum(pair["kind"] is not None for pair in nearby)
        )
        row["aromatic_ring_pairs"] = nearby
    if len(reports) == len(poses):
        ranked = sorted(eligible, key=lambda model: -rows[model]["aromatic_contacts"])
        result["orders"]["aromatic"] = ranked + [model for model in native_order if model not in eligible]
    else:
        result["orders"]["aromatic"] = native_order.copy()
    info["diagnosed_models"] = list(reports)
    info["wall_seconds"] = perf_counter() - started
    result["aromatic_geometry"] = info
    return result


def indexed_bonds(molecule):
    atoms = tuple((atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(), atom.GetTotalNumHs())
                  for atom in molecule.GetAtoms())
    bonds = sorted((min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                    max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), bond.GetBondType())
                   for bond in molecule.GetBonds())
    return atoms, bonds


def score_aromatic_pool(pool, receptor, ccd):
    """Score one fixed pool without invoking the docking engine."""
    poses = read_poses(pool)
    rows = [{"model": pose.model, "native": pose.native_score, "native_rank": pose.native_rank}
            for pose in poses]
    return add_aromatic_order({"poses": rows, "orders": rank_orders(rows, ("native",)),
                               "abstentions": {}}, poses, receptor, ccd)


def mmff_properties(molecule):
    props = AllChem.MMFFGetMoleculeProperties(molecule, mmffVariant=MMFF_VARIANT)
    if props is None:
        raise ValueError(f"No {MMFF_VARIANT} parameters for prepared ligand")
    return props


def minimized_energy(molecule, properties, conformer_id, fixed_atoms):
    ff = AllChem.MMFFGetMoleculeForceField(molecule, properties, confId=conformer_id,
                                           nonBondedThresh=100.0)
    for index in fixed_atoms:
        ff.AddFixedPoint(index)
    ff.Initialize()
    if ff.Minimize(maxIts=500) and ff.Minimize(maxIts=1500):
        raise ValueError("MMFF minimization did not converge")
    energy = ff.CalcEnergy()
    if not isfinite(energy):
        raise ValueError("Nonfinite MMFF energy")
    return energy


def gas_reference(ligand):
    molecule = Chem.Mol(ligand)
    molecule.RemoveAllConformers()
    config = AllChem.ETKDGv3()
    config.randomSeed = GAS_EMBED_SEED
    config.numThreads = 1
    conformers = AllChem.EmbedMultipleConfs(molecule, numConfs=GAS_CONFORMERS, params=config)
    properties = mmff_properties(molecule)
    energies = []
    for conf_id in conformers:
        try:
            energies.append(minimized_energy(molecule, properties, conf_id, ()))
        except ValueError:
            continue
    if not energies:
        raise ValueError("No gas-phase MMFF reference conformer converged")
    return min(energies), len(energies)


def strain_for_pool(poses):
    started = perf_counter()
    ligand, identity = decode_ligand(poses[0])
    reference, n_reference = gas_reference(ligand)
    energies = {}
    for pose in poses:
        molecule, current = decode_ligand(pose)
        if current != identity:
            raise ValueError("Prepared ligand chemistry differs across the pose pool")
        conformer = molecule.GetConformer()
        atoms = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum() != 1]
        xyz = conformer.GetPositions()[atoms]
        energy = minimized_energy(molecule, mmff_properties(molecule), conformer.GetId(), atoms)
        if not np.array_equal(xyz, conformer.GetPositions()[atoms]):
            raise ValueError("MMFF constrained relaxation modified the frozen heavy coordinates")
        energies[pose.model] = energy
    return {"identity": identity, "gas_reference_kcal_mol": reference,
            "gas_conformers_converged": n_reference, "pose_kcal_mol": energies,
            "wall_seconds": perf_counter() - started}


def mmff_affinity_score(pose, scale, bound_kcal_mol, free_kcal_mol):
    """Replace the saved, scaled Vina intraligand contribution at unit weight."""
    if pose.native_intra is None or pose.native_unbound is None:
        raise ValueError("MMFF affinity requires saved Vina INTRA and UNBOUND energies")
    return pose.native_score + scale * (bound_kcal_mol - free_kcal_mol
                                        - (pose.native_intra - pose.native_unbound))


def mmff_affinity_annotation(poses, strain, scale):
    """Two affinity summaries for one seeded pool, without changing pose order."""
    if len({p.native_unbound for p in poses}) != 1:
        raise ValueError("MMFF affinity requires a common Vina unbound reference")
    scores = [{"model": p.model, "native_rank": p.native_rank,
               "affinity_score": mmff_affinity_score(
                   p, scale, strain["pose_kcal_mol"][p.model], strain["gas_reference_kcal_mol"])}
              for p in poses]
    first = min(scores, key=lambda r: r["native_rank"])
    best = min(scores, key=lambda r: (r["affinity_score"], r["native_rank"]))
    return {"method": "MMFF94s unit-weight gas-phase intraligand replacement",
            "reference": {"energy_kcal_mol": strain["gas_reference_kcal_mol"],
                          "conformers_attempted": GAS_CONFORMERS,
                          "conformers_converged": strain["gas_conformers_converged"],
                          "embedding_seed": GAS_EMBED_SEED},
            "vina_flexibility_scale": scale, "native_rank_one": first,
            "minimum_across_emitted_models": best, "per_model": scores}


def receptor_atoms(path):
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM ", "HETATM ")):
            ident = line[6:16], line[17:27], line[70:].strip()
            ad_type = line.split()[-1]
            position = [float(line[i:i + 8]) for i in (30, 38, 46)]
            atoms.append((ident, ad_type, position))
    if not atoms:
        raise ValueError(f"No receptor atom records in {path}")
    return atoms


def protein_clash_delta(original, alternate):
    """Use only pairs changed by a side-chain shift, excluding intraresidue pairs."""
    if [a[:2] for a in original] != [a[:2] for a in alternate]:
        raise ValueError("Alternate receptor changes PDBQT atom identity, order, charge or type")
    old = np.array([a[2] for a in original])
    new = np.array([a[2] for a in alternate])
    changed = np.any(new != old, axis=1)
    if not changed.any():
        raise ValueError("Alternate receptor has no moved atoms")
    backbone = {"N", "CA", "C", "O", "OXT"}
    if any(atom[0][0][6:].strip() in backbone and moved for atom, moved in zip(original, changed)):
        raise ValueError("Alternate receptor changes a backbone coordinate")
    # A flat heavy-atom list supports matched pairs and avoids penalizing
    # unchanged contacts or adding the same two moved atoms twice.
    heavy = np.array([typ not in RECEPTOR_HYDROGEN for _, typ, _ in original])
    moved = np.flatnonzero(changed & heavy)
    others = np.flatnonzero(heavy)
    if not moved.size:
        raise ValueError("Alternate receptor only moves hydrogens")
    residue = np.asarray([atom[0][1] for atom in original])
    valid = ((residue[moved, None] != residue[None, others])
             & (~changed[None, others] | (moved[:, None] < others[None, :])))
    radii = np.asarray([ATOM_RADII.get(atom[1], 1.7) for atom in original])
    contact = radii[moved, None] + radii[None, others] - CLASH_TOLERANCE_ANGSTROM

    def penalty(xyz):
        distance = np.linalg.norm(xyz[moved, None, :] - xyz[None, others, :], axis=-1)
        return CLASH_WEIGHT * float(np.square(np.maximum(contact - distance, 0)[valid]).sum())

    return penalty(new) - penalty(old)


def load_alternates(baseline, evidence):
    """Load the observed state probabilities and calculate clash costs once per receptor."""
    if evidence is None:
        return []
    spec = json.loads(evidence.read_text())
    p0 = float(spec["baseline_population"])
    if not 0 < p0 <= 1:
        raise ValueError("Baseline population must be in (0, 1]")
    original = receptor_atoms(baseline)
    states = []
    for state in spec["alternates"]:
        population = float(state["population"])
        if not 0 < population <= 1:
            raise ValueError("Alternate population must be in (0, 1]")
        receptor = (evidence.parent / state["receptor"]).resolve(strict=True)
        states.append({"id": state["id"], "receptor": receptor, "population": population,
                       "prior_kcal_mol": -RT_KCAL * log(population / p0),
                       "clash_kcal_mol": protein_clash_delta(original, receptor_atoms(receptor))})
    if len({s["id"] for s in states}) != len(states):
        raise ValueError("Need unique alternate state identifiers")
    if abs(p0 + sum(s["population"] for s in states) - 1.0) > 1e-6:
        raise ValueError("Baseline and alternate observed state populations must sum to 1")
    return states


def rank_orders(rows, arms):
    return {arm: [row["model"] for row in sorted(
        (row for row in rows if row[arm] is not None),
        key=lambda row: (-row[arm] if arm == "gnina" else row[arm], row["native_rank"], row["model"]))]
        for arm in arms}


def score_pool(binary, receptor, pool, *, waters=None, sidechains=None, timeout=120, satisfaction=False,
               gnina_scores_sdf=None, aromatic_ccd=None, prepared_ligand=None, ligand_bonds=None,
               mmff_affinity=False):
    """Return all model scores and state evidence without reading crystal outcomes."""
    if mmff_affinity and not prepared_ligand:
        raise ValueError("MMFF affinity requires the source --prepared-ligand topology")
    if prepared_ligand and aromatic_ccd:
        raise ValueError("Prepared topology rescoring cannot be combined with aromatic ordering")
    poses = read_poses(pool)
    coordinate_inputs = fixed_coordinate_inputs(poses, prepared_ligand, ligand_bonds)
    cnn_scores = read_gnina_scores(poses, gnina_scores_sdf) if gnina_scores_sdf else None
    reasons = {}
    try:
        [rank1] = [pose for pose in poses if pose.native_rank == 1]
        unbound = rank1.native_unbound
        if not prepared_ligand and (unbound is None or rank1.native_intra != unbound or any(
            pose.native_intra is None or pose.native_unbound != unbound for pose in poses
        )):
            raise ValueError("Missing common saved Vina unbound reference or nonzero native rank-1 intra minus unbound")
        strain = strain_for_pool(poses)
    except (RuntimeError, ValueError) as error:
        strain = None
        reasons["strain"] = str(error)
    if not waters:
        reasons["water"] = "No predeclared independent water evidence"
        if satisfaction:
            reasons["satisfied_water"] = reasons["water"]
    try:
        alternates = load_alternates(receptor, sidechains)
    except (OSError, ValueError, KeyError) as error:
        alternates = []
        reasons["sidechain"] = str(error)
    if not alternates:
        reasons.setdefault("sidechain", "No distinguishable alternate in the registered evidence"
                           if sidechains else "No predeclared independent sidechain evidence")

    rows = []
    with TemporaryDirectory(prefix="vina-physical-pool-") as work:
        ligand = Path(work) / "model.pdbqt"
        reference_terms = None
        if prepared_ligand:
            ligand.write_text(rank1.text)
            reference_terms = evaluate_coordinates(binary, ligand, receptor, waters, timeout,
                                                   satisfaction=satisfaction, **coordinate_inputs)
            unbound = reference_terms["intra"]
        for pose in poses:
            ligand.write_text(pose.text)
            original = (reference_terms if reference_terms is not None and pose == rank1 else
                        evaluate_coordinates(binary, ligand, receptor, waters, timeout,
                                             satisfaction=satisfaction, **coordinate_inputs))
            scale = original["scale"]
            if prepared_ligand and reference_terms["scale"] != scale:
                raise ValueError("Fixed coordinates changed the prepared-ligand normalization")
            native = (original["zero_reference"] - scale * unbound) if prepared_ligand else pose.native_score
            water_credit = original["water_delta"] - original["combined_delta"]
            joint_polar = original["water_delta"] if waters else original["combined_delta"]
            if strain is None:
                strain_delta = 0.0
                pose_strain = strain_score = native_relative = None
            else:
                pose_strain = strain["pose_kcal_mol"][pose.model] - strain["gas_reference_kcal_mol"]
                intra = original["intra"] if prepared_ligand else pose.native_intra
                native_relative = intra - unbound
                strain_delta = scale * (STRAIN_WEIGHT * pose_strain - native_relative)
                strain_score = native + strain_delta
            state_rows = []
            for state in alternates:
                alternative = evaluate_coordinates(binary, ligand, state["receptor"], waters, timeout, **coordinate_inputs)
                if abs(scale - alternative["scale"]) > 1e-10:
                    raise ValueError("An alternate receptor unexpectedly changed ligand normalization")
                if abs(original["intra"] - alternative["intra"]) > 1e-8:
                    raise ValueError("An alternate rigid receptor changed fixed ligand internal energy")
                receptor_delta = scale * (alternative["inter"] - original["inter"]
                                          + state["prior_kcal_mol"] + state["clash_kcal_mol"])
                alt_polar = alternative["water_delta"] if waters else alternative["combined_delta"]
                state_rows.append({**state, "receptor": str(state["receptor"]),
                                   "score": native + receptor_delta,
                                   "joint_score": native + receptor_delta + strain_delta + alt_polar,
                                   "engine": alternative})
            baseline_joint = native + strain_delta + joint_polar
            best_sidechain, sidechain_state = min(
                [(native, ""), *((alt["score"], alt["id"]) for alt in state_rows)])
            best_joint, joint_state = min(
                [(baseline_joint, ""), *((alt["joint_score"], alt["id"]) for alt in state_rows)])
            row = {"model": pose.model, "native_rank": pose.native_rank, "native": native,
                   "direction": native + original["direction_delta"],
                   "burial": native + original["burial_delta"],
                   "combined": native + original["combined_delta"],
                   "water": native + original["water_delta"] if waters else None,
                   "water_credit_delta": water_credit if waters else None,
                   "strain": strain_score, "sidechain": best_sidechain if alternates else None,
                   "joint": best_joint, "engine": original, "alternate_states": state_rows,
                   "sidechain_state": (sidechain_state or "baseline") if alternates else None,
                   "joint_state": joint_state or "baseline",
                   "strain_relative_kcal_mol": pose_strain,
                   "native_relative_intra_kcal_mol": native_relative,
                   "saved_native_intra_kcal_mol": pose.native_intra,
                   "saved_native_unbound_kcal_mol": pose.native_unbound,
                   "strain_replacement_delta": strain_delta if strain else None}
            if prepared_ligand:
                row["saved_native"] = pose.native_score
            if satisfaction:
                row.update(
                    satisfied=native + original["satisfied_delta"],
                    satisfied_water=native + original["satisfied_water_delta"] if waters else None,
                )
            if cnn_scores is not None:
                row.update(gnina_rank_terms(pose, cnn_scores[pose.model], native_score=native))
            rows.append(row)
    arms = ARMS + SATISFACTION_ARMS if satisfaction else ARMS
    if cnn_scores is not None:
        arms += ("gnina", "cnn_consensus")
    result = {"poses": rows, "strain": strain, "abstentions": reasons,
            "joint_components": [key for key, available in (
                ("direction", True), ("burial", True), ("water", waters),
                ("strain", strain), ("sidechain", alternates)
            ) if available],
            "orders": rank_orders(rows, arms)}
    if prepared_ligand:
        result["fixed_topology"] = {
            "prepared_ligand": str(prepared_ligand), "ligand_bonds": str(ligand_bonds) if ligand_bonds else None,
            "reference_model": rank1.model, "unbound": unbound,
            "saved_native_order": [pose.model for pose in sorted(poses, key=lambda pose: pose.native_rank)],
        }
    if mmff_affinity:
        if strain is None:
            reasons["mmff_affinity"] = reasons["strain"]
            result["mmff_affinity"] = None
        else:
            try:
                result["mmff_affinity"] = mmff_affinity_annotation(poses, strain, reference_terms["scale"])
            except ValueError as error:
                reasons["mmff_affinity"] = str(error)
                result["mmff_affinity"] = None
    return add_aromatic_order(result, poses, receptor, aromatic_ccd) if aromatic_ccd else result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--receptor", type=Path)
    parser.add_argument("--pool", type=Path)
    parser.add_argument("--prepared-ligand", type=Path, help="Lock each physical score to this original ligand topology")
    parser.add_argument("--ligand-bonds", type=Path, help="Use validated full Meeko adjacency in the supporting build")
    parser.add_argument("--export-ligand-bonds", action="store_true",
                        help="Export the checked adjacency for one --prepared-ligand and exit")
    parser.add_argument("--mmff-affinity", action="store_true",
                        help="Report experimental unit-weight MMFF affinity for native rank one and the emitted pool minimum")
    parser.add_argument("--waters", type=Path)
    parser.add_argument("--sidechains", type=Path)
    parser.add_argument("--satisfaction", action="store_true",
                        help="Add direct and water-competed polar satisfaction arms (requires a supporting build)")
    parser.add_argument("--gnina-scores-sdf", type=Path,
                        help="Add CNNscore and native-energy consensus from a complete same-coordinate GNINA SDF")
    parser.add_argument("--aromatic-ccd", type=Path,
                        help="Add geometry-based aromatic ordering using a component CCD CIF")
    parser.add_argument("--vinardo-binary", type=Path, help="AutoDock Vina executable for fixed-coordinate Vinardo")
    parser.add_argument("--vinardo-box", type=Path, help="Existing registered Vina search box for the pool")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--export-gnina-sdf", action="store_true",
                        help="Write the fixed MODELs to one multi-molecule SDF for GNINA")
    args = parser.parse_args()
    if args.mmff_affinity and (not args.binary or not args.receptor or not args.prepared_ligand or
                               args.export_gnina_sdf or args.export_ligand_bonds):
        parser.error("--mmff-affinity requires --binary, --receptor, --prepared-ligand and a pose pool")
    if args.output.exists():
        parser.error("--output already exists")
    if args.export_ligand_bonds:
        if not args.prepared_ligand or any((args.pool, args.binary, args.receptor, args.ligand_bonds,
                args.waters, args.sidechains, args.satisfaction, args.gnina_scores_sdf, args.aromatic_ccd,
                args.vinardo_binary, args.vinardo_box, args.export_gnina_sdf)):
            parser.error("Bond export requires only --prepared-ligand and --output")
        try:
            edges = mapped_ligand_bonds(args.prepared_ligand.read_text())
        except (OSError, ValueError) as error:
            parser.error(str(error))
        with args.output.open("x") as stream:
            stream.writelines(f"{left} {right}\n" for left, right in edges)
        return
    if not args.pool:
        parser.error("Scoring and GNINA export require --pool")
    if args.ligand_bonds and not args.prepared_ligand:
        parser.error("--ligand-bonds requires --prepared-ligand")
    if args.prepared_ligand and (not (args.binary or args.vinardo_binary) or args.export_gnina_sdf or
                                 (args.binary and args.aromatic_ccd)):
        parser.error("--prepared-ligand requires a scoring binary and cannot change the baseline for aromatic ordering")
    if args.export_gnina_sdf:
        if args.gnina_scores_sdf or args.aromatic_ccd or args.vinardo_binary or args.vinardo_box:
            parser.error("--export-gnina-sdf cannot be combined with GNINA import, aromatic scoring or Vinardo")
        export_gnina_sdf(args.pool, args.output)
        return
    if bool(args.vinardo_binary) != bool(args.vinardo_box) or (args.vinardo_binary and not args.receptor):
        parser.error("Vinardo requires --vinardo-binary, --vinardo-box and --receptor")
    needs_physics = (not (args.gnina_scores_sdf or args.aromatic_ccd or args.vinardo_binary)
                     or args.satisfaction or args.waters or args.sidechains)
    if args.aromatic_ccd and not args.receptor:
        parser.error("Aromatic scoring requires --receptor")
    if (args.binary and not args.receptor) or (needs_physics and not args.binary) or (
        args.receptor and not (args.binary or args.aromatic_ccd or args.vinardo_binary)
    ):
        parser.error("Physical scoring requires --binary and --receptor")
    if args.binary:
        result = score_pool(args.binary.resolve(strict=True), args.receptor, args.pool,
                            waters=args.waters, sidechains=args.sidechains, satisfaction=args.satisfaction,
                            gnina_scores_sdf=args.gnina_scores_sdf, aromatic_ccd=args.aromatic_ccd,
                            prepared_ligand=args.prepared_ligand, ligand_bonds=args.ligand_bonds,
                            mmff_affinity=args.mmff_affinity)
    elif args.gnina_scores_sdf:
        result = score_gnina_pool(args.pool, args.gnina_scores_sdf)
        if args.aromatic_ccd:
            add_aromatic_order(result, read_poses(args.pool), args.receptor, args.aromatic_ccd)
    elif args.aromatic_ccd:
        result = score_aromatic_pool(args.pool, args.receptor, args.aromatic_ccd)
    else:
        rows = [{"model": pose.model, "native": pose.native_score, "native_rank": pose.native_rank}
                for pose in read_poses(args.pool)]
        result = {"poses": rows, "orders": rank_orders(rows, ("native",)), "abstentions": {}}
    if args.vinardo_binary:
        add_vinardo_order(result, read_poses(args.pool), args.vinardo_binary, args.receptor, args.vinardo_box,
                          prepared_ligand=args.prepared_ligand, ligand_bonds=args.ligand_bonds)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
