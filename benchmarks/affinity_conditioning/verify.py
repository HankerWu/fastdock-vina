"""Independent FKBP AD4 map, symmetry and fixed-geometry checks on published coordinates."""
import argparse
import csv
import itertools
import json
import tarfile
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
GRID = BENCH / "cubic_grid/fkbp_ad4"
COVER = BENCH / "rigid_pose_coverage/fkbp/raw"
REC = ("1FKG", "1FKH")
ARMS = ("official", "official_high", "base", "high", "cost_pro", "edge",
        "feasible_edge", "coverage", "uniform", "uniform_work")


def members(archive, wanted):
    contents = {}
    with tarfile.open(archive, "r|xz") as tar:
        for obj in tar:
            if obj.name in wanted:
                contents[obj.name] = tar.extractfile(obj).read().decode("ascii")
    if contents.keys() != set(wanted):
        raise ValueError(f"Missing members in {archive}")
    return contents


class Ligand:
    def __init__(self, pdbqt):
        lines = pdbqt.splitlines()
        smiles, = [s.split()[2] for s in lines if s.startswith("REMARK SMILES ") and "IDX" not in s]
        self.graph = Chem.MolFromSmiles(smiles)
        pairs = [t for s in lines if s.startswith("REMARK SMILES IDX ")
                 for t in map(int, s.split()[3:])]
        serials = dict(zip(pairs[::2], pairs[1::2], strict=True))
        self.order = [serials[i + 1] for i in range(self.graph.GetNumAtoms())]
        by_id = {int(s[6:11]): s for s in lines if s.startswith(("ATOM  ", "HETATM"))}
        atoms = [by_id[i] for i in self.order]
        self.kinds = np.array([s[77:].strip() for s in atoms])
        self.q = np.array([float(s[70:76]) for s in atoms])
        self.tors, = [int(s.split()[1]) for s in lines if s.startswith("TORSDOF ")]
        self.reference = self.coords(lines)
        checked = Chem.Mol(self.graph)
        Chem.RemoveStereochemistry(checked)
        conformer = Chem.Conformer(checked.GetNumAtoms())
        conformer.Set3D(True)
        for i, position in enumerate(self.reference):
            conformer.SetAtomPosition(i, position)
        checked.AddConformer(conformer)
        Chem.AssignStereochemistryFrom3D(checked, replaceExistingTags=True)
        if Chem.MolToSmiles(checked) != Chem.MolToSmiles(self.graph):
            raise ValueError("Prepared 3D configuration differs from the recorded isomer")
        self.bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.graph.GetBonds()])

    def coords(self, lines):
        sites = {int(s[6:11]): [float(s[i:i+8]) for i in (30, 38, 46)]
                 for s in lines if s.startswith(("ATOM  ", "HETATM"))}
        if set(sites) != set(self.order):
            raise ValueError("Inconsistent particle identities")
        return np.array([sites[i] for i in self.order])

    def check_center(self, xyz):
        # A tetrahedral center with three heavy neighbors and an implicit H.
        [(idx, cip)] = Chem.FindMolChiralCenters(self.graph, includeUnassigned=True)
        neighbors = [a.GetIdx() for a in self.graph.GetAtomWithIdx(idx).GetNeighbors()]
        if len(neighbors) != 3 or cip != "S":
            raise ValueError("Unexpected pipecolate stereochemistry")
        def volume(x):
            return np.linalg.det(x[neighbors] - x[idx])
        if volume(self.reference) * volume(xyz) <= 0:
            raise ValueError("Reflected stereocenter")


class Fields:
    def __init__(self, receptor):
        kinds = ("A", "C", "N", "OA", "e", "d")
        keys = [f"maps/{receptor}/receptor.{t}.map" for t in kinds]
        files = members(GRID / f"maps-{receptor}.tar.xz", keys)
        self.values = {}
        for name, k in zip(keys, kinds, strict=True):
            body = files[name].split("\n", 6)
            spacing = float(body[3].split()[1])
            spans = np.array(list(map(int, body[4].split()[1:])))
            center = np.array(list(map(float, body[5].split()[1:])))
            lattice = (spacing, tuple(spans), tuple(center))
            if self.values and lattice != self.lattice:
                raise ValueError("Different AD4 lattices")
            self.lattice = lattice
            self.low = center - spacing * spans / 2
            self.values[k] = np.fromstring(body[6], sep=" ").reshape(tuple(spans+1), order="F")

    def score(self, ligand, xyz, official):
        spacing, spans, _ = self.lattice
        fraction = (xyz-self.low) / spacing
        u = np.clip(fraction, 0, spans)
        if np.max(np.abs(fraction-u))*spacing > 1e-8:
            raise ValueError("Atom outside map; this comparison only treats in-box poses")
        lower = np.minimum(np.floor(u).astype(int), np.array(spans)-1)
        offset = u-lower
        by_type = {}
        for k, tensor in self.values.items():
            sampled = np.zeros(len(xyz))
            for corner in itertools.product((0,1), repeat=3):
                weight = np.prod(np.where(corner, offset, 1-offset), axis=1)
                sampled += weight * tensor[tuple((lower+corner).T)]
            by_type[k] = sampled
        def cap(x):
            return x*1000/(1000+np.maximum(x, 0))
        typed = np.array([by_type[t][i] for i, t in enumerate(ligand.kinds)])
        e, d, q = by_type["e"], by_type["d"], ligand.q
        if official:
            atom = cap(typed) + q*cap(e) + abs(q)*cap(d)
        else:
            atom = cap(typed + q*e + abs(q)*d)
        return float(np.sum(atom) + 0.2983*ligand.tors)


def active_pools():
    for arm in ARMS:
        files = {f"pools/{arm}/{rec}_{n}.pdbqt": (rec,n) for rec in REC for n in (5,6)}
        pools = members(COVER / f"{arm}.tar.xz", files)
        yield "coverage", arm, 1047349, {files[k]: v for k,v in pools.items()}
    keys = {f"pools/{n}_{r}_{s}_ad4.pdbqt": (s,r,n) for s in (131071,196613,262147)
            for r in REC for n in (5,6)}
    pools = members(GRID / "pools.tar.xz", keys)
    for seed in (131071,196613,262147):
        yield "official_legacy", "official", seed, {(r,n):pools[k] for k,(s,r,n) in keys.items() if s==seed}


def modes(pdbqt):
    pieces = pdbqt.split("ENDMDL")
    if pieces[-1].strip() or len(pieces) != 10:
        raise ValueError("Expected nine complete models")
    return [part.splitlines() for part in pieces[:-1]]


def main(compare, output):
    files = members(BENCH / "cubic_grid/human_panel/affinity-pools.tar.xz",
                    ["affinity/5.pdbqt", "affinity/6.pdbqt"])
    lig5, lig6 = (Ligand(files[f"affinity/{i}.pdbqt"]) for i in (5,6))
    # Full, non-unique matches expose both methyl exchange and phenyl reflection.
    mappings = np.array(lig6.graph.GetSubstructMatches(lig5.graph, useChirality=True, uniquify=False))
    deleted = {tuple(sorted(set(range(33))-set(ids))) for ids in mappings}
    if len(mappings) != 4 or len(deleted) != 1:
        raise ValueError("Unexpected scaffold degeneracy")
    remove = deleted.pop()
    delete_serials = sorted(lig6.order[i] for i in remove)
    stripped = Chem.RWMol(lig6.graph)
    for i in reversed(remove):
        stripped.RemoveAtom(int(i))
    Chem.SanitizeMol(stripped)
    if Chem.MolToSmiles(stripped) != Chem.MolToSmiles(lig5.graph):
        raise ValueError("After hydrogen completion, the deletion is not 5")
    fields = {rec:Fields(rec) for rec in REC}
    expected = {}
    with compare.open() as stream:
        for row in csv.DictReader(stream):
            key = (row["panel"],row["arm"],int(row["seed"]),row["receiver"],int(row["rank"]))
            if key in expected:
                raise ValueError(f"Duplicate row {key}")
            expected[key] = row
    error = symmetry = printed = 0.
    length_error, checked = [], 0
    b0,b1=lig5.bonds.T
    ideal = np.linalg.norm(lig5.reference[b0]-lig5.reference[b1],axis=1)
    for panel, arm, seed, pools in active_pools():
        official = arm.startswith("official")
        for rec in REC:
            for rank, (a,b) in enumerate(zip(modes(pools[rec,5]),modes(pools[rec,6]),strict=True), 1):
                p5,p6=lig5.coords(a),lig6.coords(b)
                lig5.check_center(p5)
                lig6.check_center(p6)
                S = fields[rec].score
                s5,s6 = S(lig5,p5,official),S(lig6,p6,official)
                if rank == 1:
                    first5=s5
                variants=[]
                worst_bond_delta=[]
                for match in mappings:
                    xyz = p6[match]
                    lig5.check_center(xyz)
                    variants.append(S(lig5,xyz,official))
                    worst_bond_delta.append(np.max(np.abs(np.linalg.norm(xyz[b0]-xyz[b1],axis=1)-ideal)))
                symmetry = max(symmetry, np.ptp(variants))
                length_error.append(min(worst_bond_delta))
                printed_values=[float(next(s for s in lines if s.startswith("REMARK VINA RESULT")).split()[3])
                                for lines in (a,b)]
                printed=max(printed, abs(s5-printed_values[0]), abs(s6-printed_values[1]))
                row=expected.pop((panel,arm,seed,rec,rank))
                for name, value in {"s5_native":s5, "s6_native":s6, "s5_projected":variants[0],
                                    "substitution":s6-variants[0],"placement":variants[0]-first5,
                                    "net":s6-first5}.items():
                    error = max(error, abs(float(row[name])-value))
                checked+=1
    if expected or error>5e-8 or symmetry>1e-10 or printed>0.01:
        raise ValueError(f"QA failed rows={len(expected)}, manual={error}, symmetry={symmetry}, print={printed}")
    common = mappings[0]
    charge_shift = lig6.q[common] - lig5.q
    # A published atom position is rounded by at most 0.0005 Å per component.
    # A distance between two endpoints contributes at most sqrt(3)*0.001 Å,
    # and two independently rounded distances double this bound.
    rounding_bound = 2 * np.sqrt(3.0) * 0.001
    if min(length_error) <= rounding_bound:
        raise ValueError("Cannot exclude exact fixed-bond reachability at file precision")
    report = {
        "comparison_rows": checked,
        "native_modes": 2 * checked,
        "scaffold_embeddings": len(mappings),
        "projected_stereo_checks": len(mappings) * checked,
        "removed_pdbqt_serials": delete_serials,
        "formal_charges": [Chem.GetFormalCharge(l.graph) for l in (lig5, lig6)],
        "formulas": [rdMolDescriptors.CalcMolFormula(l.graph) for l in (lig5, lig6)],
        "printed_partial_charge_sums": [float(l.q.sum()) for l in (lig5, lig6)],
        "retained_charge_without_retyping": float(lig6.q[common].sum()),
        "retained_atoms_with_changed_partial_charge": int(np.count_nonzero(charge_shift)),
        "largest_retained_atom_partial_charge_change": float(max(abs(charge_shift))),
        "retained_atoms_with_changed_ad4_type": int(np.sum(lig5.kinds != lig6.kinds[common])),
        "max_error_six_score_columns": error,
        "max_symmetry_score_spread": symmetry,
        "max_native_score_print_discrepancy": printed,
        "range_of_best_mapping_max_bond_length_change_A": [min(length_error), max(length_error)],
        "coordinate_rounding_bound_for_two_bond_lengths_A": rounding_bound,
    }
    output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("poses", type=Path)
    parser.add_argument("report", type=Path)
    args=parser.parse_args()
    main(args.poses,args.report)
