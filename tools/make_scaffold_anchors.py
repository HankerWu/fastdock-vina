"""Write all graph-supported point assignments from one independently aligned donor SDF."""

import argparse
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS
from score_physical_pool import mapped_ligand_bonds, serial_index_pairs


def mcs_parameters():
    params = rdFMCS.MCSParameters()
    params.MaximizeBonds = False
    params.Timeout = 60
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrderExact
    params.AtomCompareParameters.MatchFormalCharge = True
    params.AtomCompareParameters.MatchChiralTag = True
    params.AtomCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.RingMatchesRingOnly = True
    params.BondCompareParameters.CompleteRingsOnly = True
    return params


def stereo_constraints(molecule):
    """Use RDKit's atom/bond neighbor order, including each alkene's defining atoms."""
    tetra = {}
    for atom in molecule.GetAtoms():
        tag = atom.GetChiralTag()
        if tag in (
            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        ):
            neighbors = tuple(n.GetIdx() for n in atom.GetNeighbors())
            tetra[atom.GetIdx()] = (tag, neighbors)

    double = {}
    for bond in molecule.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in (
            Chem.BondStereo.STEREOZ,
            Chem.BondStereo.STEREOE,
            Chem.BondStereo.STEREOCIS,
            Chem.BondStereo.STEREOTRANS,
        ):
            continue
        start, end = bond.GetBeginAtom(), bond.GetEndAtom()
        left, right = bond.GetStereoAtoms()
        sides = [
            (
                ref,
                frozenset(
                    n.GetIdx()
                    for n in atom.GetNeighbors()
                    if n.GetIdx() != other.GetIdx()
                ),
            )
            for ref, atom, other in ((left, start, end), (right, end, start))
        ]
        opposite = stereo in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
        double[start.GetIdx(), end.GetIdx()] = (opposite, *sides)
    return tetra, double


def stereo_compatible(mapping, query_constraints, donor_constraints):
    """Reject an observable inversion; partial or unspecified stereo supplies no veto.

    A tetrahedral center needs both specified labels and the complete correspondence
    between explicit neighbors. With three explicit neighbors, the implicit fourth
    position has the same place in both RDKit neighbor orders. An alkene needs two
    specified geometries and at least one corresponding substituent on each end.
    """
    qtetra, qdouble = query_constraints
    dtetra, ddouble = donor_constraints
    for q, (qtag, qneighbors) in qtetra.items():
        donor_center = dtetra.get(mapping.get(q))
        if donor_center is None:
            continue
        dtag, dneighbors = donor_center
        mapped = [mapping.get(i) for i in qneighbors]
        if len(mapped) != len(dneighbors) or set(mapped) != set(dneighbors):
            continue
        order = [dneighbors.index(i) for i in mapped]
        odd = sum(a > b for j, a in enumerate(order) for b in order[j + 1 :]) % 2
        if bool(odd) != (qtag != dtag):
            return False

    for (start, end), (qopposite, qleft, qright) in qdouble.items():
        ends = (mapping.get(start), mapping.get(end))
        if ends in ddouble:
            dopposite, dleft, dright = ddouble[ends]
        elif ends[::-1] in ddouble:
            dopposite, dright, dleft = ddouble[ends[::-1]]
        else:
            continue
        switched = False
        for (qref, qside), (dref, dside) in ((qleft, dleft), (qright, dright)):
            substituent = next(
                (q for q in sorted(qside) if mapping.get(q) in dside), None
            )
            if substituent is None:
                break
            switched ^= (substituent != qref) != (mapping[substituent] != dref)
        else:
            if switched != (qopposite != dopposite):
                return False
    return True


class StereoMatch(rdFMCS.MCSFinalMatchCheck):
    """Check a growing FMCS embedding against the supplied molecule-specific stereo."""

    def __init__(self, constraints):
        super().__init__()
        self.constraints = constraints

    def __call__(self, left, right, atom_pairs, _bond_pairs, _parameters):
        return stereo_compatible(
            dict(atom_pairs),
            self.constraints[left.GetIntProp("_stereo_source")],
            self.constraints[right.GetIntProp("_stereo_source")],
        )


def all_correspondences(query, donor):
    params = mcs_parameters()
    params.StoreAll = True
    query_stereo = stereo_constraints(query)
    donor_stereo = stereo_constraints(donor)
    molecules = [query, donor]
    if query_stereo[1] or donor_stereo[1]:
        # RDKit's tetrahedral FMCS final check can reject even an alkene matched
        # to itself. Use explicit local orientation checks for these inputs.
        params.AtomCompareParameters.MatchChiralTag = False
        params.FinalMatchChecker = StereoMatch((query_stereo, donor_stereo))
        molecules = [Chem.Mol(mol) for mol in molecules]
        for index, mol in enumerate(molecules):
            mol.SetIntProp("_stereo_source", index)
    result = rdFMCS.FindMCS(molecules, params)
    if result.canceled:
        raise ValueError("Incomplete maximum common substructure search")
    if result.numAtoms < 6:
        raise ValueError("The two molecules have fewer than six common scaffold atoms")
    mappings = set()
    # FMCS uses the charge constraint while constructing the common graph.
    # Some resulting SMARTS encode only the element. Recheck every embedding pair.
    qcharges = [atom.GetFormalCharge() for atom in query.GetAtoms()]
    dcharges = [atom.GetFormalCharge() for atom in donor.GetAtoms()]
    embedding_pairs = 0
    for pattern in result.degenerateSmartsQueryMolDict.values():
        matches = [
            mol.GetSubstructMatches(
                pattern, useChirality=True, uniquify=False, maxMatches=100001
            )
            for mol in (query, donor)
        ]
        if any(not group or len(group) >= 100001 for group in matches):
            raise ValueError("Incomplete enumeration of substructure embeddings")
        embedding_pairs += len(matches[0]) * len(matches[1])
        if embedding_pairs > 1000000:
            raise ValueError("Too many paired embeddings for exact enumeration")
        for query_atoms in matches[0]:
            for donor_atoms in matches[1]:
                pairs = tuple(zip(query_atoms, donor_atoms, strict=True))
                if all(
                    qcharges[q] == dcharges[d] for q, d in pairs
                ) and stereo_compatible(dict(pairs), query_stereo, donor_stereo):
                    mappings.add(tuple(sorted(pairs)))
    if not mappings:
        raise ValueError("No valid assignment from selected donor")
    return sorted(mappings)


def query_graph(text):
    mapped_ligand_bonds(text)
    lines = text.splitlines()
    (smiles,) = [
        line.split()[2]
        for line in lines
        if line.startswith("REMARK SMILES ")
        and not line.startswith("REMARK SMILES IDX ")
    ]
    graph = Chem.MolFromSmiles(smiles)
    serials = {
        i - 1: serial for i, serial in serial_index_pairs(lines, "REMARK SMILES IDX ")
    }
    return graph, serials


def format_assignments(mappings, serials, coordinates):
    blocks = []
    for mapping in mappings:
        rows = []
        for q, d in mapping:
            xyz = coordinates[d]
            if len(xyz) != 3 or not all(math.isfinite(x) for x in xyz):
                raise ValueError("Incomplete or nonfinite donor coordinates")
            rows.append(f"{serials[q]} " + " ".join(f"{v:.8f}" for v in xyz))
        blocks.append("\n".join(rows))
    return "\nOR\n".join(blocks) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "query",
        type=Path,
        help="Prepared single-ligand PDBQT with Meeko SMILES and serial map",
    )
    parser.add_argument(
        "donor",
        type=Path,
        help="Single bound reference ligand SDF, aligned into the docking receptor frame",
    )
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    query, serials = query_graph(args.query.read_text())
    records = Chem.SDMolSupplier(str(args.donor), removeHs=True)
    if len(records) != 1 or records[0] is None:
        parser.error("donor must contain one valid structure")
    donor = records[0]
    conformer = donor.GetConformer()
    if not conformer.Is3D():
        parser.error("donor must contain a 3D experimental conformation")
    maps = all_correspondences(query, donor)
    args.output.write_text(format_assignments(maps, serials, conformer.GetPositions()))
    print(f"{len(maps)} assignments on {len(maps[0])} common heavy atoms")


if __name__ == "__main__":
    main()
