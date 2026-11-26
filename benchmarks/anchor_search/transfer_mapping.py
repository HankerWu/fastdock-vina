"""Select one independent scaffold donor using only the prepared ligand graphs."""

from rdkit import Chem
from rdkit.Chem import rdFMCS


def canonical_match(mol, pattern):
    matches = mol.GetSubstructMatches(
        pattern, uniquify=False, useChirality=True, maxMatches=100001
    )
    if not matches or len(matches) > 100000:
        raise ValueError("Missing or excessive scaffold mappings")
    return min(matches)


def choose_reference(query_id, molecules):
    """Return (donor ID, matched query atoms, matched donor atoms) or None."""
    query = molecules[query_id]
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
    best = None
    for donor_id, donor in sorted(molecules.items()):
        if donor_id == query_id:
            continue
        result = rdFMCS.FindMCS([query, donor], params)
        if result.canceled:
            raise ValueError(f"MCS timed out for {query_id}/{donor_id}")
        if result.numAtoms < 6:
            continue
        overlap = 2 * result.numAtoms / (query.GetNumAtoms() + donor.GetNumAtoms())
        order = (-result.numAtoms, -overlap, donor_id)
        if best is None or order < best[0]:
            pattern = Chem.MolFromSmarts(result.smartsString)
            best = (
                order,
                (
                    donor_id,
                    canonical_match(query, pattern),
                    canonical_match(donor, pattern),
                ),
            )
    return best[1] if best is not None else None
