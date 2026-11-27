"""Prepare the prespecified cognate panel from PDB and CCD depositor files."""
import argparse
from collections import defaultdict
from importlib.metadata import version
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "pair_gradient/pose"))
import prepare as historical

SYSTEMS = {
    "4BFP": ("TNKS2", "SWY", "A"),
    "4I1R": ("MALT1", "LZU", "A"),
    "2AA2": ("NR3C2", "AS4", "A"),
    "2YIY": ("MAP3K7", "YIY", "A"),
}


PPAR = {
    "2PRG": ("PPARG", "BRL", "A"),
    "3TKM": ("PPARD", "GW0", "A"),
}


def deposition_audit(pdb, component, chain, meta):
    lines = pdb.read_text().splitlines()
    ligand = [s for s in lines if s.startswith("HETATM") and s[17:20].strip() == component and s[21] == chain and s[76:78].strip() not in {"H", "D"}]
    ligand_instances = sorted({s[21:27] for s in ligand})
    if len(ligand_instances) != 1:
        raise ValueError(f"Ambiguous ligand instances {ligand_instances}")
    reference = np.array([historical.base.coords(s) for s in ligand])
    groups = defaultdict(list)
    for s in lines:
        if not s.startswith(("ATOM  ", "HETATM")) or s[76:78].strip() in {"H", "D"}:
            continue
        if s in ligand:
            continue
        group = s[17:20].strip() if s.startswith("HETATM") else "protein"
        groups[s[21], group].append(s)
    environment = []
    for (sid, name), atoms in sorted(groups.items()):
        xyz = np.array([historical.base.coords(s) for s in atoms])
        distance = np.linalg.norm(xyz[:, None, :] - reference[None, :, :], axis=2)
        environment.append({"chain": sid, "component": name,
                            "min_distance_A": round(float(distance.min()), 3),
                            "heavy_atoms_within_4A": int((distance.min(axis=1) <= 4).sum()),
                            "retained_in_prepared_receptor": sid == chain and (name == "protein" or name in meta.get("metals", [])) if meta.get("receptor_status") == "ready" else None})
    return {"experimental_resolution": next(s.strip() for s in lines if s.startswith("REMARK   2 RESOLUTION.")),
            "ligand_instances": ligand_instances,
            "chemical_name": [s[15:].strip() for s in lines if s.startswith("HETNAM") and s[11:14].strip() == component],
            "ligand_LINK_records": [s for s in lines if s.startswith("LINK  ") and (s[17:20].strip() == component or s[47:50].strip() == component)],
            "environment": environment, "reference_atom_graph_maps": len(meta.get("mapped", [])),
            "mapped_heavy_atoms": len(meta.get("serial", [])), "reference_alternates": meta.get("n_altloc"),
            "deposition_DOI": f"https://doi.org/10.2210/pdb{pdb.stem}/pdb",
            "pdb_url": f"https://files.rcsb.org/download/{pdb.name}",
            "ccd_url": f"https://files.rcsb.org/ligands/download/{component}.cif",
            "publication_DOI": next((s[19:].strip() for s in lines if s.startswith("JRNL        DOI")), None)}


def main(work, systems):
    historical.CASES.update(systems)
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    result = {}
    for pdb, (target, component, chain) in systems.items():
        start = time.monotonic()
        try:
            historical.prepare(work, pdb)
            meta = json.loads((work / pdb / "preparation.json").read_text())
        except (ValueError, RuntimeError) as exc:
            meta = dict(target=target, ligand=component, receptor_status=str(exc))
            (work / pdb / "preparation.json").write_text(json.dumps(meta) + "\n")
        result[pdb] = {k: meta.get(k) for k in ("target", "ligand", "pdb_sha256", "ccd_sha256", "receptor_status", "closest_omitted_heavy_atom",
                                                 "closest_repaired_ligand_distance", "metals", "canonical_smiles")}
        result[pdb].update(deposition_audit(work / "raw" / (pdb + ".pdb"), component, chain, meta))
        result[pdb]["prep_wall_s"] = time.monotonic() - start
    (work / "sources.json").write_text(json.dumps(dict(systems=result, python=sys.version.split()[0],
        packages={p: version(p) for p in ("rdkit", "meeko", "pdbfixer", "OpenMM", "gemmi", "numpy", "scipy")}), indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("--ppar", action="store_true", help="Generate the separately registered second-stage PPAR panel")
    a = p.parse_args()
    main(a.work.resolve(), PPAR if a.ppar else SYSTEMS)
