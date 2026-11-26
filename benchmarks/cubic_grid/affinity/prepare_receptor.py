"""Select the ligand-bound chain and run the previously validated repository preparer."""

import hashlib
import json
import sys
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial.distance import cdist

repo, pb, out = [Path(p).resolve() for p in sys.argv[1:]]
sys.path.insert(0, str(repo / "benchmarks/torsion_kinase_affinity"))
from prepare_receptors import prepare

(out / "raw").mkdir(parents=True, exist_ok=True)
audit = (
    json.loads((out / "structure_audit.json").read_text())
    if (out / "structure_audit.json").exists()
    else {}
)
for target, pdb, ccd, chain in [
    ("SYK", "6ZCY", "QF8", "A"),
    ("FAK", "6YT6", "PKE", "B"),
]:
    if os.environ.get("PREP_TARGET") and os.environ["PREP_TARGET"] != target:
        continue
    d = pb / f"{pdb}_{ccd}"
    p, s = next(d.glob("*_protein.pdb")), next(d.glob("*_ligand.sdf"))
    lines = p.read_text().splitlines()
    molecule = Chem.SDMolSupplier(str(s), removeHs=True)[0]
    points = molecule.GetConformer().GetPositions()
    # Temporary ligand records supply the original crystal coordinates to the
    # existing preparation code for box centering and repair-distance audits.
    records = [
        l[:21] + "A" + l[22:]
        for l in lines
        if l.startswith("ATOM  ") and l[21] == chain
    ]
    # 6ZCY Leu404 has only N, CA, C, O and CB observed. The standard
    # side-chain repair puts the CD2 carbon 0.98 A from CB. Preserve the
    # observed heavy atoms as Ala, beyond 10 A from the cocrystal ligand.
    trunc = {}
    if pdb == "6ZCY":
        leu = [l for l in records if l[22:26].strip() == "404"]
        assert {l[12:16].strip() for l in leu} == {"N", "CA", "C", "O", "CB"}
        xyz = [[float(l[k : k + 8]) for k in (30, 38, 46)] for l in leu]
        distance = float(cdist(xyz, points).min())
        assert distance > 10
        trunc["404"] = {"original": "LEU", "distance_A": distance}
        records = [l[:17] + "ALA" + l[20:] if l in leu else l for l in records]
    # Complete the observed C terminus analytically in its peptide plane.
    # PDBFixer can otherwise place OXT in contact with the carbonyl oxygen.
    end_id = records[-1][22:26]
    terminal = {l[12:16].strip(): l for l in records if l[22:26] == end_id}
    if "OXT" not in terminal:

        def coord(line):
            return np.array([float(line[k : k + 8]) for k in (30, 38, 46)])

        c = coord(terminal["C"])
        bonds = [coord(terminal[a]) - c for a in ("O", "CA")]
        direction = sum(b / np.linalg.norm(b) for b in bonds)
        pos = c - 1.25 * direction / np.linalg.norm(direction)
        records.append(
            terminal["O"][:12]
            + " OXT"
            + terminal["O"][16:30]
            + "".join(f"{v:8.3f}" for v in pos)
            + terminal["O"][54:]
        )
    for i, (v, atom) in enumerate(zip(points, molecule.GetAtoms()), 1):
        elem = atom.GetSymbol().upper()
        records.append(
            f"HETATM{i:5} {elem:>4} {ccd:>3} A 999    "
            + "".join(f"{a:8.3f}" for a in v)
            + f"  1.00  0.00          {elem:>2}"
        )
    (out / "raw" / f"{pdb}.pdb").write_text("\n".join(records + ["END", ""]))
    omit = defaultdict(list)
    for l in lines:
        if l.startswith(("HETATM", "ATOM  ")) and l[76:78].strip() not in ("H", "D"):
            if l.startswith("HETATM") or l[21] != chain:
                label = f"{l[:6].strip()}:{l[21]}:{l[17:20].strip() if l.startswith('HETATM') else 'protein'}"
                omit[label].append([float(l[k : k + 8]) for k in (30, 38, 46)])
    ring_repair = None
    try:
        prepare(out, pdb, ccd)
    except subprocess.CalledProcessError:
        if pdb != "6ZCY":
            raise
        # The source only resolves Phe382 through CB. Its reconstructed ring
        # contacts Lys402. Keep all observed coordinates and pick the least
        # clashing ring rotamer solely from the protein structure.
        from repair_phe import place_ring, place_distal

        protein = out / "prepared" / (pdb + "-protein.pdb")
        lines_after = protein.read_text().splitlines()
        ring_repair = place_ring(lines_after, records)
        if not os.environ.get("REPLAY_REJECTED_SYK_REPAIR"):
            ring_repair["Lys386"] = place_distal(
                lines_after, 386, "CG", "CD", ("CE", "NZ")
            )
            ring_repair["Ser379"] = place_distal(lines_after, 379, "CA", "CB", ("OG",))
        protein.write_text("\n".join(lines_after) + "\n")
        with (out / "prepared" / (pdb + ".prep.log")).open("a") as f:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "meeko.cli.mk_prepare_receptor",
                    "--read_pdb",
                    protein,
                    "-p",
                    out / "prepared" / (pdb + ".pdbqt"),
                ],
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True,
            )
    audit[target] = dict(
        pdb=pdb,
        chain=chain,
        component=ccd,
        additional_truncations=trunc,
        ring_repair=ring_repair,
        source_pdb_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        source_sdf_sha256=hashlib.sha256(s.read_bytes()).hexdigest(),
        nearest_omitted_A={v: float(cdist(a, points).min()) for v, a in omit.items()},
        prep=json.loads((out / "prepared" / f"{pdb}.json").read_text()),
    )
    (out / "structure_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
