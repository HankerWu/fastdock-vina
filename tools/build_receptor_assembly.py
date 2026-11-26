"""Build a rigid receptor from deposited PDB REMARK 350 assembly transformations."""

import argparse
import json
from pathlib import Path
import string

import gemmi
import numpy as np

CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits


def coordinates(line):
    return [float(line[i : i + 8]) for i in (30, 38, 46)]


def build(deposition, receptor, assembly_id, output):
    structure = gemmi.read_structure(str(deposition))
    if str(deposition).lower().endswith((".cif", ".mmcif")):
        raise ValueError("Use the PDB deposition with author-chain REMARK 350 records")
    assembly = next((a for a in structure.assemblies if a.name == assembly_id), None)
    if assembly is None:
        raise ValueError(f"assembly {assembly_id} missing in {deposition}")
    lines = receptor.read_text().splitlines()
    if any(l.startswith(("ROOT", "BRANCH", "BEGIN_RES", "MODEL")) for l in lines):
        raise ValueError("expected one rigid, asymmetric-unit receptor")
    if any(
        not l.startswith(("ATOM  ", "HETATM", "REMARK", "TER", "END")) and l.strip()
        for l in lines
    ):
        raise ValueError("unrecognized rigid receptor input record")
    originals = {}
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            originals.setdefault(line[21], []).append(line)
    if not originals:
        raise ValueError("empty receptor")
    # BIOMT applies only in deposition coordinates. Require common observed Cα positions.
    for chain, atoms in originals.items():
        source_chain = structure[0].find_chain(chain)
        if not source_chain:
            raise ValueError(f"chain {chain} absent in deposition")
        source_ca = {
            (r.seqid.num, r.seqid.icode.strip()): a.pos
            for r in source_chain
            for a in r
            if a.name == "CA" and a.altloc in ("\0", "A")
        }
        distances = [
            np.linalg.norm(np.asarray(coordinates(l)) - list(source_ca[k]))
            for l in atoms
            if l[12:16].strip() == "CA"
            for k in [(int(l[22:26]), l[26].strip())]
            if k in source_ca
        ]
        if not distances or np.max(distances) > 0.05:
            raise ValueError(
                f"chain {chain}: input and deposited Cα coordinates differ"
            )
    operations = []
    for group in assembly.generators:
        for source in group.chains:
            if source not in originals:
                raise ValueError(f"assembly chain {source} absent in receptor")
            for op in group.operators:
                matrix = np.asarray(op.transform.mat.tolist())
                offset = np.asarray(list(op.transform.vec))
                if not np.allclose(
                    matrix @ matrix.T, np.eye(3), atol=1e-4
                ) or not np.isclose(np.linalg.det(matrix), 1.0, atol=1e-4):
                    raise ValueError(f"invalid rigid transform for {source}")
                if any(
                    c == source
                    and np.allclose(matrix, m, atol=1e-6, rtol=0)
                    and np.allclose(offset, t, atol=1e-5, rtol=0)
                    for c, m, t in operations
                ):
                    continue
                operations.append((source, matrix, offset))
    if not operations or len(operations) > len(CHAIN_IDS):
        raise ValueError("unsupported number of PDBQT component chains")
    serial = 0
    out = []
    metadata = []
    identity = lambda m, t: (
        np.allclose(m, np.eye(3), atol=1e-6, rtol=0) and np.allclose(t, 0, atol=1e-6)
    )
    # Preserve source-chain names for identity copies and use unused characters for the rest.
    used = {c for c, m, t in operations if identity(m, t)}
    spare = iter(c for c in CHAIN_IDS if c not in used)
    for source, m, t in operations:
        dest = source if identity(m, t) else next(spare)
        metadata.append(
            dict(source=source, chain=dest, rotation=m.tolist(), translation=t.tolist())
        )
        for line in originals[source]:
            serial += 1
            if serial > 99999:
                raise ValueError("PDBQT atom serial overflow")
            xyz = m @ np.asarray(coordinates(line)) + t
            cols = "".join(f"{x:8.3f}" for x in xyz)
            if len(cols) != 24:
                raise ValueError("PDBQT coordinate column overflow")
            out.append(
                line[:6]
                + f"{serial:5d}"
                + line[11:21]
                + dest
                + line[22:30]
                + cols
                + line[54:]
            )
    output.write_text("\n".join(out) + "\n")
    return dict(
        assembly=assembly_id,
        description=assembly.oligomeric_details,
        copies=metadata,
        atoms=serial,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("deposition", type=Path)
    p.add_argument("receptor", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--assembly", required=True)
    p.add_argument("--metadata", type=Path)
    a = p.parse_args()
    result = build(a.deposition, a.receptor, a.assembly, a.output)
    if a.metadata:
        a.metadata.write_text(json.dumps(result, indent=2) + "\n")
