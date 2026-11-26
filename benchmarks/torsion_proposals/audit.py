"""Verify each receptor-frame displacement via SIFTS, independent graph matching and SciPy alignment."""

import argparse
import csv
import json
import sys
from pathlib import Path
from urllib.request import urlopen

import gemmi
import numpy as np


def mappings(work, code):
    file = work / "raw" / (code + "_pdbe.cif")
    if not file.exists():
        url = f"https://www.ebi.ac.uk/pdbe/coordinates/files/{code.lower()}_updated.cif"
        file.write_bytes(urlopen(url, timeout=60).read())
    block = gemmi.cif.read_file(str(file)).sole_block()
    fields = (
        "group_PDB",
        "auth_asym_id",
        "label_atom_id",
        "label_alt_id",
        "pdbx_sifts_xref_db_acc",
        "pdbx_sifts_xref_db_num",
        "label_comp_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
    )
    ca = {}
    for group, chain, atom, alt, accession, num, res, x, y, z in block.find(
        "_atom_site.", fields
    ):
        if (
            group == "ATOM"
            and chain == "A"
            and atom == "CA"
            and alt in (".", "A")
            and num not in ("?", ".")
        ):
            ca[(accession, num, res)] = np.array([float(x), float(y), float(z)])
    if not ca:
        raise ValueError("empty SIFTS")
    return ca


def audit(repo, work):
    sys.path.insert(0, str(repo.resolve() / "benchmarks/paired_crossdocking"))
    import independent_cross as ind

    md = json.loads((work / "metadata.json").read_text())
    ind.COMP = {c: v["ligand"] for c, v in md.items()}
    data = {}
    for code in ind.COMP:
        _, center, refs, serials, extra = ind.complex_data(work, code)
        data[code] = (mappings(work, code), center, refs, serials, extra)
    # These two BRD9 depositions use different UniProt isoform accessions.
    # Match contiguous SIFTS-numbered residue strings; all unique shared
    # seven-residue stretches must imply the same translation.
    first, second = (data[c][0] for c in ("6Y7L", "6YQW"))

    def sevenmers(ca):
        sequence = {int(n): r for (_, n, r) in ca}
        windows = {}
        for i in sequence:
            if all(j in sequence for j in range(i, i + 7)):
                w = tuple(sequence[j] for j in range(i, i + 7))
                windows.setdefault(w, []).append(i)
        return windows

    w1, w2 = sevenmers(first), sevenmers(second)
    offsets = {
        w2[k][0] - w1[k][0]
        for k in w1.keys() & w2.keys()
        if len(w1[k]) == len(w2[k]) == 1
    }
    if len(offsets) != 1:
        raise ValueError("ambiguous BRD9 isoform translation")
    offset = offsets.pop()
    data["6Y7L"] = (
        {("Q9H8M2", str(int(n) + offset), r): v for (_, n, r), v in first.items()},
        *data["6Y7L"][1:],
    )
    groups = {}
    for q in data:
        for r in data:
            if md[q]["target"] == md[r]["target"]:
                groups[q, r] = ind.fitted(data[q], data[r])
    checks = []
    for f in sorted(work.glob("modes-*.csv")):
        rows = list(csv.DictReader(f.open()))
        counts = {}
        max_rms = max_ca = 0.0
        for row in rows:
            q, r, cap, arm, seed = [
                row[k] for k in ("ligand", "receiver", "cap", "arm", "seed")
            ]
            path = work / f"runs_{cap}" / f"{q}_{r}_{seed}_{arm}.pdbqt"
            if path not in counts:
                blocks = [
                    s
                    for s in path.read_text().split("ENDMDL")
                    if "REMARK VINA RESULT:" in s
                ]
                counts[path] = blocks
            block = counts[path][int(row["rank"]) - 1]
            atom = {
                int(line[6:11]): ind.xyz(line)
                for line in block.splitlines()
                if line.startswith(("ATOM", "HETATM"))
            }
            points = np.array([atom[i] for i in data[q][3]])
            refs, _, n_ca, backbone = groups[q, r]
            rms = float(
                np.sqrt(np.mean(np.sum((refs - points) ** 2, axis=2), axis=1)).min()
            )
            max_rms = max(max_rms, abs(float(row["rmsd"]) - rms))
            max_ca = max(max_ca, abs(float(row["ca_rmsd"]) - backbone))
            if n_ca != int(row["n_ca"]):
                raise ValueError((q, r, "mismatched SIFTS coverage", n_ca, row["n_ca"]))
        checks.append(
            dict(
                table=f.name,
                cap=cap,
                arm=arm,
                runs=len(counts),
                modes=len(rows),
                max_rmsd_diff=max_rms,
                max_ca_diff=max_ca,
            )
        )
    (work / "audit.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(json.dumps(checks, indent=2))
    if any(v["max_rmsd_diff"] > 1e-5 for v in checks):
        raise ValueError("RMSD mismatch")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    audit(a.repo, a.work)
