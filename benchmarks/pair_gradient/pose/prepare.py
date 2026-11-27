"""Prepare frozen ADORA2A and renin cases using the historical non-crystal ligand recipe."""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "benchmarks/paired_crossdocking"))
import run as base

CASES = {
    "3EML": ("ADORA2A", "ZMA", "A"),
    "3UZA": ("ADORA2A", "T4G", "A"),
    "1RNE": ("REN", "C60", "A"),
    "2V0Z": ("REN", "C41", "C"),
    "3G6Z": ("REN", "A7T", "A"),
    "3SFC": ("REN", "S53", "B"),
}


def prepare(work, pdb):
    target, component, chain = CASES[pdb]
    out = work / pdb
    raw = work / "raw"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    structure, ccd = raw / (pdb + ".pdb"), raw / (component + ".cif")
    meta = {"target": target, "pdb": pdb, "ligand": component, "chain": chain}
    meta["pdb_sha256"] = base.fetch(
        structure, f"https://files.rcsb.org/download/{pdb}.pdb"
    )
    meta["ccd_sha256"] = base.fetch(
        ccd, f"https://files.rcsb.org/ligands/download/{component}.cif"
    )
    lines = structure.read_text().splitlines()
    selected = [
        line[:21] + "A" + line[22:]
        for line in lines
        if line.startswith(("ATOM  ", "HETATM")) and line[21] == chain
    ]
    refs = base.dep_records(selected, component)
    crystal = np.asarray(list(refs[0].values()))
    meta.update(base.ligand(ccd, refs, out / "ligand.pdbqt"))
    meta["center"] = crystal.mean(axis=0).tolist()
    meta["ca"] = {
        line[22:27].strip() + ":" + line[17:20]: base.coords(line)
        for line in selected
        if line.startswith("ATOM  ")
        and line[16] in " A"
        and line[12:16].strip() == "CA"
    }
    omitted = [
        (
            float(np.linalg.norm(crystal - base.coords(line), axis=1).min()),
            line[12:27].strip(),
        )
        for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
        and line[76:78].strip() not in {"H", "D"}
        and line[17:20].strip() != "HOH"
        and not (
            line[21] == chain
            and (line.startswith("ATOM  ") or line[17:20].strip() == component)
        )
    ]
    meta["closest_omitted_heavy_atom"] = min(omitted) if omitted else None
    if omitted and min(omitted)[0] < 8:
        meta["receptor_status"] = "omitted nonwater atom within 8 A"
    else:
        try:
            # The historical Meeko reader can merge insertion-coded residues.
            # Use a bijective ordinal ID in preparation only, preserving original
            # deposition numbering for cross-receptor backbone correspondence.
            residue_ids = dict.fromkeys(
                line[22:27] for line in selected if line.startswith("ATOM  ")
            )
            ordinal = {key: i for i, key in enumerate(residue_ids, 1)}
            rec_lines = [
                line[:22] + f"{ordinal[line[22:27]]:4d} " + line[27:]
                if line.startswith("ATOM  ")
                else line
                for line in selected
            ]
            meta["receptor_ordinal_to_deposition"] = list(residue_ids)
            meta.update(base.receptor(rec_lines, out, pdb, crystal))
            distance = meta["closest_repaired_ligand_distance"]
            meta["receptor_status"] = (
                "ready"
                if distance is None or distance >= 8
                else "rebuilt atom within 8 A"
            )
        except (RuntimeError, ValueError, OSError) as exc:
            meta["receptor_status"] = f"receptor preparation failed: {exc}"
        except base.subprocess.CalledProcessError:
            meta["receptor_status"] = "Meeko receptor preparation failed; see prep.log"
    (out / "preparation.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(pdb, meta["receptor_status"], meta["closest_omitted_heavy_atom"], flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("ids", nargs="*", default=list(CASES))
    args = p.parse_args()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    for code in args.ids:
        prepare(args.work, code)
