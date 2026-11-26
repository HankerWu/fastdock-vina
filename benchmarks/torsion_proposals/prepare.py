"""Prepare a reciprocal human crystal panel; use only molecular graph to initialize each query."""

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

PANELS = {
    "BACE1": {"2B8L": "5HA", "2ZJM": "F1M"},
    "EPHX2": {"3ANS": "S38", "1VJ5": "CIU", "3WKE": "AUB"},
    "CHEK1": {"2E9P": "77A", "2YWP": "A42"},
    "BRD9": {"6Y7L": "QMG", "6YQW": "82I"},
}


def main(repo, work):
    sys.path.insert(0, str(repo.resolve() / "benchmarks" / "paired_crossdocking"))
    import run as base

    raw = work / "raw"
    raw.mkdir(exist_ok=True, parents=True)

    def prepare_one(entry):
        target, pdb, ligand = entry
        out = work / pdb
        out.mkdir(exist_ok=True)
        row = {"target": target, "pdb": pdb, "ligand": ligand}
        try:
            pdbfile = raw / (pdb + ".pdb")
            ccd = raw / (ligand + ".cif")
            row["pdb_sha256"] = base.fetch(
                pdbfile, f"https://files.rcsb.org/download/{pdb}.pdb"
            )
            row["ccd_sha256"] = base.fetch(
                ccd, f"https://files.rcsb.org/ligands/download/{ligand}.cif"
            )
            ls = pdbfile.read_text().splitlines()
            refs = base.dep_records(ls, ligand)
            if not refs[0]:
                raise ValueError("no query in chain A")
            obs = np.array(list(refs[0].values()))
            center = obs.mean(axis=0)
            excluded = [
                (
                    np.linalg.norm(obs - base.coords(line), axis=1).min(),
                    line[17:27].strip(),
                )
                for line in ls
                if line.startswith(("HETATM", "ATOM  "))
                and line[17:20] != "HOH"
                and (
                    line[21] != "A"
                    or (line.startswith("HETATM") and line[17:20].strip() != ligand)
                )
                and line[76:78].strip() not in ("H", "D")
            ]
            row["closest_omission"] = min(excluded) if excluded else None
            if excluded and min(excluded)[0] < 8:
                raise ValueError("nonwater omitted within 8 A")
            lig = base.ligand(ccd, refs, out / "ligand.pdbqt")
            rec = base.receptor(ls, out, pdb, obs)
            ca = {
                line[22:27].strip() + ":" + line[17:20]: base.coords(line)
                for line in ls
                if line.startswith("ATOM  ")
                and line[21] == "A"
                and line[12:16].strip() == "CA"
                and line[16] in " A"
            }
            row.update(lig, **rec, ca=ca, center=center.tolist(), status="ready")
        except Exception as e:
            row.update(status="excluded", exception=repr(e))
        print(
            pdb,
            row["status"],
            row.get("exception", ""),
            row.get("closest_omission"),
            flush=True,
        )
        (out / "preparation.json").write_text(json.dumps(row, indent=2) + "\n")
        return row

    cases = [
        (target, pdb, lig)
        for target, entries in PANELS.items()
        for pdb, lig in entries.items()
    ]
    with ThreadPoolExecutor(3) as pool:
        rows = list(pool.map(prepare_one, cases))
    meta = {r["pdb"]: r for r in rows if r["status"] == "ready"}
    if "1VJ5" in meta:
        meta.pop("3WKE", None)
    for target in PANELS:
        if sum(r["target"] == target for r in meta.values()) != 2:
            raise ValueError(
                f"Missing structural pair for {target}; inspect preparation.json"
            )
    (work / "metadata.json").write_text(json.dumps(meta) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("repo", type=Path)
    p.add_argument("work", type=Path)
    a = p.parse_args()
    main(a.repo, a.work)
