"""Assemble the two mature CatD protein fragments using the deposited A/B coordinates."""

import argparse
import json
from pathlib import Path
import sys
import time
from urllib.request import urlretrieve
import numpy as np
from prepare_transfer import CASES, REPO

sys.path.insert(0, str(REPO / "benchmarks/paired_crossdocking"))
import run as base


def prepare(work, pdb):
    started = time.perf_counter()
    d = work / pdb
    meta = json.loads((d / "preparation.json").read_text())
    lines = (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
    xyz = np.array(meta["mapped"][0])
    _, ccd, chain = CASES[pdb]
    ligid = next(
        l[22:26].strip()
        for l in lines
        if l.startswith("HETATM") and l[17:20] == ccd and l[21] == chain
    )
    omitted = [
        l
        for l in lines
        if l.startswith(("ATOM  ", "HETATM"))
        and l[76:78].strip() not in ("H", "D")
        and l[17:20].strip() != "HOH"
        and not (l.startswith("ATOM  ") and l[21] in ("A", "B"))
        and not (l[17:20] == ccd and l[21] == chain and l[22:26].strip() == ligid)
    ]
    nearest = min(
        (float(np.linalg.norm(xyz - base.coords(l), axis=1).min()), l[12:27].strip())
        for l in omitted
    )
    record = {
        "receptor_chains": ["A", "B"],
        "nearest_omitted_nonwater": nearest,
        "status": "omitted nonwater atom within 8 A",
    }
    if nearest[0] >= 8:
        reports, atoms = {}, []
        try:
            for ch in ("A", "B"):
                group = [l for l in lines if l.startswith("ATOM  ") and l[21] == ch]
                ordinals = dict(
                    zip(
                        dict.fromkeys(l[22:27] for l in group), range(1, len(group) + 1)
                    )
                )
                transformed = [
                    l[:21] + "A" + f"{ordinals[l[22:27]]:4d} " + l[27:] for l in group
                ]
                sub = d / ch
                sub.mkdir(exist_ok=True)
                report = base.receptor(transformed, sub, pdb, xyz)
                near = report["closest_repaired_ligand_distance"]
                if near is not None and near < 8:
                    raise ValueError(f"Repaired {ch} atom {near:.3f} A from ligand")
                reports[ch] = report
                for l in (sub / "receptor.pdbqt").read_text().splitlines():
                    if l.startswith(("ATOM  ", "HETATM")):
                        atoms.append(
                            l[:6] + f"{len(atoms) + 1:5d}" + l[11:21] + ch + l[22:]
                        )
            (d / "receptor.pdbqt").write_text("\n".join(atoms) + "\n")
            sdf = work / "raw" / (pdb + "-bound.sdf")
            if not sdf.exists():
                urlretrieve(
                    f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id=A&auth_seq_id={ligid}&encoding=sdf&copy_all_categories=false",
                    sdf,
                )
            record.update(status="ready", reports=reports)
        except (ValueError, RuntimeError, base.subprocess.CalledProcessError) as e:
            record["status"] = str(e)
    record["extra_preparation_s"] = time.perf_counter() - started
    (d / "assembly.json").write_text(json.dumps(record, indent=2) + "\n")
    print(pdb, record["status"], nearest, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    work = p.parse_args().work.resolve()
    for pdb in ("6QBG", "6QCB", "6QBH"):
        prepare(work, pdb)
