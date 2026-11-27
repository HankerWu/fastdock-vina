"""Build and evaluate the preregistered ATAD2 cognate and cross-docking panel."""

import argparse
import json
import runpy
import sys
from hashlib import sha256
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
sys.path.insert(0, str(BENCH / "ligand_hydroxyl/transfer"))
import run_transfer

run_transfer.GROUPS = {"ATAD2": {"5A82": "YEJ", "6YB4": "OJH"}}
from check_transfer import verify


def prepare(work):
    raw = work / "raw"
    raw.mkdir(exist_ok=True, parents=True)
    source = "https://files.rcsb.org/download/6YB4.cif"
    h = run_transfer.cross.fetch(raw / "6YB4.cif", source)
    structure = gemmi.read_structure(str(raw / "6YB4.cif"))
    if len(structure) != 1 or [chain.name for chain in structure[0]] != ["AAA"]:
        raise ValueError("Unexpected 6YB4 model or deposition chain")
    # The unique author chain exceeds the one-column PDB field.
    structure[0][0].name = "A"
    text = structure.make_pdb_string()
    if "REMARK   2 RESOLUTION." not in text:
        text = (
            f"REMARK   2 RESOLUTION.    {structure.resolution:.2f} ANGSTROMS.\n" + text
        )
    (raw / "6YB4.pdb").write_text(text)
    run_transfer.prepare(work)
    meta = json.loads((work / "metadata.json").read_text())
    sources = meta["6YB4"]["sources"]
    generated = sources.pop("https://files.rcsb.org/download/6YB4.pdb")
    sources[source] = h
    meta["6YB4"]["gemmi_generated_pdb_sha256"] = generated
    meta["6YB4"]["author_chain_to_pdb"] = {"AAA": "A"}
    (work / "metadata.json").write_text(json.dumps(meta, sort_keys=True) + "\n")


def score(work, report):
    sys.path.insert(0, str(HERE.parent))
    from coupled import CoupledAmides
    from evaluate import ARMS, quantiles
    from model import AmideReceptor

    count = runpy.run_path(str(BENCH / "terminal_rotor_score/independent/verify.py"))[
        "independently_count"
    ]
    geom = pd.read_csv(work / "geometry.csv")
    report.mkdir(exist_ok=True, parents=True)
    sites, poses, cells = [], [], []
    metadata = json.loads((work / "metadata.json").read_text())
    receptors = {}
    for (q, r, cap, seed, engine), frame in geom.groupby(
        ["query", "receiver", "cap", "seed", "engine"]
    ):
        key = {
            "target": "ATAD2",
            "query": q,
            "receiver": r,
            "kind": "self" if q == r else "cross",
            "cap": cap,
            "seed": seed,
            "engine": engine,
        }
        if r not in receptors:
            rec = AmideReceptor((work / r / "receptor.pdbqt").read_text())
            receptors[r] = (rec, CoupledAmides(rec))
        rec, joint = receptors[r]
        query = (work / q / "ligand.pdbqt").read_text()
        scale = 1 / (1 + 0.05846 * count(query.encode()))
        blocks = [
            b
            for b in (work / "runs" / f"{q}_{r}_{cap}_{seed}_{engine}.pdbqt")
            .read_text()
            .split("ENDMDL")
            if b.strip()
        ]
        pool = []
        for rr, b in zip(frame.sort_values("model").itertuples(), blocks, strict=True):
            p = dict(
                **key,
                model=rr.model,
                rmsd_A=rr.rmsd,
                native=rr.dry,
                atom_record_sha256=sha256(
                    "\n".join(
                        l for l in b.splitlines() if l.startswith(("ATOM  ", "HETATM"))
                    ).encode()
                ).hexdigest(),
                reference_unbound=float(
                    next(
                        l for l in b.splitlines() if l.startswith("REMARK UNBOUND:")
                    ).split()[-1]
                ),
                inter_plus_intra=float(
                    next(
                        l
                        for l in b.splitlines()
                        if l.startswith("REMARK INTER + INTRA:")
                    ).split()[-1]
                ),
                n_ca=rr.n_ca,
                receptor_fit_A=rr.ca_rmsd,
            )
            if int(b.split("MODEL")[1].split()[0]) != rr.model:
                raise ValueError("model identity")
            if engine == "pro":
                terms = rec.score(query, b, detail=True)
                terms["joint"] = joint.correction(terms["sites"])
                p["scale"] = scale
                for site, apo, dl, change in terms["sites"]:
                    sites.append(
                        dict(
                            **key,
                            model=rr.model,
                            residue=site,
                            apo_delta=apo,
                            ligand_delta=dl,
                            ensemble_change=change,
                        )
                    )
                for arm in ARMS[1:] + ["joint"]:
                    p[arm] = rr.dry + scale * terms[arm]
            poses.append(p)
            pool.append(p)
        for arm in ARMS + ["joint"] if engine == "pro" else ["native"]:
            ranked = (
                sorted(pool, key=lambda p: (p[arm], p["model"]))
                if arm != "native"
                else pool
            )
            rms = [p["rmsd_A"] for p in ranked]
            cells.append(
                dict(
                    **key,
                    arm=arm,
                    selected=ranked[0]["model"],
                    top1_A=rms[0],
                    best5_A=min(rms[:5]),
                    oracle_A=min(rms),
                    regret_A=rms[0] - min(rms),
                    count=len(rms),
                    delta_A=rms[0] - pool[0]["rmsd_A"],
                )
            )
    for name, data in [("poses", poses), ("cells", cells), ("sites", sites)]:
        pd.DataFrame(data).to_csv(report / f"{name}.csv", index=False)
    metrics = []
    for (engine, arm), g in pd.DataFrame(cells).groupby(["engine", "arm"]):
        for field in ["top1_A", "best5_A", "regret_A", "oracle_A", "delta_A"]:
            metrics.append(
                dict(
                    engine=engine, arm=arm, field=field, n=len(g), **quantiles(g[field])
                )
            )
    pd.DataFrame(metrics).to_csv(report / "summary.csv", index=False)
    # Diagnostic heavy coordinates from the deposition. Hydrogens only set donor
    # topology here. The potential does not use their coordinates. Never rank with these rows.
    diagnostics = []
    for query_id, (rec, _) in receptors.items():
        ligand = (work / query_id / "ligand.pdbqt").read_text()
        for mapping, coordinates in enumerate(metadata[query_id]["mapped"]):
            by_serial = dict(
                zip(metadata[query_id]["serial"], coordinates, strict=True)
            )
            lines = []
            for line in ligand.splitlines():
                if (
                    line.startswith(("ATOM  ", "HETATM"))
                    and int(line[6:11]) in by_serial
                ):
                    line = (
                        line[:30]
                        + "".join(f"{v:8.3f}" for v in by_serial[int(line[6:11])])
                        + line[54:]
                    )
                lines.append(line)
            terms = rec.score(ligand, "\n".join(lines), detail=True)
            for resid, apo, dl, dc in terms["sites"]:
                endpoints = next(
                    site.endpoints for site in rec.sites if site.residue == resid
                )
                d = np.linalg.norm(
                    endpoints[:, None] - np.asarray(coordinates)[None, :], axis=2
                )
                diagnostics.append(
                    dict(
                        query=query_id,
                        mapping=mapping,
                        residue=resid,
                        min_to_O_A=d[0].min(),
                        min_to_N_A=d[1].min(),
                        apo_delta=apo,
                        ligand_delta=dl,
                        ensemble_correction=dc,
                    )
                )
    pd.DataFrame(diagnostics).to_csv(report / "crystal_sites.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("step", choices=["prepare", "dock", "check", "score"])
    p.add_argument("work", type=Path)
    p.add_argument("--pro", type=Path)
    p.add_argument("--upstream", type=Path)
    p.add_argument("--report", type=Path)
    a = p.parse_args()
    w = a.work.resolve()
    if a.step == "prepare":
        prepare(w)
    elif a.step == "dock":
        run_transfer.dock(
            w, {"pro": a.pro.resolve(), "upstream": a.upstream.resolve()}, 2
        )
    elif a.step == "check":
        verify(w)
    else:
        score(w, a.report)
