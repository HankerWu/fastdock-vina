"""Prepare and dock the two prospective panels, or rescore archived candidate pools."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import runpy
import sys
import tarfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
NEW = {
    "3SMR": ("WDR5", "NP7", "A"),
    "4EY6": ("ACHE", "GNT", "A"),
    "4EY7": ("ACHE", "E20", "A"),
    "4NYX": ("CREBBP", "2O4", "A"),
    "4RVR": ("BAZ2B", "3WQ", "A"),
    "4J0S": ("BRD4", "1H3", "A"),
    "4QB3": ("BRD4", "30M", "A"),
}
ARCHIVES = (
    "placement-confirm/dpp4_inputs.tar.xz",
    "placement-holdout/inputs.tar.xz",
    "tangent/inputs.tar.xz",
)
OLD = ("1NAX", "3IMY", "5EU1", "5F1H", "1P8D", "1UPV", "3VJK", "4A5S")
SEEDS = (71993, 72019)
ARMS = ("native", "ensemble", "minimum", "uniform", "flip_all", "coupled")


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    with path.open("w", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def unpack(path, folder):
    with tarfile.open(path) as tar:
        tar.extractall(folder, filter="data")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("restore", "prepare", "dock", "score"))
    p.add_argument("work", type=Path)
    p.add_argument("--repo", type=Path, default=HERE.parents[2])
    p.add_argument("--model-dir", type=Path, default=HERE.parent)
    p.add_argument("--output", type=Path, default=HERE / "results")
    p.add_argument("--binary", type=Path)
    a = p.parse_args()
    work = a.work.resolve()
    bench = a.repo.resolve() / "benchmarks"
    sys.path.insert(0, str(bench / "search_backtracking/pose"))
    import backtracking_panel as panel
    from experiment import geometry
    import experiment

    sys.path.insert(0, str(a.model_dir.resolve()))
    if a.action == "score":
        from model import AmideReceptor
        from coupled import CoupledAmides
    count = runpy.run_path(str(bench / "terminal_rotor_score/independent/verify.py"))[
        "independently_count"
    ]
    out = a.output
    out.mkdir(parents=True, exist_ok=True)
    all_cells, all_poses, coverage, sites, costs, summary, closures, ecdf = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for study, ids in [("historical", OLD), ("novel", tuple(NEW))]:
        d = work / study
        (d / "runs").mkdir(parents=True, exist_ok=True)
        if a.action == "restore":
            if study == "historical":
                for name in ARCHIVES:
                    unpack(bench / "hunt_overlap/pose" / name, d)
                unpack(HERE / "transfer_evidence.tar.xz", work)
            continue
        if a.action == "prepare":
            if study == "historical":
                for name in ARCHIVES:
                    unpack(bench / "hunt_overlap/pose" / name, d)
            else:
                panel.CASES = NEW
                for q in NEW:
                    panel.prepare_input(d, q)
            continue
        meta = {q: json.loads((d / q / "preparation.json").read_text()) for q in ids}
        for q, info in meta.items():
            source = json.loads((d / q / "source.json").read_text())
            row = {
                k: source[k]
                for k in (
                    "pdb",
                    "target",
                    "ccd",
                    "chain",
                    "resolution_A",
                    "formal_charge",
                    "structure_source",
                    "ccd_source",
                    "mate_heavy_pairs_8A",
                    "closest_mate_pair",
                )
            }
            nearest = info["closest_omitted_heavy_atom"]
            row.update(
                pdb_sha256=info["pdb_sha256"],
                ccd_sha256=info["ccd_sha256"],
                study=study,
                status=info["receptor_status"],
                closest_omitted=nearest,
                altlocs=info["n_altloc"],
                mappings=len(info["mapped"]),
            )
            coverage.append(row)
            if (
                info["receptor_status"] == "ready"
                and not (d / q / "receptor.pdbqt").is_file()
            ):
                raise FileNotFoundError(d / q / "receptor.pdbqt")
        experiment.CASES = ids
        pairs = geometry(d, meta)
        if a.action == "dock":
            if not a.binary:
                p.error("dock needs --binary")
            cells, poses = [], []
            for q, r in pairs:
                for seed in SEEDS:
                    row, modes = panel.paired.execute(
                        d, meta, q, r, seed, 12000, "baseline", a.binary, False
                    )
                    cells.append(row)
                    poses.extend(modes)
            write(d / "baseline_cells.csv", cells)
            write(d / "baseline_modes.csv", poses)
            continue
        groups = defaultdict(list)
        for row in read(d / "baseline_modes.csv"):
            groups[row["key"]].append(row)
        receptor = {}
        coupled = {}
        scales = {
            q: 1 / (1 + 0.05846 * count((d / q / "ligand.pdbqt").read_bytes()))
            for q in ids
        }
        source_cells = read(d / "baseline_cells.csv")
        requested = {(q, r, s) for q, r in pairs for s in SEEDS}
        if {
            (c["query"], c["receiver"], int(c["seed"])) for c in source_cells
        } != requested or len(requested) != len(source_cells):
            raise ValueError("Docked cells differ from design")
        for cell in source_cells:
            q, r = cell["query"], cell["receiver"]
            h = {
                k: cell[k]
                for k in ("key", "target", "query", "receiver", "seed", "cap", "kind")
            }
            h.update(study=study)
            if r not in receptor:
                t = time.perf_counter()
                receptor[r] = AmideReceptor((d / r / "receptor.pdbqt").read_text())
                costs.append(
                    dict(
                        study=study,
                        receiver=r,
                        key="",
                        stage="receptor",
                        seconds=time.perf_counter() - t,
                        n=1,
                    )
                )
                t = time.perf_counter()
                coupled[r] = CoupledAmides(receptor[r])
                costs.append(
                    dict(
                        study=study,
                        receiver=r,
                        key="",
                        stage="coupled_receptor",
                        seconds=time.perf_counter() - t,
                        n=1,
                    )
                )
            ligand = (d / q / "ligand.pdbqt").read_text()
            pool = []
            archive = (d / "runs" / (cell["key"] + ".pdbqt")).read_bytes()
            digest = hashlib.sha256(archive).hexdigest()
            elapsed = 0.0
            coupled_elapsed = 0.0
            for block, g in zip(
                archive.decode().split("ENDMDL")[:-1],
                groups.pop(cell["key"]),
                strict=True,
            ):
                t = time.perf_counter()
                terms = receptor[r].score(ligand, block, detail=True)
                elapsed += time.perf_counter() - t
                t = time.perf_counter()
                terms["coupled"] = coupled[r].correction(terms["sites"])
                coupled_elapsed += time.perf_counter() - t
                rank = int(g["mode"])
                if q == "4EY6":
                    # Both physical break carbons are in the original graph. G0 sites are not.
                    from rdkit import Chem

                    cg = [
                        int(line[6:11])
                        for line in block.splitlines()
                        if line.startswith(("ATOM  ", "HETATM"))
                        and line.split()[-1] == "CG0"
                    ]
                    if len(cg) != 2:
                        raise ValueError("Expected one galanthamine closure bond")
                    serial = meta[q]["serial"]
                    i, j = [serial.index(i) for i in cg]
                    xyz = panel.paired.atoms(block.splitlines())
                    mol = Chem.MolFromSmiles(meta[q]["canonical_smiles"])
                    if mol.GetBondBetweenAtoms(i, j) is None:
                        raise ValueError(
                            "Physical break sites are not bonded in the ligand graph"
                        )
                    coords = np.asarray([xyz[k][0] for k in serial])
                    ref = np.asarray(meta[q]["mapped"])
                    refdist = np.linalg.norm(ref[:, i] - ref[:, j], axis=1)
                    conf = Chem.Conformer(len(serial))
                    for n, point in enumerate(coords):
                        conf.SetAtomPosition(n, point)
                    observed = Chem.Mol(mol)
                    observed.AddConformer(conf)
                    Chem.AssignStereochemistryFrom3D(observed)
                    assigned = dict(Chem.FindMolChiralCenters(observed))
                    specified = Chem.FindMolChiralCenters(mol)
                    mismatch = sum(assigned.get(n) != label for n, label in specified)
                    closures.append(
                        dict(
                            **h,
                            mode=rank,
                            carbon_1=cg[0],
                            carbon_2=cg[1],
                            predicted_A=float(np.linalg.norm(coords[i] - coords[j])),
                            reference_min_A=float(min(refdist)),
                            reference_max_A=float(max(refdist)),
                            specified_centers=len(specified),
                            stereo_mismatches=mismatch,
                        )
                    )
                native = float(g["score"])
                scale = scales[q]
                remarks = {
                    l.split(":", 1)[0]: float(l.split(":", 1)[1].split()[0])
                    for l in block.splitlines()
                    if l.startswith(("REMARK INTER + INTRA:", "REMARK UNBOUND:"))
                }
                reconstructed = scale * (
                    remarks["REMARK INTER + INTRA"] - remarks["REMARK UNBOUND"]
                )
                if abs(native - reconstructed) > 0.0016:
                    raise ValueError(
                        f"Native normalization does not reproduce printed terms: {h}, {rank}, {native}, {reconstructed}"
                    )
                pose = dict(
                    **h,
                    mode=rank,
                    rmsd_A=float(g["rmsd_A"]),
                    clashes=int(g["clash_pairs"]),
                    native=native,
                    scale=scale,
                    active_sites=terms["active_sites"],
                )
                for arm in ARMS[1:]:
                    pose[arm] = native + scale * terms[arm]
                pool.append(pose)
                for name, apo, dl, c in terms["sites"]:
                    sites.append(
                        dict(
                            **h,
                            mode=rank,
                            site=name,
                            apo_delta=apo,
                            ligand_delta=dl,
                            correction=c,
                        )
                    )
            costs.append(
                dict(
                    study=study,
                    receiver=r,
                    key=cell["key"],
                    stage="rescore",
                    seconds=elapsed,
                    n=len(pool),
                )
            )
            costs.append(
                dict(
                    study=study,
                    receiver=r,
                    key=cell["key"],
                    stage="additional_coupling",
                    seconds=coupled_elapsed,
                    n=len(pool),
                )
            )
            all_poses.extend(pool)
            for arm in ARMS:
                order = sorted(pool, key=lambda row: (row[arm], row["mode"]))
                best = order[0]
                oracle = min(z["rmsd_A"] for z in pool)
                native1 = pool[0]["rmsd_A"]
                all_cells.append(
                    dict(
                        **h,
                        method=arm,
                        n=len(pool),
                        selected=best["mode"],
                        top1_A=best["rmsd_A"],
                        delta_A=best["rmsd_A"] - native1,
                        runner_up_gap=order[1][arm] - best[arm]
                        if len(order) > 1
                        else "",
                        reordered_positions=sum(
                            z["mode"] != i for i, z in enumerate(order, 1)
                        ),
                        best5_A=min(z["rmsd_A"] for z in order[:5]),
                        oracle_A=oracle,
                        regret_A=best["rmsd_A"] - oracle,
                        clashes=best["clashes"],
                        ca_count=cell["ca_count"],
                        ca_rmsd_A=cell["ca_rmsd_A"],
                        docking_s=cell["seconds"],
                        hunt_evals=cell["hunt_evals"],
                        second_evals=cell["second_evals"],
                        explicit_evals=cell["explicit_evals"],
                        pool_sha256=digest,
                    )
                )
        if groups:
            raise ValueError("Unmatched candidate coordinates")
    if a.action != "score":
        return
    for study in ("historical", "novel"):
        studyrows = [z for z in all_cells if z["study"] == study]
        for arm in ARMS:
            armrows = [z for z in studyrows if z["method"] == arm]
            for group in ["all", "cognate", "cross"] + sorted(
                {c["target"] for c in armrows}
            ):
                group_rows = [
                    z for z in armrows if group in ("all", z["kind"], z["target"])
                ]
                if not group_rows:
                    continue
                h = dict(
                    study=study,
                    method=arm,
                    group=group,
                    n=len(group_rows),
                    changed=sum(z["selected"] != 1 for z in group_rows),
                )
                for field in ["top1_A", "delta_A", "best5_A", "oracle_A", "regret_A"]:
                    x = [z[field] for z in group_rows]
                    q = np.quantile(
                        x, [0.25, 0.5, 0.75, 0.9, 0.95], method="averaged_inverted_cdf"
                    )
                    h.update(
                        {
                            field + "_" + k: v
                            for k, v in zip(
                                ("mean", "q25", "median", "q75", "p90", "p95"),
                                (np.mean(x), *q),
                            )
                        }
                    )
                    if group == "all" and arm in ("native", "ensemble"):
                        values, nobs = np.unique(x, return_counts=True)
                        for value, frac in zip(values, np.cumsum(nobs) / len(x)):
                            ecdf.append(
                                dict(
                                    study=study,
                                    method=arm,
                                    metric=field,
                                    error_A=value,
                                    fraction=frac,
                                )
                            )
                summary.append(h)
    for name, records in [
        ("cells", all_cells),
        ("modes", all_poses),
        ("receptor_coverage", coverage),
        ("site_terms", sites),
        ("cost", costs),
        ("summary", summary),
        ("galanthamine_geometry", closures),
        ("ecdf", ecdf),
    ]:
        write(out / (name + ".csv"), records)


if __name__ == "__main__":
    main()
