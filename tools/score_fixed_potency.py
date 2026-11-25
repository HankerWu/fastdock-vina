"""Freeze physical scores on saved native top poses before reading matched Ki values."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
import json
from math import isfinite
from pathlib import Path
import re
import subprocess
from time import perf_counter

from rdkit import Chem
from scipy.stats import spearmanr
import yaml

from explicit_rbfe import configuration, digest, ligand_molecules, now, put_csv, put_json
from score_physical_pool import ARMS, decode_ligand, read_poses, score_pool


def build_cases(sources):
    molecules = ligand_molecules(sources["historical_ligands"])
    measurements = yaml.safe_load(sources["measurements"].read_text())
    if set(molecules) != set(measurements):
        raise ValueError("The ligand identity and measurement populations differ")
    for name, molecule in molecules.items():
        assay = measurements[name]["measurement"]
        config = sources["config"]
        if (assay["type"] != "ki" or assay["unit"] != "uM" or assay.get("relation", "=") != "="
                or assay["doi"] != config["measurement_doi"]
                or not assay["comment"].startswith(config["measurement_locator"])
                or Chem.MolToSmiles(molecule) != Chem.MolToSmiles(Chem.MolFromSmiles(measurements[name]["smiles"]))):
            raise ValueError(f"{name}: unexpected identity or assay metadata")
        score_record = sources["docked"] / name / "complete.json"
        yield dict(name=name, identity=Chem.MolToSmiles(molecule), reference=assay["comment"],
                   source_docking_record=str(score_record))


def freeze_pose(case, directory):
    dock = json.loads(Path(case["source_docking_record"]).read_text())
    command = dock["command"]
    receptor = Path(command[command.index("--receptor") + 1]).resolve(strict=True)
    pool = Path(dock["poses"]).resolve(strict=True)
    poses = read_poses(pool)
    [native] = [pose for pose in poses if pose.native_rank == 1]
    if native.native_intra is None or native.native_intra != native.native_unbound or any(
            p.native_unbound != native.native_unbound for p in poses):
        raise ValueError(f"{pool}: expected one shared intra reference matching native rank 1")
    _, identity = decode_ligand(native)
    if identity != case["identity"]:
        raise ValueError(f"{pool}: pose identity does not match the experimental series")
    blocks = re.findall(r"^MODEL\s+(\d+)\s*\n(.*?)^ENDMDL\s*$", pool.read_text(), flags=re.M | re.S)
    [block] = [text for model, text in blocks if int(model) == native.model]
    fixed = directory / (case["name"] + ".pdbqt")
    with fixed.open("x") as stream:
        stream.write(f"MODEL 1\n{block}ENDMDL\n")
    [check] = read_poses(fixed)
    if (check.native_score, check.native_rank, check.text, check.native_intra, check.native_unbound) != (
        native.native_score, native.native_rank, native.text, native.native_intra, native.native_unbound
    ):
        raise ValueError("Native score, coordinates or saved inner energy changed")
    return dict(case, source_pool=str(pool), source_model=native.model, receptor=str(receptor),
                fixed_pose=str(fixed), reference_score=native.native_score,
                historical_search_wall_seconds=dock["wall_seconds"])


def run_score(case, binary, inputs, waters, sidechains, timeout):
    start = perf_counter()
    record = dict(case)
    try:
        record = freeze_pose(case, inputs)
        result = score_pool(binary, Path(record["receptor"]), Path(record["fixed_pose"]),
                            waters=waters, sidechains=sidechains, timeout=timeout)
        [pose] = result["poses"]
        record.update(status="ok", scores={arm: pose[arm] for arm in ARMS},
                      abstentions=result["abstentions"], score_payload=result)
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as error:
        record.update(status="failed", reason=str(error))
    record["wall_seconds"] = perf_counter() - start
    return record


def freeze(config_path, binary, output, reuse_path=None):
    root = config_path.resolve().parent
    config = json.loads(config_path.read_text())
    sources = configuration(root / config["rbfe_config"])
    cases = list(build_cases(sources))
    if len(cases) != config["expected_ligands"]:
        raise ValueError("The selected series has the wrong cardinality")
    waters, sidechains = ((root / config[k]).resolve() for k in ("waters", "sidechains"))
    water_evidence = waters if waters.is_file() else None
    waters = waters if water_evidence and waters.stat().st_size else None
    sidechains = sidechains if sidechains.is_file() else None
    binary = binary.resolve(strict=True)
    binary_id = digest(binary)
    code_id = digest(Path(__file__).with_name("score_physical_pool.py"))
    config_id = digest(config_path)
    water_id = digest(water_evidence) if water_evidence else None
    sidechain_id = digest(sidechains) if sidechains else None
    identity = dict(configuration_sha256=config_id, binary_sha256=binary_id,
                    score_physical_pool_sha256=code_id, waters_sha256=water_id, sidechains_sha256=sidechain_id)
    reused = {}
    if reuse_path:
        previous = json.loads(reuse_path.read_text())
        if any(previous.get(k) != v for k, v in identity.items()):
            raise ValueError("Reuse requires identical config, code, engine and evidence")
        reused = {r["name"]: r for r in previous["ligands"] if r["status"] == "ok"}
    output.mkdir()
    inputs = output / "fixed_poses"
    inputs.mkdir()
    start = perf_counter()
    pending = [case for case in cases if case["name"] not in reused]
    with ThreadPoolExecutor(max_workers=config["score_workers"]) as runner:
        scored = list(runner.map(lambda case: run_score(case, binary, inputs, waters, sidechains,
                                                       config["score_timeout_seconds"]), pending))
    complete = {**reused, **{row["name"]: row for row in scored}}
    rows = [complete[case["name"]] for case in cases]
    if digest(binary) != binary_id or digest(Path(__file__).with_name("score_physical_pool.py")) != code_id:
        raise RuntimeError("Scorer changed while the run was active")
    put_json(output / "predictions.json", dict(frozen_utc=now(), **identity,
        binary=str(binary), prepared_waters=str(water_evidence) if water_evidence else None,
        prepared_sidechains=str(sidechains) if sidechains else None,
        reuse_from=str(reuse_path) if reuse_path else None, reused_predictions=len(reused),
        wall_seconds=perf_counter()-start+(previous["wall_seconds"] if reuse_path else 0), ligands=rows))
    print([(row["name"], row["status"]) for row in rows], flush=True)


def compare(config_path, frozen_path, output):
    root = config_path.resolve().parent
    config = json.loads(config_path.read_text())
    sources = configuration(root / config["rbfe_config"])
    frozen = json.loads(frozen_path.read_text())
    if digest(config_path) != frozen["configuration_sha256"]:
        raise ValueError("Physical scoring configuration differs from freeze")
    measurements = yaml.safe_load(sources["measurements"].read_text())
    rows = []
    for prediction in frozen["ligands"]:
        name = prediction["name"]
        measured = float(measurements[name]["measurement"]["value"])*1000
        if not isfinite(measured) or measured <= 0:
            raise ValueError(f"{name} has no valid positive Ki")
        rows.append(dict(name=name, Ki_nM=measured, status=prediction["status"],
                         reference=prediction["reference"], wall_seconds=prediction["wall_seconds"],
                         **prediction.get("scores", {arm: None for arm in ARMS})))
    pairs = []
    for first, second in combinations(rows, 2):
        low, high = sorted((first, second), key=lambda row: row["Ki_nM"])
        ratio = high["Ki_nM"]/low["Ki_nM"]
        if ratio < config["minimum_ki_ratio"]:
            continue
        row = dict(lower_Ki=low["name"], higher_Ki=high["name"], ratio=ratio)
        for arm in ARMS:
            gap = None if low[arm] is None or high[arm] is None else low[arm] - high[arm]
            row[f"{arm}_gap"] = gap
            row[f"{arm}_correct"] = None if gap is None else 0.5 if gap == 0 else int(gap < 0)
        pairs.append(row)
    stats = {}
    for arm in ARMS:
        valid = [r for r in rows if r[arm] is not None]
        key = f"{arm}_correct"
        comparable = [p for p in pairs if p[key] is not None]
        rho = float(spearmanr([r["Ki_nM"] for r in valid], [r[arm] for r in valid]).statistic) if len(valid)>1 else None
        stats[arm] = dict(ligand_coverage=len(valid), pair_coverage=len(comparable),
                          correct=sum(p[key] for p in comparable), ties=sum(p[key]==0.5 for p in comparable),
                          gains_vs_native=sum(p["native_correct"] is not None and p[key] > p["native_correct"] for p in comparable),
                          losses_vs_native=sum(p["native_correct"] is not None and p[key] < p["native_correct"] for p in comparable),
                          spearman=rho if rho is not None and isfinite(rho) else None)
    by_name = {row["name"]: row for row in rows}
    edges = json.loads((root/config["rbfe_selection"]).read_text())["selected"]
    rbfe_rows = []
    for edge in edges:
        a, b = [by_name["lig_" + edge[x]] for x in ("a", "b")]
        rbfe_rows.append(dict(a=a["name"], b=b["name"], A_Ki_nM=a["Ki_nM"], B_Ki_nM=b["Ki_nM"],
                              **{f"{k}_B_minus_A": b[k]-a[k] if a[k] is not None and b[k] is not None else None
                                 for k in ARMS}))
    output.mkdir()
    put_csv(output/"by_ligand.csv", rows)
    put_csv(output/"by_pair.csv", pairs)
    put_csv(output/"rbfe_edges.csv", rbfe_rows)
    put_json(output/"summary.json", dict(frozen_utc=frozen["frozen_utc"], analysed_utc=now(),
        failures=sum(row["status"] != "ok" for row in rows), eligible_pairs=len(pairs),
        prediction_wall_seconds=frozen["wall_seconds"], engine_call_seconds=sum(
            row["score_payload"]["poses"][0]["engine"]["wall_seconds"]
            + sum(alt["engine"]["wall_seconds"] for alt in row["score_payload"]["poses"][0]["alternate_states"])
            for row in frozen["ligands"] if row["status"]=="ok"),
        score_process_wall_seconds=sum(row["wall_seconds"] for row in rows), arms=stats))
    print(stats, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="stage", required=True)
    for name in ("freeze", "compare"):
        action = actions.add_parser(name)
        action.add_argument("--config", type=Path, required=True)
        action.add_argument("--output", type=Path, required=True)
    actions.choices["freeze"].add_argument("--binary", type=Path, required=True)
    actions.choices["freeze"].add_argument("--reuse-predictions", type=Path)
    actions.choices["compare"].add_argument("--frozen", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "freeze":
        freeze(args.config, args.binary, args.output, args.reuse_predictions)
    else:
        compare(args.config, args.frozen, args.output)


if __name__ == "__main__":
    main()
