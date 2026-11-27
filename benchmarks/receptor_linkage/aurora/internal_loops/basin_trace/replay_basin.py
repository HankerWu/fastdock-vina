"""Run label-blind searches, then measure coordinates using independently supplied crystal data."""
import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import replay
from summarize import save

TAGS = ["start", "gradient", "adopted", "hunt_return", "mc_accept", "second_return",
        "saved", "pooled", "explicit_input", "explicit_output"]
PHASES = ["other", "hunt", "second", "explicit"]
SEEDS = (61231, 45101, 45102, 45103, 45104)


def read_trace(path):
    with path.open("rb") as stream:
        n = int.from_bytes(stream.read(4), "little")
        if n < 1 or n > 512:
            raise ValueError("Expected a little-endian geometry stream")
        serials = np.fromfile(stream, "<u4", n)
    dtype = np.dtype([(k, "<u4") for k in ("traj", "step", "phase", "event", "local")] +
                     [("energy", "<f8"), ("xyz", "<f8", (n, 3))])
    return serials.tolist(), np.memmap(path, dtype=dtype, mode="r", offset=4*(n+1))


def stage_labels(tr):
    names = np.asarray(["unused", *TAGS])
    labels = names[tr["event"]].astype("U24")
    for p, phase in enumerate(PHASES[1:], 1):
        mask = tr["phase"] == p
        gradient = tr["event"] == 2
        labels[mask & gradient & (tr["local"] == 1)] = phase + "_start"
        labels[mask & gradient & (tr["local"] > 1)] = phase + "_trial_rejected"
        labels[mask & (tr["event"] == 3)] = phase + "_trial_adopted"
    return labels


def distances(tr, refs):
    rmsd = np.full(len(tr), np.inf)
    for i in range(0, len(tr), 1024):
        stop = min(len(tr), i + 1024)
        for reference in refs:
            delta = tr["xyz"][i:stop] - reference
            rmsd[i:stop] = np.minimum(rmsd[i:stop], np.sqrt((delta**2).sum(axis=2).mean(axis=1)))
    if not np.all(np.isfinite(rmsd)):
        raise ValueError("Nonfinite ligand coordinates")
    return rmsd


COUNTERS = ("Search trajectory:", "Hunt chain:", "Pocket start:", "Search counters:")


def counters(path):
    return [line for line in path.read_text().splitlines() if line.startswith(COUNTERS)]


def summarize(path, work, drug, seed, graph, dest, production, log):
    source_serials, tr = read_trace(path)
    _, _, target_serials, refs = replay.parameters(work, drug)
    if sorted(source_serials) != sorted(target_serials):
        raise ValueError("Observer and PDBQT serials do not match")
    in_trace_order = [target_serials.index(s) for s in source_serials]
    refs = refs[:, in_trace_order]
    distances_A = distances(tr, refs)
    labels = stage_labels(tr)
    metadata = {"drug": drug, "graph": graph, "seed": seed}
    case = f"{drug}_{graph}_{seed}"
    # Every evaluated coordinate has exactly one classification and agrees with CLI counters.
    for phase in PHASES[1:]:
        expected = sum(map(int, re.findall(rf"\b{phase}_evals=(\d+)", log.read_text())))
        measured = np.count_nonzero((tr["phase"] == PHASES.index(phase)) & np.isin(tr["event"], [2, 3]))
        if measured != expected or not expected:
            raise ValueError(("Incomplete objective calls", phase, measured, expected))
    rows = []
    witnesses = set()
    for stage in sorted(set(labels)):
        mask = labels == stage
        for trajectory in [*sorted(set(tr["traj"][mask])), "all"]:
            hits = np.flatnonzero(mask if trajectory == "all" else mask & (tr["traj"] == trajectory))
            best = int(hits[np.argmin(distances_A[hits])])
            witnesses.add(best)
            row = dict(**metadata, trajectory=trajectory, stage=stage, count=len(hits),
                       min_A=float(distances_A[best]), min_row=best,
                       min_tr=int(tr["traj"][best]), min_step=int(tr["step"][best]),
                       min_local=int(tr["local"][best]))
            for threshold in [2, 3]:
                near = hits[distances_A[hits] <= threshold]
                row[f"le{threshold}"] = len(near)
                row[f"first_le{threshold}_row"] = int(near[0]) if len(near) else ""
                if len(near):
                    witnesses.add(int(near[0]))
            rows.append(row)
    # Output serialization rounds each Cartesian component to 0.001 Å. Each published
    # result must correspond to one complete, unrounded final observer geometry.
    finals = tr["xyz"][labels == "explicit_output"]
    published = np.array(replay.geom.mode_coords(production, source_serials))
    deltas = np.max(np.abs(published[:, None] - finals), axis=(2, 3)).min(axis=1)
    if max(deltas) > 0.000501:
        raise ValueError(("Final coordinates not reproduced", deltas))
    # Keep the two MC steps around first local-return entry so that both entry
    # and its parent can be independently examined, along with the last actual
    # accepted parent if the immediately preceding MC step was rejected.
    for threshold in (2, 3):
        hits = np.flatnonzero((labels == "hunt_return") & (distances_A <= threshold))
        if len(hits):
            first = hits[0]
            same = tr["traj"] == tr["traj"][first]
            parents = np.flatnonzero(same & (tr["step"] < tr["step"][first]) & (labels == "mc_accept"))
            steps = [tr["step"][first], tr["step"][first] - 1]
            if len(parents):
                steps.append(tr["step"][parents[-1]])
            witnesses.update(np.flatnonzero(same & np.isin(tr["step"], steps)))
    chosen = np.asarray(sorted(witnesses))
    # One compact exact partition per run, including every gradient attempt.
    # For full raw coordinates retain the binary stream or replay its recorded digest.
    payload = {k: tr[k] for k in ("traj", "step", "phase", "event", "local")}
    np.savez_compressed(dest / (case + ".npz"), **payload, rmsd_A=distances_A,
                        witness_index=chosen, witness_energy=tr["energy"][chosen],
                        witness_xyz=tr["xyz"][chosen], serials=source_serials)
    return rows, float(max(deltas)), len(tr)


def run(work, past, binary, out, scratch):
    out.mkdir(exist_ok=True, parents=True)
    scratch.mkdir(exist_ok=True, parents=True)
    center = json.loads((work / "prepared/2X6E.json").read_text())["center"]
    box = [arg for dim, x in zip("xyz", center, strict=True)
           for arg in (f"--center_{dim}", x, f"--size_{dim}", 26)]
    receptor = work / replay.PROTEINS["1OL7_released"]
    rows, runs = [], []
    for drug in ("at9283", "danusertib"):
        for seed in SEEDS:
            for graph in (("complete", "inferred") if drug == "at9283" and seed == 61231 else ("complete",)):
                cohort = "historical" if seed == 61231 else "fresh"
                arm = "both_graph" if graph == "complete" else "both"
                previous = past/cohort/"blind"/f"{arm}_{drug}_{seed}.pdbqt"
                case = f"{drug}_{graph}_{seed}"
                raw = scratch / (case + ".bin")
                pose = scratch / (case + ".pdbqt")
                log = scratch / (case + ".log")
                ligand = work / "prepared" / (drug + ".pdbqt")
                bonds = work / "raw" / ("at9283_chemical.bonds" if drug == "at9283" else "danusertib_native.bonds")
                args = [binary, "--receptor", receptor, "--ligand", ligand, *box,
                        "--seed", seed, "--exhaustiveness", 2, "--cpu", 1, "--num_modes", 9,
                        "--energy_range", 100, "--max_evals", 25000, "--out", pose]
                if graph == "complete":
                    args += ["--ligand_bonds", bonds]
                args += ["--pose_selection", "late", "--search_counters"]
                # Diagnostic references are first loaded inside summarize, after process exit.
                os.environ["FD_BASIN_TRACE"] = str(raw)
                cost = replay.invoke(args, log)
                if pose.read_bytes() != previous.read_bytes() or counters(log) != counters(previous.with_suffix(".log")):
                    raise ValueError(("Archived production control mismatch", case))
                summary, final_delta, events = summarize(raw, work, drug, seed, graph, out, pose, log)
                rows.extend(summary)
                runs.append(dict(drug=drug, graph=graph, seed=seed, events=events, **cost,
                                 evals=sum(map(int, re.findall(r"\b(?:hunt|second|explicit)_evals=(\d+)", log.read_text()))),
                                 pool_sha256=replay.digest(pose), raw_sha256=replay.digest(raw),
                                 receptor_sha256=replay.digest(receptor), ligand_sha256=replay.digest(ligand),
                                 bonds_sha256=replay.digest(bonds) if graph == "complete" else "",
                                 max_coordinate_rounding_A=final_delta))
                save(out/"stages.csv", rows)
                save(out/"runs.csv", runs)
                print(case, events, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("work", "segment_run", "binary", "results", "scratch"):
        parser.add_argument(key, type=Path)
    run(*(p.resolve() for p in vars(parser.parse_args()).values()))
