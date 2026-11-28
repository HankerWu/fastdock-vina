"""Reconstruct native merge decisions and audit them against real engine membership."""
import argparse
import json
import lzma
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stage_local_hunt"))
from coverage import specs, native, run_panel

PREVIOUS = HERE.parent / "stage_local_hunt/coverage"


def read_cells():
    return {row["key"]: row for row in json.loads((PREVIOUS / "costs.json").read_text())}


def live(tmp, work, builds):
    work.mkdir(parents=True, exist_ok=False)
    baselines = read_cells()
    costs = []
    for _, job in specs(tmp):
        _, inputs, q, r, _, seed, _ = job
        for arm in ("stable", "candidate"):
            name = f"reciprocal-{q['pdb']}-{r['pdb']}-{seed}-{arm}"
            trace = tmp / (name + ".tsv")
            binary = builds / arm / "build/vina"
            os.environ["FD_MERGE_TRACE"] = str(trace)
            commands = []
            execute = run_panel.execute

            def capture(command, *args):
                commands.append(command.copy())
                return execute(command, *args)

            run_panel.execute = capture
            try:
                row, _ = run_panel.dock((arm, binary, work, job))
            finally:
                run_panel.execute = execute
            assert len(commands) == 1
            out = work / (name + ".tsv.xz")
            if trace.exists():
                out.write_bytes(lzma.compress(trace.read_bytes()))
            row.update(binary_sha256=run_panel.digest(binary), command=commands[0],
                       ligand_sha256=run_panel.digest(inputs / q["pdb"] / "ligand.pdbqt"),
                       receptor_sha256=run_panel.digest(inputs / r["pdb"] / "receptor.pdbqt"),
                       trace_sha256=run_panel.digest(out) if out.exists() else None)
            costs.append(row)
            (work / "costs.json").write_text(json.dumps(costs, indent=2) + "\n")
            assert row["status"] == "ok"
            assert all(row[k] == baselines[name][k] for k in
                       ("pool_sha256", "hunt_evals", "second_evals", "explicit_evals", "explicit_calls"))
    assert len(costs) == 12


def load_trace(path):
    poses, members = [], []
    for line in lzma.decompress(path.read_bytes()).decode().splitlines():
        event, *args = line.split("\t")
        if event == "P":
            radius, cap = float(args[0]), int(args[1])
        elif event == "C":
            t, pos, e, state, *xyz = args
            poses.append(dict(t=int(t), pos=int(pos), e=float(e), state=state,
                              xyz=np.array(xyz, dtype=float).reshape(-1, 3)))
        elif event == "M":
            members.append(list(map(int, args)))
        elif event == "I":
            final = list(map(int, args))
        elif event == "G":
            sorted_final = list(map(int, args))
        else:
            raise ValueError(event)
    assert len(poses) == len(members)
    return radius, cap, poses, members, final, sorted_final


def attach_crystal_labels(tmp, baseline, key, job, poses, final):
    _, inputs, query, _, ref, _, _ = job
    serials = [int(l[6:11]) for l in (inputs / query["pdb"] / "ligand.pdbqt").read_text().splitlines()
               if l.startswith(("ATOM  ", "HETATM")) and l[77:].strip() not in {"HS", "HD", "H"}]
    order = [serials.index(i) for i in query["serial"]]
    old = PREVIOUS / (key + ".xz")
    assert run_panel.digest(old) == baseline["compressed_sha256"]
    archive = tmp / "old.trace"
    archive.write_bytes(lzma.decompress(old.read_bytes()))
    paths, explicit = native(archive)
    assert [poses[j]["state"] for j in final] == [f[3] for f in explicit["I"]]
    for t, path in enumerate(paths):
        retained = next(event["M"] for event in reversed(path) if event["M"])
        saved = {d[3]: d for step in path for d in step["D"]}
        group = [p for p in poses if p["t"] == t]
        assert [p["state"] for p in group] == [r[3] for r in retained]
        for i, p in enumerate(group):
            row = saved[p["state"]]
            assert p["pos"] == i and p["e"] == float(row[2])
            # The old crystal-labelled observer used PDBQT record order, which
            # is distinct from the order in the native overlap vector.
            cart = np.fromstring(row[8], dtype=int, sep=" ").reshape(-1, 3) / 1000
            err = np.sqrt(((cart[order][None] - ref) ** 2).sum(axis=-1).mean(axis=1)).min()
            p.update(step=int(row[1]), crystal=float(err))


def replay_cell(key, poses, snapshots, final, global_expected, radius, cap):
    energies = [p["e"] for p in poses]
    coords = np.array([p["xyz"] for p in poses])
    distances = np.sqrt(((coords[:, None] - coords) ** 2).sum(axis=-1).mean(axis=-1))
    near = {i for i,p in enumerate(poses) if p["crystal"] <= 2.0}
    comparisons, arrivals, retained, events = [], [], [], []
    for i, pose in enumerate(poses):
        before = retained.copy()
        e = energies[i]
        effects = []
        blocker = None

        def record(j, check, condition):
            comparisons.append(dict(key=key, arrival=i, resident=j, comparison=check,
                                    candidate_e=e, resident_e=energies[j],
                                    index_rmsd_A=distances[i,j] if "overlap" in check else "",
                                    condition_true=int(condition)))
            return condition

        if len(retained) == cap and record(retained[-1], "capacity_energy_ge", e >= energies[retained[-1]]):
            action, blocker = "capacity_admission", retained[-1]
        else:
            pos = 0
            while pos < len(retained):
                j = retained[pos]
                if record(j, "insert_before_lower_priority", e < energies[j]):
                    break
                if record(j, "better_or_tied_overlap", distances[i,j] < radius):
                    blocker = j
                    break
                pos += 1
            if blocker is None:
                action = "insert"
                tail = []
                for j in retained[pos:]:
                    if record(j, "worse_overlap", distances[i,j] < radius):
                        effects.append((j, "overlap_eviction"))
                    else:
                        tail.append(j)
                retained = retained[:pos] + [i] + tail
                if len(retained) > cap:
                    effects.append((retained.pop(), "capacity_eviction"))
            else:
                action = "overlap_admission"
        if blocker is not None:
            events.append((i, blocker, i, action))
        for victim, reason in effects:
            events.append((i, i, victim, reason))
        assert retained == snapshots[i], (key, i, retained, snapshots[i])
        arrivals.append(dict(key=key, arrival=i, trajectory=pose["t"]+1, local_rank=pose["pos"]+1,
                             save_step=pose["step"], native_e=e, crystal_A=pose["crystal"],
                             before=" ".join(map(str,before)), action=action,
                             blocker="" if blocker is None else blocker,
                             evictions=" ".join(f"{k}/{reason}" for k,reason in effects),
                             after=" ".join(map(str,retained))))
    assert retained == final
    # Offline stable, ascending energy across all trajectories, not re-refined.
    global_ids = []
    for i in sorted(range(len(poses)), key=lambda j: energies[j]):
        if all(distances[i,j] >= radius for j in global_ids):
            global_ids.append(i)
        if len(global_ids) == cap:
            break
    assert global_ids == global_expected
    losses = []
    for when, rival, victim, reason in events:
        losses.append(dict(key=key, arrival=when, lost=victim, rival=rival, reason=reason,
                           lost_e=energies[victim], rival_e=energies[rival],
                           diagnostic_pair_distance_A=distances[victim,rival],
                           lost_crystal_A=poses[victim]["crystal"], rival_crystal_A=poses[rival]["crystal"],
                           rival_is_final=int(rival in final),
                           rival_later_lost_at=next((k for k,_,v,_ in events if k>when and v==rival),"")))
    def index_list(ids):
        return " ".join(map(str,ids))
    best_id = min(final, key=lambda j: poses[j]["crystal"])
    assert best_id in snapshots[-1] and all(row["lost"] != best_id for row in losses)
    assert set(range(len(poses))) - set(final) == {row["lost"] for row in losses}
    assert len(losses) + len(final) == len(poses)
    summary = dict(key=key, inputs=len(poses), radius_A=radius, capacity=cap,
                   native=index_list(final), global_order=index_list(global_ids),
                   native_oracle_A=poses[best_id]["crystal"],
                   global_order_oracle_A=min(poses[j]["crystal"] for j in global_ids),
                   terminal_near=index_list(sorted(near)), native_near=index_list(sorted(near&set(final))),
                   global_order_near=index_list(sorted(near&set(global_ids))),
                   rejected=len(losses), comparisons=len(comparisons))
    return summary, comparisons, arrivals, losses


def analyze(tmp, raw, out):
    out.mkdir(parents=True, exist_ok=False)
    costs = {c["key"]: c for c in json.loads((raw / "costs.json").read_text())}
    baselines = read_cells()
    tables = {s: [] for s in ("cells","comparisons","arrivals","losses")}
    for _, job in specs(tmp):
        _, _, q, r, _, seed, _ = job
        for arm in ("stable", "candidate"):
            key = f"reciprocal-{q['pdb']}-{r['pdb']}-{seed}-{arm}"
            packed = raw / (key + ".tsv.xz")
            assert run_panel.digest(packed) == costs[key]["trace_sha256"]
            assert run_panel.digest(raw / (key + ".pdbqt")) == baselines[key]["pool_sha256"] == costs[key]["pool_sha256"]
            radius, cap, poses, snapshots, final, global_final = load_trace(packed)
            attach_crystal_labels(tmp, baselines[key], key, job, poses, final)
            summary, comparisons, arrivals, losses = replay_cell(key, poses, snapshots, final,
                                                                 global_final, radius, cap)
            tables["cells"].append(summary)
            tables["comparisons"].extend(comparisons)
            tables["arrivals"].extend(arrivals)
            tables["losses"].extend(losses)
    for name, rows in tables.items():
        run_panel.csvfile(out / (name + ".csv"), rows)
    assert len(tables["cells"]) == 12


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run","analyze"))
    parser.add_argument("input", type=Path, help="observer build directory for run; raw directory for analyze")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        if args.mode == "run":
            live(Path(tmp), args.output.resolve(), args.input.resolve())
        else:
            analyze(Path(tmp), args.input.resolve(), args.output.resolve())
