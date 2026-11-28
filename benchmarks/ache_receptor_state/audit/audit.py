"""Verify receptor-substitution results from raw records, including atom identity and energy receipts."""

import argparse
import hashlib
import json
import lzma
import re
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
BENCH = STUDY.parent
PRIOR = BENCH / "water_rank_envelope/transfer"
sys.path.insert(0, str(BENCH / "cubic_grid/human_panel"))
from verify import reference  # independent deposited-PDB/CCD graph builder

COMPONENTS = {"4EY5": "HUP", "4EY6": "GNT", "4EY7": "E20"}
ARMS = ("pro24k", "pro48k", "official24k")
COUNTERS = (
    "hunt_calls",
    "hunt_evals",
    "second_calls",
    "second_evals",
    "explicit_calls",
    "explicit_evals",
)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def xyz(text):
    return {
        int(l[6:11]): [float(l[i : i + 8]) for i in (30, 38, 46)]
        for l in text.splitlines()
        if l.startswith(("ATOM  ", "HETATM"))
    }


def transform(source, target):
    src, sc, _ = source
    dst, tc, _ = target
    keys = sorted(
        k
        for k in src.keys() & dst.keys()
        if np.linalg.norm(src[k] - sc) < 14 and np.linalg.norm(dst[k] - tc) < 14
    )
    a, b = [np.array([z[k] for k in keys]) for z in (src, dst)]
    rot, _ = Rotation.align_vectors(b - b.mean(axis=0), a - a.mean(axis=0))

    def apply(pos):
        return rot.apply(pos - a.mean(axis=0)) + b.mean(axis=0)

    return apply, {"shared_ca": len(keys), "rmsd_A": rmsd(apply(a), b)}


def all_coordinates(text):
    sites = xyz(text)
    return np.array([sites[k] for k in sorted(sites)])


def rmsd(a, b):
    return float(np.linalg.norm(a - b) / np.sqrt(len(a)))


class Frame:
    def __init__(self, work):
        raw = work / "raw"
        raw.mkdir()
        (raw / "4EY5.pdb").write_bytes(lzma.decompress((PRIOR / "4EY5.pdb.xz").read_bytes()))
        (raw / "HUP.cif").write_bytes((PRIOR / "HUP.cif").read_bytes())
        with tarfile.open(BENCH / "receptor_amide_states/transfer/transfer_evidence.tar.xz") as t:
            for fname in ("4EY6.pdb", "4EY7.pdb", "GNT.cif", "E20.cif"):
                (raw / fname).write_bytes(t.extractfile("novel/raw/" + fname).read())
        self.struct = {r: reference(work, r, c) for r, c in COMPONENTS.items()}
        self.moves = {}
        self.fits = {}
        for q in COMPONENTS:
            for r in COMPONENTS:
                self.moves[q, r], self.fits[q + "_" + r] = transform(self.struct[q], self.struct[r])
        ligand = (PRIOR / "HUP.pdbqt").read_text()
        (smiles,) = re.findall(r"^REMARK SMILES (\S+)$", ligand, re.MULTILINE)
        self.graph = Chem.MolFromSmiles(smiles)
        self.identity = Chem.MolToSmiles(self.graph)
        tokens = [
            int(i)
            for l in ligand.splitlines()
            if l.startswith("REMARK SMILES IDX")
            for i in l.split()[3:]
        ]
        self.order = [v for k, v in sorted(zip(tokens[::2], tokens[1::2], strict=True))]
        truth = self.struct["4EY5"][2]
        mappings = truth.GetSubstructMatches(self.graph, uniquify=False, useChirality=True)
        assert len(mappings) == 1 and len(self.order) == 18
        self.crystal = truth.GetConformer().GetPositions()[list(mappings[0])]
        self.targets = {r: self.moves["4EY5", r](self.crystal) for r in COMPONENTS}
        self.template = ligand
        self.input = self.heavy(ligand)
        excluded = []
        for line in (raw / "4EY5.pdb").read_text().splitlines():
            if not line.startswith(("HETATM", "ATOM  ")) or line[76:78].strip() in ("H", "D"):
                continue
            if (
                line[17:20] in ("HOH", "DOD")
                or line[21] == "A"
                and (line.startswith("ATOM  ") or line[17:20] == "HUP")
            ):
                continue
            pos = np.array([float(line[i : i + 8]) for i in (30, 38, 46)])
            excluded.append(
                (float(np.linalg.norm(self.crystal - pos, axis=1).min()), line[12:27].strip())
            )
        self.cognate_guard = min(excluded)

    def heavy(self, model):
        coords = xyz(model)
        x = np.array([coords[i] for i in self.order])
        conf = Chem.Conformer(len(x))
        for i, position in enumerate(x):
            conf.SetAtomPosition(i, position)
        mol = Chem.Mol(self.graph)
        mol.AddConformer(conf)
        Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
        assert Chem.MolToSmiles(mol) == self.identity
        return x


def compare(new, path, key, columns, tol=1e-8):
    """Match identities strictly before comparing numbers."""
    other = pd.read_csv(path)
    l = new.set_index(key).sort_index()
    r = other.set_index(key).sort_index()
    if len(l) != len(r) or not l.index.equals(r.index) or not l.index.is_unique:
        raise ValueError(f"Independent keys disagree with {path}")
    error = {}
    for col in columns:
        d = np.abs(l[col].to_numpy(dtype=float) - r[col].to_numpy(dtype=float))
        if not np.all(d <= tol):
            raise ValueError(f"{path}:{col} differs up to {d.max()}")
        error[col] = float(d.max())
    return error


def counters(text):
    ans = {}
    for n in COUNTERS:
        vals = re.findall(r"\b" + n + r"=(\d+)", text)
        assert len(vals) == (1 if n.startswith("explicit") else 3)
        ans[n] = sum(map(int, vals))
    return ans


def distribution(v):
    return dict(
        zip(
            ("mean", "q1", "median", "q3", "p90", "p95"),
            [float(np.mean(v)), *np.quantile(v, [0.25, 0.5, 0.75, 0.90, 0.95]).tolist()],
            strict=True,
        )
    )


def blind(frame, out, path, compare_primary):
    modes, cells, production, replays = [], [], {}, {}
    previous = {
        arm: (PRIOR / ("" if arm == "pro24k" else arm) / "dry.pdbqt").read_bytes() for arm in ARMS
    }
    key = ["receiver", "arm", "seed"]
    issued_centers = []
    for p in sorted((path / "raw").glob("*/*/dry.pdbqt")):
        arm, seed = p.parent.name.rsplit("_", 1)
        seed, rec = int(seed), p.parents[1].name
        k = {"receiver": rec, "arm": arm, "seed": seed}
        cmd = json.loads((p.parent / "command.json").read_text())
        site = (
            frame.struct[rec][1]
            if arm != "pro24k_mapped"
            else frame.moves["4EY6", "4EY7"](frame.struct["4EY6"][1])
        )
        center = [float(cmd[cmd.index("--center_" + a) + 1]) for a in "xyz"]
        issued_centers.append(float(np.max(np.abs(center - site))))
        assert all(float(cmd[cmd.index("--size_" + a) + 1]) == 24 for a in "xyz")
        assert int(cmd[cmd.index("--seed") + 1]) == seed and max(issued_centers) < 1e-10
        assert not (arm.endswith("mapped") and rec != "4EY7")
        if rec == "4EY6" and seed == 1103501:
            assert p.read_bytes() == previous[arm]
        local = []
        for rank, mode in enumerate(p.read_text().split("ENDMDL")[:-1], 1):
            score = float(re.search(r"REMARK VINA RESULT:\s+(\S+)", mode)[1])
            distance = rmsd(frame.heavy(mode), frame.targets[rec])
            local.append(distance)
            modes.append(dict(**k, mode=rank, rmsd_A=distance, score=score))
        assert 0 < len(local) <= 9
        oracle = min(local)
        vals = dict(
            **k,
            modes=len(local),
            top1_A=local[0],
            best5_A=min(local[:5]),
            oracle_A=oracle,
            oracle_mode=local.index(oracle) + 1,
            regret_A=local[0] - oracle,
            serious=local[0] > 5,
        )
        # Compare counted replay to release files before using observer counts.
        for run in (path / "replay" / rec / p.parent.name).iterdir():
            assert (run / "dry.pdbqt").read_bytes() == p.read_bytes()
            replays[rec, arm, seed, run.name] = sha((run / "dry.pdbqt").read_bytes())
        counter_log = (
            path / "replay" / rec / p.parent.name / "counter/dry.log"
            if arm == "official24k"
            else p.with_suffix(".log")
        )
        vals.update(counters(counter_log.read_text()))
        production[rec, arm, seed] = sha(p.read_bytes())
        cells.append(vals)
    assert len(cells) == 21 and len(replays) == 27
    modes, cells = pd.DataFrame(modes), pd.DataFrame(cells)
    cells.to_csv(out / "blind_cells.csv", index=False)
    group = []
    for (rec, arm), g in cells.groupby(["receiver", "arm"]):
        group.append(
            dict(
                receiver=rec,
                arm=arm,
                cells=len(g),
                severe=int(g.serious.sum()),
                **distribution(g.top1_A),
            )
        )
    pd.DataFrame(group).to_csv(out / "blind_groups.csv", index=False)
    pairs = []
    for arm in ARMS:
        for seed in (1103501, 1103503, 1103507):
            subset = cells[(cells.arm == arm) & (cells.seed == seed)].set_index("receiver")
            for col in ("top1_A", "oracle_A", "regret_A"):
                pairs.append(
                    {
                        "arm": arm,
                        "seed": seed,
                        "metric": col,
                        "delta_4EY7_minus_4EY6": float(
                            subset.loc["4EY7", col] - subset.loc["4EY6", col]
                        ),
                    }
                )
    pd.DataFrame(pairs).to_csv(out / "blind_receiver_differences.csv", index=False)
    audit = {
        "production_pools": len(cells),
        "modes": len(modes),
        "replay_byte_matches": len(replays),
        "max_command_center_error_A": max(issued_centers),
        "original_pool_byte_matches": 3,
        "zero_pools": int((cells.modes == 0).sum()),
        "severe": int(cells.serious.sum()),
        "min_published_rmsd_A": cells.groupby("receiver").oracle_A.min().to_dict(),
        "counters_from_matched_outputs": {c: int(cells[c].sum()) for c in COUNTERS},
    }
    original = pd.read_csv(path / "runs.csv")
    for row in original.itertuples():
        assert row.output_sha256 == production[row.receiver, row.arm, row.seed]
    observed = pd.read_csv(path / "replay/replays.csv")
    for row in observed.itertuples():
        assert row.output_sha256 == replays[row.receiver, row.arm, row.seed, row.binary]
    if compare_primary:
        audit["mode_max_abs_error"] = compare(
            modes, path / "modes.csv", key + ["mode"], ["rmsd_A", "score"]
        )
        audit["cell_max_abs_error"] = compare(
            cells,
            path / "cells.csv",
            key,
            ["top1_A", "best5_A", "oracle_A", "oracle_mode", "regret_A", "modes"],
        )
        contrasts = pd.read_csv(path / "paired.csv")
        indexed = cells.set_index(key)
        pair_errors = []
        for r in contrasts.itertuples():
            a, b = indexed.loc[(r.receiver1, r.arm1, r.seed)], indexed.loc[(r.receiver0, r.arm0, r.seed)]
            for col in ("top1_A", "best5_A", "oracle_A", "regret_A"):
                pair_errors.append(abs(float(a[col] - b[col]) - getattr(r, col)))
        if not pair_errors or max(pair_errors) > 1e-8:
            raise ValueError("Paired contrasts disagree")
        audit["paired_rows_checked"] = len(contrasts)
        audit["max_paired_metric_error_A"] = max(pair_errors)
        audit["max_group_stat_error_A"] = 0.0
        for r in pd.read_csv(path / "groups.csv").itertuples():
            g = cells[(cells.receiver == r.receiver) & (cells.arm == r.arm)]
            for stat, val in distribution(g[r.measure]).items():
                audit["max_group_stat_error_A"] = max(
                    audit["max_group_stat_error_A"], abs(val - getattr(r, stat))
                )
        assert audit["max_group_stat_error_A"] < 1e-8
        ref = pd.concat(
            (original[original.arm != "official24k"], observed[observed.binary == "counter"])
        )
        for row in cells.itertuples():
            recorded = ref[
                (ref.arm == row.arm) & (ref.receiver == row.receiver) & (ref.seed == row.seed)
            ].iloc[0]
            assert all(recorded[n] == getattr(row, n) for n in COUNTERS)
    return audit


def score_log(text):
    component = lambda name: float(re.search(name + r"\s*:\s*(\S+)", text)[1])
    score = component("Estimated Free Energy of Binding")
    parts = [
        component(x)
        for x in (
            "Final Intermolecular Energy",
            "Final Total Internal Energy",
            "Torsional Free Energy",
            "Unbound System's Energy",
        )
    ]
    assert parts[1:] == [0.0, 0.0, 0.0] and abs(parts[0] - score) < 1e-10
    return score


def diagnostic(frame, out, directory):
    rows = []
    pools = {
        arm: (PRIOR / ("" if arm == "pro24k" else arm) / "dry.pdbqt").read_text().split("ENDMDL")[:-1]
        for arm in ARMS
    }
    originals = {
        f"{arm}_{i:02d}": model for arm, models in pools.items() for i, model in enumerate(models, 1)
    }
    oracle = {
        arm: 1 + int(np.argmin([rmsd(frame.heavy(m), frame.targets["4EY6"]) for m in pool]))
        for arm, pool in pools.items()
    }
    assert oracle == {"pro24k": 1, "pro48k": 2, "official24k": 4}
    rot, _ = Rotation.align_vectors(
        frame.crystal - frame.crystal.mean(0), frame.input - frame.input.mean(0)
    )
    projected = rot.apply(all_coordinates(frame.template) - frame.input.mean(0)) + frame.crystal.mean(0)
    rigid_error = []
    path = directory / "raw.tar.xz"
    assert frame.cognate_guard[0] < 8 and "DMS A 605" in frame.cognate_guard[1]
    with tarfile.open(path) as archive:
        assert not any(n.startswith("4EY5/") and n.endswith("/input.pdbqt") for n in archive.getnames())
        for rec in ("4EY6", "4EY7"):
            for name in [*originals, "projection"]:

                def raw(file, directory=f"{rec}/{name}"):
                    return archive.extractfile(f"{directory}/{file}").read().decode()

                start = frame.heavy(raw("input.pdbqt"))
                end = frame.heavy(raw("output.pdbqt"))
                assert raw("observed.pdbqt") == raw("output.pdbqt")
                expected = (
                    frame.moves["4EY5", rec](projected)
                    if name == "projection"
                    else frame.moves["4EY6", rec](all_coordinates(originals[name]))
                )
                # A final PDBQT writer rounds each coordinate to .001 Å.
                err = float(np.max(np.abs(expected - all_coordinates(raw("input.pdbqt")))))
                assert err < 0.000501, (rec, name, err)
                rigid_error.append(err)
                obs = re.findall(r"REFINE_OBSERVE\s+(\S+)\s+(\S+)\s+(\d+)", raw("observer.log"))
                assert obs
                rows.append(
                    {
                        "receiver": rec,
                        "name": name,
                        "initial": score_log(raw("initial.log")),
                        "final": score_log(raw("final.log")),
                        "start_A": rmsd(start, frame.targets[rec]),
                        "end_A": rmsd(end, frame.targets[rec]),
                        "calls": sum(int(v[2]) for v in obs),
                    }
                )
    df = pd.DataFrame(rows)
    df.to_csv(out / "diagnostic_endpoints.csv", index=False)
    gaps, projection = [], []
    for rec, sub in df.groupby("receiver"):
        keyed = sub.set_index("name")
        projected = keyed.loc["projection"]
        old = keyed.drop("projection")
        best = old.final.min()
        projection.append(
            {
                "receiver": rec,
                "rmsd_A": projected.end_A,
                "score": projected.final,
                "frozen_endpoint_best": float(best),
                "delta_score": float(projected.final - best),
                "retained": projected.end_A <= 2,
                "wins": projected.final < best - 0.01,
            }
        )
        for arm, idx in oracle.items():
            rank1 = keyed.loc[f"{arm}_01"]
            near = keyed.loc[f"{arm}_{idx:02d}"]
            gaps.append(
                {
                    "receiver": rec,
                    "arm": arm,
                    "oracle_mode": idx,
                    "start_gap": near.initial - rank1.initial,
                    "final_gap": near.final - rank1.final,
                }
            )
    pd.DataFrame(gaps).to_csv(out / "diagnostic_gaps.csv", index=False)
    pd.DataFrame(projection).to_csv(out / "diagnostic_projection.csv", index=False)
    primary = df.rename(
        columns={
            "name": "start",
            "initial": "score_in",
            "final": "score_out",
            "start_A": "rmsd_in",
            "end_A": "rmsd_out",
            "calls": "gradient_evals",
        }
    )
    primary_error = compare(
        primary,
        directory / "endpoints.csv",
        ["receiver", "start"],
        ["rmsd_in", "rmsd_out", "score_in", "score_out", "gradient_evals"],
    )
    same_state_violations = 0
    for arm, pool in pools.items():
        printed = [float(re.search(r"REMARK VINA RESULT:\s+(\S+)", m)[1]) for m in pool]
        same_state = df[(df.receiver == "4EY6") & df.name.str.startswith(arm + "_")].sort_values("name")
        vals = same_state.initial.to_numpy()
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                if printed[j] - printed[i] > 0.01 and vals[j] < vals[i] - 0.01:
                    same_state_violations += 1
    assert same_state_violations == 0
    return {
        "excluded_cognate_contact": frame.cognate_guard,
        "independent_vs_primary": primary_error,
        "same_state_score_order_inversions_gt_point01": same_state_violations,
        "raw_archive_sha256": sha(path.read_bytes()),
        "successes": len(df),
        "scored_coordinates": 2 * len(df),
        "observed_gradient_calls": int(df.calls.sum()),
        "max_transformed_start_component_error_A": max(rigid_error),
        "original_oracle_id": oracle,
    }


def main(out, selection, compare_primary):
    wall, cpu = time.perf_counter(), time.process_time()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        frame = Frame(Path(tmp))
        result = {"proper_protein_fits": frame.fits}
        if selection in ("blind", "both"):
            result["blind"] = blind(frame, out, STUDY / "blind/results", compare_primary)
        if selection in ("diagnostic", "both"):
            result["diagnostic"] = diagnostic(frame, out, STUDY / "diagnostic")
    result.update(wall_s=time.perf_counter() - wall, cpu_s=time.process_time() - cpu)
    (out / f"{selection}_receipt.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--part", choices=["both", "blind", "diagnostic"], default="both")
    p.add_argument(
        "--raw-only",
        action="store_true",
        help="compute from records before comparison tables are available",
    )
    args = p.parse_args()
    main(args.output, args.part, not args.raw_only)
