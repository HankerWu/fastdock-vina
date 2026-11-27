"""Replay the two preregistered assay panels with frozen input topology and native top poses."""
import argparse
import csv
import lzma
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import re
from statistics import mean, median
import subprocess
import tarfile
from tempfile import TemporaryDirectory

from scipy.stats import kendalltau, spearmanr

BENCH = Path(__file__).resolve().parents[1]
WEIGHT = 0.05846


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def save(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def archive_files(path, pattern):
    with tarfile.open(path, "r|xz") as tar:
        return {m.name: tar.extractfile(m).read() for m in tar if m.isfile() and re.fullmatch(pattern, m.name)}


def assay_inputs():
    gal = BENCH / "pair_gradient/galectin_external"
    prep = archive_files(gal / "raw.tar.xz", r"prepared/.+\.pdbqt")
    pools = archive_files(BENCH / "silent_torsion/affinity/primary/raw.tar.xz", r"runs/.+_parent\.pdbqt")
    for row in read(gal / "manifest.csv"):
        key = row["compound"]
        for rec in ("1KJL", "3T1L"):
            yield ("Gal3", key, float(row["kd_uM"]), rec, 214039,
                   prep[f"prepared/{rec}.pdbqt"], prep[f"prepared/{key}.pdbqt"], pools[f"runs/{rec}_{key}_parent.pdbqt"])
    fxa = archive_files(BENCH / "frozen_core/affinity/fxa_frozen.tar.xz", r".+\.pdbqt")
    for row in read(BENCH / "affinity_partition/fxa_series.csv"):
        if row["included"] != "1":
            continue
        key = row["ligand"]
        for rec in ("2P16", "1FJS"):
            for seed in (271828, 314159):
                yield ("FXa", key, float(row["ki_nm"]), rec, seed,
                       fxa[f"{rec}/{rec}.pdbqt"], fxa[f"{rec}/{key}.pdbqt"], fxa[f"{rec}/{key}_{seed}.pdbqt"])


    s6k = BENCH / "affinity_reference_ensemble/results"
    files = archive_files(s6k / "s6k_evidence.tar.xz", r"(?:prepared/[^/]+|pools/[^/]+_pro)\.pdbqt")
    for row in read(s6k / "s6k_panel.csv"):
        if row["decision"] != "included":
            continue
        key = row["cid"]
        for seed in (271828, 314159):
            yield ("S6K1", key, float(row["kd_nM"]), "3A60", seed,
                   files["prepared/3A60.pdbqt"], files[f"prepared/{key}.pdbqt"], files[f"pools/3A60_{key}_{seed}_pro.pdbqt"])


def terms(binary, receptor, ligand, pose):
    output = subprocess.check_output([binary, "--cpu", "1", "--seed", "214039", "--score_only",
                                      "--receptor", receptor, "--ligand", ligand, "--pose_coordinates", pose], text=True)
    output = output.split("Estimated Free Energy of Binding")[-1]
    values = re.findall(r":\s+(-?\d+\.\d+)\s+\(kcal/mol\)", output)
    if len(values) != 10:
        raise ValueError(output)
    return list(map(float, values))


def count(probe, ligand):
    out = subprocess.check_output([probe, ligand], text=True)
    return float(next(csv.DictReader(out.splitlines()))["num_tors"])


def replay(args, work):
    cells, counts = [], {}
    for target, key, measured, rec, seed, receptor_bytes, ligand_bytes, pool in assay_inputs():
        receptor, ligand, pose = [work / name for name in ("receptor.pdbqt", "ligand.pdbqt", "pose.pdbqt")]
        receptor.write_bytes(receptor_bytes)
        ligand.write_bytes(ligand_bytes)
        first = pool.decode().split("ENDMDL")[0]
        pose.write_text(first + "ENDMDL\n")
        if ligand_bytes not in counts:
            counts[ligand_bytes] = tuple(count(p, ligand) for p in [args.parent_probe, args.fixed_probe])
        old, new = counts[ligand_bytes]
        cli = [terms(b, receptor, ligand, pose) for b in (args.parent, args.fixed)]
        if cli[0][1:8] != cli[1][1:8] or cli[0][9] != cli[1][9]:
            raise ValueError(f"Potential or reference changed for {key}")
        # The first native mode supplies its own unbound intra. Thus its numerator is INTER.
        def remark(term):
            return float(re.search(r"^REMARK " + re.escape(term) + r":\s+(-?[\d.]+)", first, re.M)[1])
        if remark("INTRA") != remark("UNBOUND"):
            raise ValueError(f"Archived first mode used another reference: {key}")
        inter = remark("INTER")
        a, b = (inter / (1 + WEIGHT * n) for n in (old, new))
        if abs(a - remark("VINA RESULT")) > 0.001 or any(
            abs(e[0] - e[1] / (1 + WEIGHT * n)) > 0.001 for e, n in zip(cli, (old, new))
        ):
            raise ValueError(f"Unreconciled scaling for {key}")
        cells.append(dict(target=target, ligand=key, measurement=measured, receptor=rec, seed=seed,
                          old_rotors=old, fixed_rotors=new, saved_score=remark("VINA RESULT"),
                          saved_inter=inter, old=a, fixed=b, fixed_xyz_inter=cli[0][1],
                          fixed_xyz_intra=cli[0][4], cli_old=cli[0][0], cli_fixed=cli[1][0]))
    assert len(cells) == 234 and len({(r["target"], r["ligand"]) for r in cells}) == 89
    return cells


def aggregate(cells, field):
    scores = defaultdict(lambda: defaultdict(list))
    for row in cells:
        scores[row["ligand"]][row["receptor"]].append(row[field])
    return {key: mean(median(vals) for vals in recs.values()) for key, recs in scores.items()}


def credit(diff):
    return 1.0 if diff > 0 else (0.0 if diff < 0 else 0.5)


def summarize(cells):
    summaries, pairs = [], []
    for target in ("Gal3", "FXa", "S6K1"):
        panel = [r for r in cells if r["target"] == target]
        for rec in ["mean", *sorted({r["receptor"] for r in panel})]:
            sub = [r for r in panel if rec == "mean" or r["receptor"] == rec]
            source = {r["ligand"]: r["measurement"] for r in sub}
            for analysis, fields in (("saved_inter", ("old", "fixed")), ("cli_rounded_xyz", ("cli_old", "cli_fixed"))):
                predicted = [aggregate(sub, f) for f in fields]
                assay_pairs = []
                for x, y in combinations(sorted(source, key=source.get), 2):
                    if source[x] == source[y]:
                        continue
                    gaps = [v[y] - v[x] for v in predicted]
                    row = dict(target=target, receptor=rec, analysis=analysis, strong=x, weak=y,
                               threefold=int(source[y] >= 3 * source[x] - 1e-12), old_credit=credit(gaps[0]),
                               fixed_credit=credit(gaps[1]), old_gap=gaps[0], fixed_gap=gaps[1])
                    assay_pairs.append(row)
                pairs.extend(assay_pairs)
                for subset, group in (("all", assay_pairs), ("threefold", [p for p in assay_pairs if p["threefold"]])):
                    changes = [p["fixed_credit"] - p["old_credit"] for p in group]
                    deltas = [sum(p["fixed_credit"] - p["old_credit"] for p in group if key not in (p["strong"], p["weak"])) for key in source]
                    summary = dict(target=target, receptor=rec, analysis=analysis, subset=subset,
                                   compounds=len(source), pairs=len(group), old_credit=sum(p["old_credit"] for p in group),
                                   fixed_credit=sum(p["fixed_credit"] for p in group), gain=sum(d > 0 for d in changes),
                                   loss=sum(d < 0 for d in changes), unchanged=changes.count(0),
                                   leave_one_min=min(deltas), leave_one_max=max(deltas),
                                   pairs_within_0p001=sum(min(abs(p["old_gap"]), abs(p["fixed_gap"])) <= .001 for p in group))
                    # Rank correlations are compound level, only defined here on the entire cohort.
                    for label, values in zip(("old", "fixed"), predicted):
                        y = [values[k] for k in source]
                        summary[label + "_rho"] = spearmanr(list(source.values()), y).statistic if subset == "all" else ""
                        summary[label + "_tau"] = kendalltau(list(source.values()), y).statistic if subset == "all" else ""
                    summaries.append(summary)
    return summaries, pairs


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("parent", "fixed", "parent_probe", "fixed_probe", "out"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        data = replay(a, Path(tmp))
    stats, pairs = summarize(data)
    save(a.out / "poses.csv", data)
    save(a.out / "summary.csv", stats)
    raw = a.out / "pairs.csv"
    save(raw, pairs)
    raw.with_suffix(".csv.xz").write_bytes(lzma.compress(raw.read_bytes()))
    raw.unlink()
