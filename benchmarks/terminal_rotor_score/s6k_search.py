"""End-to-end S6K1 transfer under the archived two- and four-trajectory budgets."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import re
from statistics import median
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory

from scipy.stats import kendalltau, spearmanr

from assay import BENCH, archive_files, credit, read, save


def run(args, work):
    bench = BENCH / "affinity_reference_ensemble/results"
    inputs = archive_files(bench / "s6k_evidence.tar.xz", r"prepared/(?:[^/]+\.pdbqt|3A60.json)")
    for name, content in inputs.items():
        (work / Path(name).name).write_bytes(content)
    center = json.loads((work / "3A60.json").read_text())["center"]
    records = {r["cid"]: float(r["kd_nM"]) for r in read(bench / "s6k_panel.csv") if r["decision"] == "included"}
    raw = work / "raw"
    raw.mkdir()

    def process(job):
        key, seed, nt = job
        out = []
        for arm in (("old", "fixed") if (int(key) + seed) % 2 == 0 else ("fixed", "old")):
            binary = args.parent if arm == "old" else args.fixed
            stem = raw / f"{key}_{seed}_{nt}_{arm}"
            pool, log, cost = [stem.with_suffix(s) for s in (".pdbqt", ".log", ".cost")]
            cmd = [binary, "--receptor", work / "3A60.pdbqt", "--ligand", work / f"{key}.pdbqt",
                   "--cpu", "1", "--exhaustiveness", str(nt), "--max_evals", "25000", "--num_modes", "9",
                   "--energy_range", "100", "--seed", str(seed), "--out", pool, "--pose_selection", "late", "--search_counters"]
            for axis, v in zip("xyz", center):
                cmd += ["--center_" + axis, str(v), "--size_" + axis, "26"]
            subprocess.run([sys.executable, "-S", BENCH / "rigid_flip/process_cost.py", log, cost, *cmd], check=True)
            text = pool.read_text()
            xyz = "\n".join(s for s in text.splitlines() if s.startswith(("ATOM ", "HETATM", "ENDMDL")))
            trace = "\n".join(s for s in log.read_text().splitlines() if s.startswith(("Search trajectory:", "Search counters:", "Native merge:", "Hunt chain:", "Pocket start:")))
            gradient_calls = sum(map(int, re.findall(r"\b(?:hunt_evals|second_evals|explicit_evals)=(\d+)", trace)))
            assert gradient_calls > 0, stem
            out.append(dict(compound=key, seed=seed, trajectories=nt, arm=arm, kd_nM=records[key],
                            first=float(re.search(r"REMARK VINA RESULT:\s+(-?[\d.]+)", text)[1]),
                            evaluations=gradient_calls, xyz_sha256=hashlib.sha256(xyz.encode()).hexdigest(),
                            trace_sha256=hashlib.sha256(trace.encode()).hexdigest(), **json.loads(cost.read_text())))
        assert out[0]["xyz_sha256"] == out[1]["xyz_sha256"] and out[0]["trace_sha256"] == out[1]["trace_sha256"], job
        return out

    with ThreadPoolExecutor(args.workers) as ex:
        rows = [r for result in ex.map(process, ((k, s, n) for k in records for s in (271828, 314159) for n in (2, 4))) for r in result]
    save(args.out / "s6k_search.csv", rows)
    upstream = [dict(compound=r["ligand"], first=float(r["score"]), arm="upstream_archived", trajectories=2)
                for r in read(bench / "s6k_docks.csv") if r["arm"] == "upstream"]
    stats = []
    for arm, nt in sorted({(r["arm"], r["trajectories"]) for r in rows + upstream}):
        cells = [r for r in rows + upstream if r["arm"] == arm and r["trajectories"] == nt]
        pred = {k: median(r["first"] for r in cells if r["compound"] == k) for k in records}
        pairs = [(k, j) for k in records for j in records if records[j] >= 3 * records[k]]
        stats.append(dict(arm=arm, trajectories=nt, compounds=len(records), pairs=len(pairs),
                          correct=sum(credit(pred[j] - pred[k]) for k, j in pairs),
                          rho=spearmanr(list(records.values()), list(pred.values())).statistic,
                          tau=kendalltau(list(records.values()), list(pred.values())).statistic,
                          evals=sum(r["evaluations"] for r in cells) if arm != "upstream_archived" else "",
                          cpu_s=sum(r["cpu_s"] for r in cells) if arm != "upstream_archived" else "",
                          wall_s=sum(r["wall_s"] for r in cells) if arm != "upstream_archived" else "",
                          peak_rss_KiB=max(r["peak_rss_KiB"] for r in cells) if arm != "upstream_archived" else ""))
    save(args.out / "s6k_search_summary.csv", stats)
    with tarfile.open(args.out / "s6k_search.tar.xz", "w:xz") as archive:
        archive.add(raw, arcname="raw")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("parent", "fixed", "out"):
        parser.add_argument(arg, type=Path)
    parser.add_argument("--workers", type=int, default=2)
    a = parser.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as d:
        run(a, Path(d))
