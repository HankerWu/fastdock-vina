"""Re-score the identical baseline poses under both binaries; never optimize or align a ligand."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import itertools
import json
from pathlib import Path
import re
import subprocess
import tempfile

from experiment import csvout

LABELS = {
    "inter": "Final Intermolecular Energy",
    "intra": "Final Total Internal Energy",
    "side_rigid": "Flex   - Receptor",
    "lig_rigid": "Ligand - Receptor",
    "lig_side": "Ligand - Flex side chains",
    "isolated_score": "Estimated Free Energy of Binding",
}


def run(job):
    inp, poses, arms, key, results, scorer = job
    cohort, q, r, seed = key
    pool = poses / "raw" / ("-".join((cohort, q, r, seed, "baseline")) + ".pdbqt")
    meta = json.loads((inp / r / "metadata.json").read_text())
    output = []
    with tempfile.TemporaryDirectory() as folder:
        lig, side = [Path(folder) / x for x in ("ligand.pdbqt", "side.pdbqt")]

        def compare(origin, rank, rmsd):
            values = dict(
                cohort=cohort,
                query=q,
                receiver=r,
                seed=seed,
                scorer=scorer,
                source=origin,
                rank=rank,
                rmsd=rmsd,
            )
            for arm, binary in arms.items():
                cmd = [
                    binary,
                    "--scoring",
                    scorer,
                    "--receptor",
                    inp / r / "partition_rigid.pdbqt",
                    "--ligand",
                    lig,
                    "--flex",
                    side,
                    "--score_only",
                    "--cpu",
                    "1",
                ]
                for ax, val in zip("xyz", meta["center"]):
                    cmd.extend(("--center_" + ax, str(val), "--size_" + ax, "40"))
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode:
                    raise RuntimeError(proc.stdout + proc.stderr)
                for k, label in LABELS.items():
                    values[arm + "_" + k] = float(
                        re.search(re.escape(label) + r"\s*:\s*([-+.\d]+)", proc.stdout)[
                            1
                        ]
                    )
                # A single common reference and torsion multiplier preserve ordering at fixed ligand.
                values[arm + "_sum"] = values[arm + "_inter"] + values[arm + "_intra"]
            return values

        if q == r and seed == "62483":
            # Position the free generated conformer as a rigid body by centroid translation into the box.
            # No experimental orientation or torsion is transferred. Use the experimental side-chain input.
            text = (inp / q / "ligand.pdbqt").read_text().splitlines()
            import numpy as np

            heavy = [
                l
                for l in text
                if l.startswith(("ATOM  ", "HETATM"))
                and l[77:].strip() not in ("H", "HD")
            ]
            xyz = lambda l: np.array([float(l[i : i + 8]) for i in (30, 38, 46)])
            d = np.asarray(meta["center"]) - np.mean([xyz(l) for l in heavy], axis=0)
            lines = [
                l[:30] + "".join(f"{v:8.3f}" for v in xyz(l) + d) + l[54:]
                if l.startswith(("ATOM  ", "HETATM"))
                else l
                for l in text
            ]
            lig.write_text("\n".join(lines) + "\n")
            side.write_bytes((inp / r / "partition_flex.pdbqt").read_bytes())
            output.append(compare("same_input", 0, None))
        for rank, block in enumerate(pool.read_text().split("ENDMDL")[:-1], 1):
            first, rest = block.split("\nROOT\n", 1)[1].split("BEGIN_RES", 1)
            lig.write_text("ROOT\n" + first)
            side.write_text("BEGIN_RES" + rest)
            row = compare("baseline_pose", rank, results[rank])
            if scorer == "vina":
                # The published Vina objective retains the original box wall and double-precision pose.
                archived = float(
                    re.search(r"REMARK INTER \+ INTRA:\s*([-+.\d]+)", block)[1]
                )
                row["archived_sum"] = archived
                row["archived_plus_correction"] = (
                    archived + row["fixed_sum"] - row["baseline_sum"]
                )
            output.append(row)
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for f in ("inputs", "experiment", "baseline", "fixed"):
        p.add_argument(f, type=Path)
    p.add_argument("--jobs", type=int, default=1)
    a = p.parse_args()
    w, o = a.inputs.resolve(), a.experiment.resolve()
    arms = {k: getattr(a, k).resolve() for k in ("baseline", "fixed")}
    groups = {}
    for row in csv.DictReader((o / "modes.csv").open()):
        if row["arm"] == "baseline" and row["cohort"] != "assay":
            key = tuple(row[k] for k in ("cohort", "query", "receiver", "seed"))
            groups.setdefault(key, {})[int(row["rank"])] = float(row["rmsd"])
    jobs = [
        (w, o, arms, k, vals, s)
        for (k, vals), s in itertools.product(groups.items(), ("vina", "vinardo"))
    ]
    with ThreadPoolExecutor(a.jobs) as workers:
        rows = [v for group in workers.map(run, jobs) for v in group]
    csvout(o / "frozen.csv", rows)
    ranks = []
    for k, scorer in itertools.product(groups, ("vina", "vinardo")):
        candidates = [
            r
            for r in rows
            if r["source"] == "baseline_pose"
            and r["scorer"] == scorer
            and tuple(r[x] for x in ("cohort", "query", "receiver", "seed")) == k
        ]
        best = {
            a: min(candidates, key=lambda r: (r[a + "_sum"], r["rank"])) for a in arms
        }
        ranks.append(
            dict(
                cohort=k[0],
                query=k[1],
                receiver=k[2],
                seed=k[3],
                scorer=scorer,
                baseline_rank=best["baseline"]["rank"],
                fixed_rank=best["fixed"]["rank"],
                baseline_rmsd=best["baseline"]["rmsd"],
                fixed_rmsd=best["fixed"]["rmsd"],
            )
        )
        if scorer == "vina":
            anchor = min(
                candidates, key=lambda r: (r["archived_plus_correction"], r["rank"])
            )
            original = min(candidates, key=lambda r: (r["archived_sum"], r["rank"]))
            ranks[-1].update(
                published_baseline_rank=original["rank"],
                anchored_fixed_rank=anchor["rank"],
                published_baseline_rmsd=original["rmsd"],
                anchored_fixed_rmsd=anchor["rmsd"],
            )
    csvout(o / "frozen_rerank.csv", ranks)


if __name__ == "__main__":
    main()
