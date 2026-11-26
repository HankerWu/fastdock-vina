"""Measure the identical local-window perturbation around archived raw PDBQT candidates."""

import argparse
import csv
import io
import itertools
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ANCHOR = HERE.parent / "mapping_ensemble/results"
CELL_KEYS = ("receptor", "ligand", "seed")


def table(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def write(path, rows):
    with path.open("w", newline="") as source:
        writer = csv.DictWriter(source, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def raw_models(path):
    items = []
    for model in path.read_text().split("ENDMDL")[:-1]:

        def extract(label):
            return float(
                re.search(
                    r"^REMARK " + re.escape(label) + r":\s+([-0-9.]+)", model, re.M
                )[1]
            )

        items.append((extract("VINA RESULT"), extract("INTER"), extract("INTRA")))
    return items


def jobs(panel, base):
    """Each tuple carries a frozen receptor, ligand topology, pose pool, and per-model experimental references."""
    if panel == "hcaii":
        prepared = base / "fdv-second/prepared"
        runs = base / "fdv-map-states-declared"
        data = json.loads((prepared / "metadata.json").read_text())
        by_pool = {}
        for r in table(runs / "candidates.csv"):
            if r["slot"] in ("common", "second_reflected"):
                key = tuple(r[k] for k in (*CELL_KEYS, "slot"))
                by_pool.setdefault(key, []).append(r)
        for (rec, lig, seed, slot), rows in by_pool.items():
            name = f"{rec}_{lig}_{seed}_{slot}"
            source = (
                prepared / f"{lig}_ligand.pdbqt"
                if slot == "common"
                else runs / "inputs" / f"{lig}_reflected.pdbqt"
            )
            yield (
                (rec, lig, seed, "dual", slot, data[lig]["Ki_nM"]),
                prepared / f"{rec}_receptor.pdbqt",
                source,
                runs / (name + ".pdbqt"),
                base / "fdv-second/crystal_start" / f"{lig}.bonds",
                rows,
            )
    elif panel in ("brd4", "t4l", "alkyl"):
        if panel == "alkyl":
            prepared = base / "fdv-alkyl/prepared"
            runs = base / "fdv-alkyl/runs"
            labels = {
                "4W52": 175,
                "4W53": 102,
                "4W54": 68,
                "4W55": 18,
                "4W58": None,
                "4W59": None,
            }
            rows = table(HERE.parent / "alkyl_holdout/results/poses.csv")
        else:
            prepared = (
                base
                / ("fdv-brd4" if panel == "brd4" else "fdv-pro-fixtures")
                / "prepared"
            )
            runs = base / f"fdv-map-{panel}"
            labels = {
                k: v["kd_uM"]
                for k, v in json.loads((prepared / "manifest.json").read_text()).items()
            }
            rows = table(ANCHOR / panel / "modes.csv")
        groups = {}
        for r in rows:
            if r["arm"] not in ("native", "multi"):
                continue
            key = (
                r.get("receptor", r.get("receptor_state")),
                r["ligand"],
                r["seed"],
                r["arm"],
            )
            groups.setdefault(key, []).append(r)
        for (rec, lig, seed, arm), rs in groups.items():
            name = f"{rec}_{lig}_{seed}_{arm}"
            yield (
                (rec, lig, seed, arm, "same", labels[lig]),
                prepared / f"{rec}_receptor.pdbqt",
                prepared / f"{lig}_ligand.pdbqt",
                runs / (name + ".pdbqt"),
                None,
                rs,
            )
    elif panel in ("trypsin", "bovine_ca"):
        if panel == "trypsin":
            runs = prepared = base / "fdv-ionpair/trypsin"
            labels = {r["ligand"]: r["kd_uM"] for r in table(runs / "series.csv")}
            records = table(HERE.parent / "directional_ion_pair/results/modes.csv")
            rms = {
                (r["receptor"], r["ligand"], r["seed"], r["arm"], r["model"]): r[
                    "heavy_rms"
                ]
                for r in records
                if r["target"] == "trypsin"
            }
            receptors = ("3PTB", "1BTY")
            arms = ("native", "pair")
        else:
            runs = base / "fdv-proton/runs"
            prepared = base / "fdv-proton/prepared"
            labels = {
                r["ligand"]: r["measured_kd_nm"]
                for r in table(HERE.parent / "proton_linkage/fluorinated.csv")
            }
            receptors = ("1V9E", "6SKT")
            arms = ("native", "point")
            rms = {}
        for rec, lig, seed, arm in itertools.product(
            receptors, labels, ("314159", "271828"), arms
        ):
            pool = runs / f"{rec}_{lig}_{seed}_{arm}.pdbqt"
            ref = [
                {"rmsd": rms.get((rec, lig, seed, arm, str(i + 1)))}
                for i, _ in enumerate(raw_models(pool))
            ]
            receptor = prepared / (
                rec + (".pdbqt" if panel == "trypsin" else "_receptor.pdbqt")
            )
            yield (
                (rec, lig, seed, arm, "same", labels[lig]),
                receptor,
                prepared / (lig + ".pdbqt"),
                pool,
                None,
                ref,
            )
    else:
        labels = {
            r["ligand"]: r["ki_nm"]
            for r in table(HERE.parent / "affinity_partition/fxa_series.csv")
            if r["included"] == "1"
        }
        for rec, lig, seed in itertools.product(
            ("1FJS", "2P16"), labels, ("314159", "271828")
        ):
            runs = base / ("fdv-fxa1" if rec == "1FJS" else "fdv-fxa")
            arm = "native" if panel == "fxa" else "pair"
            path = (
                runs / f"{lig}_{seed}.pdbqt"
                if panel == "fxa"
                else base / "fdv-ionpair/fxa" / f"{rec}_{lig}_{seed}_pair.pdbqt"
            )
            yield (
                (rec, lig, seed, arm, "same", labels[lig]),
                runs / (rec + ".pdbqt"),
                runs / (lig + ".pdbqt"),
                path,
                None,
                [{} for _ in raw_models(path)],
            )


def sample(job, binary, count):
    (rec, lig, seed, arm, slot, label), receptor, source, pool, bonds, reference = job
    cmd = [
        binary,
        "--score_only",
        "--receptor",
        receptor,
        "--ligand",
        source,
        "--rigid_window_pool",
        pool,
        "--rigid_window_samples",
        str(count),
        "--cpu",
        "1",
        "--verbosity",
        "0",
    ]
    if bonds:
        cmd.extend(["--ligand_bonds", bonds])
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode:
        raise ValueError(f"{pool}: {proc.stderr}")
    header = "model,inter,intra,slope,"
    stats = list(csv.DictReader(io.StringIO(header + proc.stdout.split(header, 1)[1])))
    rows = []
    for values, (aff, inter, intra), ref in zip(
        stats, raw_models(pool), reference, strict=True
    ):
        if ref.get("model", ref.get("rank", values["model"])) != values["model"]:
            raise ValueError(f"Pose order mismatch for {pool}")
        rows.append(
            dict(
                receptor=rec,
                ligand=lig,
                seed=seed,
                arm=arm,
                slot=slot,
                measurement=label,
                rmsd=ref.get("rmsd"),
                sign_ok=ref.get("same_sign_as_bound", 1),
                source_affinity=aff,
                source_inter=float(ref.get("inter", inter)),
                source_intra=float(ref.get("intra", intra)),
                **values,
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "panel",
        choices=(
            "hcaii",
            "brd4",
            "t4l",
            "alkyl",
            "fxa",
            "fxa_pair",
            "trypsin",
            "bovine_ca",
        ),
    )
    p.add_argument(
        "base", type=Path, help="Directory containing the archived fdv-* trees"
    )
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--samples", type=int, default=2048)
    a = p.parse_args()
    with ThreadPoolExecutor(a.jobs) as workers:
        rows = list(
            itertools.chain.from_iterable(
                workers.map(
                    lambda j: sample(j, a.binary.resolve(), a.samples),
                    jobs(a.panel, a.base),
                )
            )
        )
    write(a.output, rows)
    print(
        a.panel,
        len(rows),
        max(abs(float(r["inter"]) - r["source_inter"]) for r in rows),
        flush=True,
    )
