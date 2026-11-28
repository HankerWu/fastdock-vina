"""Graph-only huperzine docking in the two independently measured receptors."""
import argparse
import csv
import hashlib
import json
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BENCH = HERE.parents[1]
TRANSFER = BENCH / "water_rank_envelope/transfer"
ARCHIVE = BENCH / "receptor_amide_states/transfer/transfer_evidence.tar.xz"
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import run as chemistry

SEEDS = (1103501, 1103503, 1103507)
ARMS = ("official24k", "pro24k", "pro48k")
RECEIVERS = {"4EY6": "GNT", "4EY7": "E20"}
OFFICIAL_SHA = "f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def tick():
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return np.array([time.perf_counter(), time.process_time() + children.ru_utime + children.ru_stime])


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def pocket(rows, comp):
    ligand, = chemistry.dep_records(rows, comp)
    ca = {l[22:27] + l[17:20]: chemistry.coords(l) for l in rows
          if l.startswith("ATOM") and l[21] == "A" and l[12:16] == " CA "}
    return dict(ca=ca, center=np.mean(list(ligand.values()), axis=0))


def load_inputs(work):
    past = json.loads((TRANSFER / "manifest.json").read_text())
    if sha(ARCHIVE.read_bytes()) != past["source_archive_sha256"]:
        raise ValueError("Changed receptor archive")
    generated = chemistry.ligand(TRANSFER / "HUP.cif", [], work / "HUP.pdbqt")
    hup = (work / "HUP.pdbqt").read_bytes()
    if generated["mapped"] or hup != (TRANSFER / "HUP.pdbqt").read_bytes() or sha(hup) != past["graph_only_pdbqt_sha256"]:
        raise ValueError("Graph-only HUP preparation changed")
    identities = {"receptor_archive": past["source_archive_sha256"], "HUP.pdbqt": sha(hup),
                  "HUP.cif": sha((TRANSFER / "HUP.cif").read_bytes())}
    centers, pockets = {}, {}
    with tarfile.open(ARCHIVE) as archive:
        for pdb, ccd in RECEIVERS.items():
            def get(member):
                data = archive.extractfile("novel/" + member).read()
                identities[member] = sha(data)
                return data
            dep = get(f"raw/{pdb}.pdb")
            info = json.loads(get(f"{pdb}/preparation.json"))
            receptor = get(f"{pdb}/receptor.pdbqt")
            if sha(dep) != info["pdb_sha256"] or info["receptor_status"] != "ready":
                raise ValueError(f"Unverified receptor {pdb}")
            (work / f"{pdb}.pdbqt").write_bytes(receptor)
            pockets[pdb] = pocket(dep.decode().splitlines(), ccd)
            centers[pdb] = info["center"]
            np.testing.assert_allclose(pockets[pdb]["center"], centers[pdb], rtol=0, atol=1e-10)
    rot, trans, residual, count = chemistry.alignment(pockets["4EY6"], pockets["4EY7"])
    centers["4EY7_mapped"] = (pockets["4EY6"]["center"] @ rot + trans).tolist()
    return dict(hashes=identities, centers_A=centers, noncognate_fit=dict(
        shared_ca=count, rmsd_A=residual, rotation=rot.tolist(), translation_A=trans.tolist()))


def cli(seed, arm, center):
    args = ["--receptor", "receptor.pdbqt", "--ligand", "HUP.pdbqt", "--cpu", "1", "--seed", str(seed),
            "--exhaustiveness", "3", "--max_evals", "48000" if arm == "pro48k" else "24000",
            "--num_modes", "9", "--energy_range", "100"]
    if arm != "official24k":
        args.append("--search_counters")
    args += ["--out", "dry.pdbqt"]
    for axis, value in zip("xyz", center, strict=True):
        args += ["--center_" + axis, str(value), "--size_" + axis, "24"]
    return args


def run_job(exe, args, work, directory):
    directory.mkdir(parents=True, exist_ok=True)
    pose = work / "dry.pdbqt"
    pose.unlink(missing_ok=True)
    save_json(directory / "command.json", ["vina", *args])
    clock = tick()
    run = subprocess.run([str(exe), *args], cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall, cpu = tick() - clock
    (directory / "dry.log").write_bytes(run.stdout)
    if pose.exists():
        shutil.copyfile(pose, directory / "dry.pdbqt")
    digest = sha(pose.read_bytes()) if pose.exists() else None
    text = run.stdout.decode(errors="replace")
    totals = {}
    for key in ("hunt_calls", "hunt_evals", "second_calls", "second_evals", "explicit_calls", "explicit_evals"):
        nums = re.findall(r"\b" + key + r"=(\d+)", text)
        totals[key] = sum(map(int, nums)) if nums else None
    return dict(status="ok" if run.returncode == 0 and digest else "failed", returncode=run.returncode,
                wall_s=wall, cpu_s=cpu, output_sha256=digest, **totals)


def main(out, official, pro):
    out.mkdir(parents=True, exist_ok=False)
    clock = tick()
    binaries = {"official": sha(official.read_bytes()), "pro": sha(pro.read_bytes())}
    if binaries["official"] != OFFICIAL_SHA:
        raise ValueError("Expected release executable")
    runs = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        start = tick()
        inputs = load_inputs(work)
        inputs["prepare_wall_cpu_s"] = (tick() - start).tolist()
        save_json(out / "inputs.json", inputs)
        for rec in RECEIVERS:
            shutil.copyfile(work / (rec + ".pdbqt"), work / "receptor.pdbqt")
            for arm in (*ARMS, *(("pro24k_mapped",) if rec == "4EY7" else ())):
                center_key = "4EY7_mapped" if arm.endswith("_mapped") else rec
                for seed in SEEDS:
                    folder = out / "raw" / rec / f"{arm}_{seed}"
                    rec_arm = "official" if arm == "official24k" else "pro"
                    exe = official if rec_arm == "official" else pro
                    row = run_job(exe, cli(seed, arm, inputs["centers_A"][center_key]), work, folder)
                    old_match = None
                    if rec == "4EY6" and seed == SEEDS[0]:
                        old_file = TRANSFER / ("" if arm == "pro24k" else arm) / "dry.pdbqt"
                        old_match = row["output_sha256"] == sha(old_file.read_bytes())
                    runs.append(dict(receiver=rec, arm=arm, seed=seed, **row, historical_match=old_match))
                    table(out / "runs.csv", runs)
    save_json(out / "search_receipt.json", dict(binary_sha256=binaries, total_wall_cpu_s=(tick()-clock).tolist(),
                                                successful=sum(row["status"] == "ok" for row in runs),
                                                attempted=len(runs)))
    if any(row["status"] != "ok" for row in runs):
        raise RuntimeError("One or more registered jobs failed; inspect raw logs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("official", type=Path)
    parser.add_argument("pro", type=Path)
    args = parser.parse_args()
    main(args.out, args.official.resolve(), args.pro.resolve())
