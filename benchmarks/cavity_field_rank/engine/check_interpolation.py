"""Independent SciPy interpolation of exported native maps on every frozen pose.

Usage: python check_interpolation.py POSE_WORK PRO_BINARY CHECK_WORK SUMMARY.json
The Vina .map writer rounds field samples to four significant digits. Use the
full precision center from the original run because .map also rounds its header.
"""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.stats import spearmanr


# src/lib/atom_constants.h. These panels use conventional organic ligand types.
CHANNELS = "C_H C_P N_P N_D N_A N_DA O_P O_D O_A O_DA S_P P_P F_H Cl_H Br_H I_H Si At Met_D".split()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cpu():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + child.ru_utime + child.ru_stime


def read_maps(work, cfg, binary, dest):
    dest.mkdir(exist_ok=True, parents=True)
    cmd = [binary, "--ligand", work/cfg["ligand"], "--receptor", work/cfg["receptor"],
           "--cpu", 1, "--seed", 1, "--randomize_only", "--out", dest/"unused_randomized_probe.pdbqt",
           "--write_maps", dest/"channel"]
    for a, center in zip("xyz", cfg["center"], strict=True):
        cmd.extend(["--center_" + a, center, "--size_" + a, cfg["size"]])
    with (dest/"export.log").open("w") as log:
        subprocess.run(list(map(str, cmd)), check=True, stdout=log, stderr=subprocess.STDOUT)
    maps = {}
    for path in dest.glob("channel.*.map"):
        with path.open() as f:
            header = [next(f).split() for _ in range(6)]
            step = float(header[3][1])
            nodes = np.array(header[4][1:], int) + 1
            samples = np.loadtxt(f)
        maps[path.name.split(".")[1]] = samples.reshape(nodes[::-1]).T
    return maps, step, nodes


def score(coords, types, maps, step, low, nodes):
    lattice = (coords - low) / step
    clipped = np.clip(lattice, 0, nodes-1)
    wall = np.abs(lattice-clipped).sum(axis=(1, 2)) * step * 1e6
    energies = np.empty(coords.shape[:2])
    for t in np.unique(types):
        mask = types == t
        query = clipped[:, mask].reshape(-1, 3).T
        values = map_coordinates(maps[CHANNELS[t]], query, order=1, mode="nearest", prefilter=False)
        energies[:, mask] = values.reshape(len(coords), -1)
    repulsive = energies > 0
    energies[repulsive] /= 1 + energies[repulsive] / 10
    return energies.sum(axis=1) + wall


def run(work, binary, tmp):
    start, start_cpu = time.perf_counter(), cpu()
    rows, map_sha = [], {}
    configs = json.loads((work/"inputs.json").read_text())["cells"]
    for name, cfg in configs.items():
        dest = tmp / name
        maps, step, nodes = read_maps(work, cfg, binary, dest)
        low = np.array(cfg["center"]) - (nodes-1) * step / 2
        map_sha[name] = {p.name: sha(p) for p in sorted(dest.glob("channel.*.map"))}
        for source in sorted((work/"banks").glob(name+"_*.source.npz")):
            stem = source.name.removesuffix(".source.npz")
            with np.load(source) as inp, np.load(source.with_name(stem+".npz")) as saved:
                serial = inp["serials"]
                x = inp["xyz"]
                coords = np.einsum("kaj,hj->kha", saved["matched_rotation"], x-x.mean(axis=0))
                coords += saved["center"][:, None, :]
            idx = np.argsort(serial)
            for arm in ["true", "perm"]:
                root = work / "raw" / f"{stem}_{arm}"
                mapping = np.array([line.split()[1:] for line in root.with_suffix(".inject").read_text().splitlines()
                                    if line.startswith("field_type ")], int)
                np.testing.assert_array_equal(mapping[:, 0], serial[idx])
                expected = np.roll(mapping[:,1], -(len(serial)//2)) if arm == "perm" else mapping[:,1]
                np.testing.assert_array_equal(expected, mapping[:, 2])
                independent = score(coords[:, idx], expected, maps, step, low, nodes)
                native = np.genfromtxt(root.with_suffix(".csv"), delimiter=",", names=True)
                np.testing.assert_array_equal(native["index"], np.arange(1024))
                winner, original = int(independent.argmin()), int(native["score"].argmin())
                rows.append(dict(case=stem, arm=arm, n=len(independent), independent_winner=winner,
                                 native_winner=original, max_abs_difference=float(np.abs(independent-native["score"]).max()),
                                 rms_difference=float(np.mean((independent-native["score"])**2)**0.5),
                                 score_rho=float(spearmanr(independent,native["score"]).statistic),
                                 changed_types=int((expected!=mapping[:,1]).sum()), heavy_atoms=len(serial)))
    return dict(rows=rows, map_sha256=map_sha, pro_sha256=sha(binary), wall_s=time.perf_counter()-start,
                cpu_s=cpu()-start_cpu)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ["pose_work","pro_binary","check_work","summary"]:
        p.add_argument(k, type=lambda arg: Path(arg).resolve())
    a = p.parse_args()
    result = run(a.pose_work, a.pro_binary, a.check_work)
    a.summary.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    if any(r["independent_winner"] != r["native_winner"] for r in result["rows"]):
        raise SystemExit("exported-map winner mismatch; see the summary")
