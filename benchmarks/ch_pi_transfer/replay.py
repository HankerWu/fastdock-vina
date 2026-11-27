"""Score native galectin and glucose pools. Chemistry and labels stay separate."""

import argparse
import resource
import sys
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
import pandas as pd

from model import WEIGHT, CHPi

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
sys.path.insert(0, str(BENCH / "galectin_affinity"))
from sugar import crystal_core, ligand_core, xyz

TARGETS = {
    "1GZW": "Gal1",
    "4Q26": "Gal1",
    "1KJL": "Gal3",
    "3T1L": "Gal3",
    "1V4S": "GCK",
    "3IDH": "GCK",
    "4ZW9": "GLUT3",
}
ARMS = ("native", "ligand_CH", "protein_CH", "bidirectional")


def blocks(path):
    data = path.read_text().strip()
    if not (data.startswith("MODEL ") and data.endswith("ENDMDL")):
        raise ValueError(f"incomplete pose pool: {path}")
    return data.split("ENDMDL")[:-1]


def printed_score(block):
    fields = [ln.split() for ln in block.splitlines() if ln.startswith("REMARK ")]
    (score,) = [float(f[3]) for f in fields if f[1:3] == ["VINA", "RESULT:"]]
    (total,) = [float(f[-1]) for f in fields if f[1:4] == ["INTER", "+", "INTRA:"]]
    (unbound,) = [float(f[-1]) for f in fields if f[1] == "UNBOUND:"]
    return score, total - unbound


def archived_scale(modes):
    # Historical engines can count half rotors. Identify the common multiplier
    # and independently replay *every* rounded pool score to file precision.
    data = np.asarray([printed_score(mode) for mode in modes])
    scores, raw = data.T
    i = np.argmax(np.abs(scores))
    nrot = round(2 * (raw[i] / scores[i] - 1) / 0.05846) / 2
    scale = 1.0 / (1 + 0.05846 * nrot)
    if nrot < 0 or np.max(np.abs(scores - raw * scale)) > 0.00151:
        raise ValueError("cannot reconstruct the archived Vina torsional normalization")
    return scale


def run(gal, glu, output):
    output.mkdir(exist_ok=True, parents=True)
    old = pd.read_csv(BENCH / "ligand_hydroxyl/glucose/results/poses.csv")
    glucose = old[old.engine == "pro"].set_index(["case", "cap", "seed", "rank"])
    recs, ligs, models = {}, {}, {}
    coords = {p: crystal_core(gal, p) for p in TARGETS if TARGETS[p].startswith("Gal")}
    ids = {str(i): ligand_core(gal, i) for i in range(5, 18)}
    records = []
    paths = sorted((gal / "pools").glob("*_native.pdbqt")) + sorted(
        (glu / "runs").glob("*_pro.pdbqt")
    )
    if len(paths) != 224:
        raise ValueError(
            f"expected 208 native galectin + 16 glucose pools, got {len(paths)}"
        )
    wall, cpu = perf_counter(), process_time()
    setup = calc = 0.0
    for path in paths:
        if path.parent.name == "pools":
            pdb, compound, seed, _ = path.stem.split("_")
            cap = 25000
            rpath = gal / "prepared" / f"{pdb}.pdbqt"
            lpath = gal / "prepared" / f"{compound}.pdbqt"
        else:
            case, cap, seed, _ = path.stem.split("_")
            pdb, compound = case.split("-")
            rpath = glu / pdb / "receptor.pdbqt"
            lpath = glu / pdb / f"{compound}.pdbqt"
        start = perf_counter()
        if rpath not in recs:
            recs[rpath] = rpath.read_text()
        if lpath not in ligs:
            ligs[lpath] = lpath.read_text()
        if (rpath, lpath) not in models:
            models[rpath, lpath] = CHPi(recs[rpath], ligs[lpath])
        model = models[rpath, lpath]
        setup += perf_counter() - start
        modes = blocks(path)
        factor = archived_scale(modes)
        for rank, mode in enumerate(modes, 1):
            start = perf_counter()
            forward, reverse = model.components(mode)
            calc += perf_counter() - start
            native, _ = printed_score(mode)
            if pdb in coords:
                atoms = {
                    int(l[6:11]): xyz(l)
                    for l in mode.splitlines()
                    if l.startswith(("ATOM  ", "HETATM"))
                }
                diff = np.asarray([atoms[i] for i in ids[compound]]) - coords[pdb]
                rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))
            else:
                ref = glucose.loc[case, int(cap), int(seed), rank]
                rmsd = ref.rmsd
                if abs(ref.scale - factor) > 1e-8 or abs(ref.native - native) > 1e-6:
                    raise ValueError("glucose label or factor does not match archive")
            normalized = WEIGHT * factor
            records.append(
                {
                    "target": TARGETS[pdb],
                    "pdb": pdb,
                    "compound": compound,
                    "cap": cap,
                    "seed": seed,
                    "mode": rank,
                    "rmsd": rmsd,
                    "metric": "core6" if pdb in coords else "heavy",
                    "factor": factor,
                    "forward_raw": forward,
                    "reverse_raw": reverse,
                    "native": native,
                    "ligand_CH": native + normalized * forward,
                    "protein_CH": native + normalized * reverse,
                    "bidirectional": native + normalized * (forward + reverse),
                }
            )
    pd.DataFrame(records).to_csv(output / "poses.csv", index=False, float_format="%.9f")
    pd.DataFrame(
        [
            {
                "pools": len(paths),
                "modes": len(records),
                "core_setup_s": setup,
                "core_scoring_s": calc,
                "wall_s": perf_counter() - wall,
                "cpu_s": process_time() - cpu,
                "peak_rss_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                / 1024,
            }
        ]
    ).to_csv(output / "cost.csv", index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("gal", "glucose", "out"):
        p.add_argument(key, type=Path)
    args = p.parse_args()
    run(args.gal, args.glucose, args.out)
