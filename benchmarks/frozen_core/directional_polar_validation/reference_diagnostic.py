"""Reference-start diagnostic, supplementary to the preregistered unguided predictions."""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "paired_crossdocking"))
import crystal_start


def energy_sum(vina, ligand, receptor, pose):
    """Keep query topology fixed; evaluate intermolecular and internal energy at supplied coordinates."""
    result = subprocess.check_output([vina, "--receptor", receptor, "--ligand", ligand, "--pose_coordinates", pose,
                                      "--score_only", "--cpu", "1", "--seed", "43103"], text=True)
    return sum(float(re.search(name + r"\s*:\s*([-\d.]+)", result)[1])
               for name in ("Final Intermolecular Energy", "Final Total Internal Energy"))


def main(work, vina):
    system = json.loads((work / "sources.json").read_text())["systems"]
    meta = {pdb: json.loads((work / pdb / "preparation.json").read_text())
            for pdb, row in system.items() if row["receptor_status"] == "ready"}
    (work / "metadata.json").write_text(json.dumps(meta))
    crystal_start.run(work, vina, pairs={(p, p) for p in meta}, compare_mappings=True)
    rows = []
    with TemporaryDirectory() as d:
        tmp = Path(d) / "pose.pdbqt"
        for pdb in meta:
            ligand, receptor = (work / pdb / (n + ".pdbqt") for n in ("ligand", "receptor"))
            ref = work / "crystal_diagnostic" / f"{pdb}-{pdb}-best_projected.pdbqt"
            query_e = energy_sum(vina, ligand, receptor, ref)
            for pool in sorted((work / "runs").glob(f"{pdb}_*_pro.pdbqt")):
                scores = []
                for block in pool.read_text().split("ENDMDL")[:-1]:
                    tmp.write_text(block.strip() + "\nENDMDL\n")
                    scores.append(energy_sum(vina, ligand, receptor, tmp))
                _, budget, seed, _ = pool.stem.split("_")
                rows.append(dict(pdb=pdb, cap=budget, seed=seed,
                                 projected_inter_plus_intra=query_e, best_blind_inter_plus_intra=min(scores),
                                 projected_minus_best_blind=query_e-min(scores),
                                 insertion_rank=1+sum(e < query_e for e in scores), modes=len(scores)))
    with (work / "reference_vs_pool.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("vina", type=Path)
    a = p.parse_args()
    main(a.work.resolve(), a.vina.resolve())
