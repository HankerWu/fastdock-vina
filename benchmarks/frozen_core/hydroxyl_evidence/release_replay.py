"""Re-run the new blind panels with the production topology correction on identical inputs."""
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent


def pools(study, root):
    if study == "esr2":
        for pose in sorted(root.glob("*/*-*-*.pdbqt")):
            pdb, budget, seed = pose.stem.split("-")
            yield pdb, budget, seed, pose
    else:
        for pose in sorted((root/"runs").glob("*_pro.pdbqt")):
            pdb, budget, seed, _ = pose.stem.split("_")
            yield pdb, budget, seed, pose


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vina", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = []
    panels = {
        "esr2": (("1QKM", "5TOA"), (12000, 60000), (52711, 52727)),
        "first_stage": (("4I1R", "2YIY"), (3000, 12000), (43103, 43117)),
        "ppar": (("2PRG", "3TKM"), (3000, 12000), (43103, 43117)),
    }
    binary = args.vina.resolve()
    for study in ("esr2","first_stage","ppar"):
        if study == "esr2":
            archive = HERE.parent / "hydroxyl_partition/esr2.tar.xz"
        else:
            archive = HERE.parent / "directional_polar_validation" / study / "raw.tar.xz"
        with tempfile.TemporaryDirectory() as work:
            root = Path(work)
            with tarfile.open(archive) as t:
                t.extractall(root, filter="data")
            out = root/"repeat.pdbqt"
            for pdb, budget, seed, expected in pools(study,root):
                site = root/pdb
                meta = json.loads((site/"preparation.json").read_text())
                cmd = [str(binary), "--receptor",str(site/"receptor.pdbqt"),"--ligand",str(site/"ligand.pdbqt"),
                       "--cpu","1","--exhaustiveness","4","--num_modes","9","--energy_range","100",
                       "--max_evals",budget,"--seed",seed,"--out",str(out)]
                if study != "esr2":
                    cmd.append("--search_counters")
                for ax,coord in zip("xyz",meta["center"]):
                    cmd.extend(["--center_"+ax,str(coord),"--size_"+ax,"24"])
                subprocess.run(cmd,capture_output=True,check=True)
                old,new = expected.read_bytes(),out.read_bytes()
                rows.append(dict(study=study,pdb=pdb,cap=budget,seed=seed,modes=new.count(b"ENDMDL"),
                                 archived_sha256=sha256(old).hexdigest(),release_sha256=sha256(new).hexdigest(),
                                 identical=old==new))
    from itertools import product
    expected_cells = {(study,pdb,cap,seed) for study, factors in panels.items()
                      for pdb,cap,seed in product(*factors)}
    observed_cells = {(r["study"],r["pdb"],int(r["cap"]),int(r["seed"])) for r in rows}
    if expected_cells != observed_cells or len(rows) != len(expected_cells):
        raise ValueError("Incomplete or duplicated release replay")
    args.output.parent.mkdir(exist_ok=True,parents=True)
    with args.output.open("w",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=rows[0],lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    if any(not r["identical"] or r["modes"] == 0 for r in rows):
        raise ValueError("Release replay differs from archived outputs")
    print(f"{sum(r['identical'] for r in rows)}/{len(rows)} matching files, {sum(r['modes'] for r in rows)} modes")


if __name__ == "__main__":
    main()
