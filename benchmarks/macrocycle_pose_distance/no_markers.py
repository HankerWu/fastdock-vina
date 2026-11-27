"""Check byte-identical poses and work counts for rigid and flexible marker-free inputs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from benchmarks.concerted_ring.run import counts, sha

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("parent","fixed","out"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    a.out.mkdir(exist_ok=True, parents=True)
    sample = ROOT / "examples" / "basic"
    cases = {"1iep":["--receptor",sample/"1iep_receptor.pdbqt"],
             "flex_fixture":["--receptor",ROOT/"benchmarks/fixtures/flex-rigid.pdbqt",
                             "--flex",ROOT/"benchmarks/fixtures/flex-residue.pdbqt"]}
    rows = {}
    for key,receptor in cases.items():
        result = {}
        for arm in ("parent","fixed"):
            dest = a.out / f"{key}_{arm}.pdbqt"
            cmd = [getattr(a,arm).resolve(), *receptor, "--ligand", sample/"1iep_ligand.pdbqt",
                   "--config", sample/"1iep_receptor.box.txt", "--num_modes","5","--cpu","1",
                   "--exhaustiveness","1","--seed","191009","--max_evals","3000","--energy_range","100",
                   "--verbosity","0","--search_counters","--out",dest]
            log = subprocess.check_output(cmd,stderr=subprocess.STDOUT,text=True)
            result[arm] = dict(hash=sha(dest), **counts(log))
        assert result["parent"] == result["fixed"]
        rows[key] = result["fixed"]
    (a.out / "parity.json").write_text(json.dumps(rows, indent=2)+"\n")
