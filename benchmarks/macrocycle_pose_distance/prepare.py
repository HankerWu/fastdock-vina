"""Prepare the two previously unused cognate structures with the frozen protocol."""
import argparse
import os
from pathlib import Path
import sys
from urllib.request import urlretrieve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "macrocycle_proposal" / "external"))
import prepare_panel

CASES = {"4EGK": ("HSP90a", "RDC", "A"), "2J9M": ("CDK2", "PY8", "A")}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    prepare_panel.CASES.update(CASES)
    for pdb, (_, component, chain) in CASES.items():
        prepare_panel.run(HERE.parents[1], work, pdb)
        lines = (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
        residue = next(l[22:26].strip() for l in lines
                       if l.startswith("HETATM") and l[17:20] == component and l[21] == chain)
        urlretrieve(f"https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id={chain}"
                    f"&auth_seq_id={residue}&encoding=sdf&copy_all_categories=false",
                    work / "raw" / (pdb + "-bound.sdf"))
