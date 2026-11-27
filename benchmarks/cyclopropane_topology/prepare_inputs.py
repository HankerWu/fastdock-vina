"""Prepare the four registered cocrystals using the common independent-query recipe."""

import os
import sys
import json
import numpy as np
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "search_backtracking/pose"))
import backtracking_panel as panel

CASES = {
    "1XMU": ("PDE4B", "ROF", "A"),
    "3WZD": ("VEGFR2", "LEV", "A"),
    "7AAD": ("PARP1", "09L", "A"),
    "7XRR": ("OX2R", "NRK", "A"),
}

if __name__ == "__main__":
    os.environ.setdefault("OPENMM_CPU_THREADS", "1")
    panel.CASES = CASES
    for pdb in sys.argv[2:] or CASES:
        work = Path(sys.argv[1]).resolve()
        started = perf_counter()
        panel.prepare_input(work, pdb)
        folder = work / pdb
        meta = json.loads((folder / "preparation.json").read_text())
        if meta["receptor_status"].startswith("omitted nonwater"):
            lines = [
                r
                for r in (work / "raw" / (pdb + ".pdb")).read_text().splitlines()
                if r.startswith(("ATOM  ", "HETATM")) and r[21] == CASES[pdb][2]
            ]
            keys = dict.fromkeys(r[22:27] for r in lines if r.startswith("ATOM  "))
            ordinal = {key: i for i, key in enumerate(keys, 1)}
            rec = [
                r[:22] + f"{ordinal[r[22:27]]:4d} " + r[27:]
                if r.startswith("ATOM  ")
                else r
                for r in lines
            ]
            meta.update(
                panel.prepare.base.receptor(
                    rec, folder, pdb, np.array(meta["mapped"][0])
                )
            )
            (folder / "preparation.json").write_text(json.dumps(meta, indent=2) + "\n")
        source = json.loads((folder / "source.json").read_text())
        source["seconds_including_final_receptor"] = perf_counter() - started
        (folder / "source.json").write_text(json.dumps(source, indent=2) + "\n")
