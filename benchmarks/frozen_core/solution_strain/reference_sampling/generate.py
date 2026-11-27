"""Generate exactly one independent eight-conformer reference library per ligand."""

import argparse
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from solvent_strain import (
    free_conformers,
    lectins,
    ligand_context,
    physical,
    proteases,
    write,
)

SEED = 1618033


def generate(out, seed=SEED):
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    with TemporaryDirectory() as temp:
        for panel in (proteases, lectins):
            for meta, initial, _ in panel(Path(temp)):
                ligand, _ = physical.decode_ligand(initial)
                context = ligand_context(ligand, "Reference")
                library = free_conformers(ligand, context, seed=seed)
                if not library:
                    raise ValueError(f"No converged conformer for {meta['compound']}")
                rows.extend(dict(compound=meta["compound"], **r) for r in library)
    if len({str(r["compound"]) for r in rows}) != 55:
        raise ValueError("Incomplete ligand panel")
    write(out / "independent_references.csv.xz", rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    generate(parser.parse_args().out)
