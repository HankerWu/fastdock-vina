"""Control for nonspecific pocket localization by anchoring the distal acetimidoylpiperidine."""

import argparse
import itertools
from pathlib import Path
from rdkit import Chem
from experiment import ROOT, SEEDS, aspsite, run, table
from make_ion_pair_anchors import GAP
from score_physical_pool import serial_index_pairs


def inverted_anchor(f, oxy):
    lines = f.read_text().splitlines()
    smiles = next(
        l.split()[2]
        for l in lines
        if l.startswith("REMARK SMILES ") and not l.startswith("REMARK SMILES IDX")
    )
    mapping = {i - 1: s for i, s in serial_index_pairs(lines, "REMARK SMILES IDX ")}
    (match,) = Chem.MolFromSmiles(smiles).GetSubstructMatches(
        Chem.MolFromSmarts("[CH3]-C(=[N+;H2])-[N;H0;R]")
    )
    ns = [mapping[n] for n in match[2:]]
    blocks = [
        "\n".join(
            f"{n} " + " ".join(f"{x:.8f}" for x in o) + f" shell {GAP}"
            for n, o in zip(ns, os)
        )
        for os in (oxy, oxy[::-1])
    ]
    return "\nOR\n".join(blocks) + "\n"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ("out", "fjs", "p16", "binary"):
        p.add_argument(k, type=Path)
    p.add_argument("--jobs", type=int, default=2)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    compounds = [
        r["ligand"]
        for r in table(ROOT / "benchmarks/affinity_partition/fxa_series.csv")
        if int(r["included"])
    ]
    sources = {}
    for rec, folder in (("1FJS", a.fjs), ("2P16", a.p16)):
        oxy = aspsite(folder / f"{rec}.pdb", "A")[1:]
        for lig in compounds:
            f = folder / f"{lig}.pdbqt"
            (a.out / f"{rec}_{lig}.inverted").write_text(inverted_anchor(f, oxy))
            sources[rec, lig] = f, folder / f"{rec}.pdbqt", folder / "box.txt"
    run(
        a.binary.resolve(),
        a.out,
        sources,
        itertools.product(("1FJS", "2P16"), compounds, SEEDS, ("inverted",)),
        a.jobs,
    )


if __name__ == "__main__":
    main()
