"""Post hoc placement control. Crystal-derived probes never enter reported predictions."""

import argparse
import csv
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from experiment import CASES, HERE, cross, timed_process, write_csv
import runpy

count_torsions = runpy.run_path(
    str(HERE.parents[1] / "terminal_rotor_score/independent/verify.py")
)["independently_count"]

sys.path.insert(0, str(HERE.parent))
from model import CHPi, WEIGHT


def crystal_template(meta, source_text, points):
    # Keep the *input* neutral atom order when assigning protonated coordinates.
    base = Chem.MolFromSmiles(meta["canonical_smiles"])
    (n,) = [a for a in base.GetAtoms() if a.GetAtomicNum() == 7]
    n.SetNumExplicitHs(n.GetTotalNumHs() + 1)
    n.SetFormalCharge(+1)
    Chem.SanitizeMol(base)
    conf = Chem.Conformer(base.GetNumAtoms())
    for i, point in enumerate(points):
        conf.SetAtomPosition(i, point)
    base.AddConformer(conf)
    m = Chem.AddHs(base, addCoords=True)
    force = AllChem.MMFFGetMoleculeForceField(m, AllChem.MMFFGetMoleculeProperties(m))
    for i in range(base.GetNumAtoms()):
        force.AddFixedPoint(i)
    force.Initialize()
    if force.Minimize(maxIts=2000):
        raise ValueError("Constrained hydrogen relaxation did not converge")
    xyz = m.GetConformer().GetPositions()
    assert np.max(np.abs(xyz[: len(points)] - points)) < 1e-7
    idx = {}
    for line in source_text.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            vals = list(map(int, line.split()[3:]))
            idx.update(zip(vals[1::2], (v - 1 for v in vals[::2])))
    out = []
    for line in source_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            pos = xyz[idx[int(line[6:11])]]
            line = line[:30] + "".join(f"{v:8.3f}" for v in pos) + line[54:]
        out.append(line)
    return "\n".join(out) + "\n"


def analyze(w, binary, out):
    meta = json.loads((w / "metadata.json").read_text())
    out.mkdir(exist_ok=True, parents=True)
    probe = w / "probes"
    probe.mkdir(exist_ok=True)
    rows = []
    for src, dst in itertools.product(CASES, repeat=2):
        template = (w / src / "protonated.pdbqt").read_text()
        model = CHPi((w / dst / "receptor.pdbqt").read_text(), template)
        scale = 1 / (1 + 0.05846 * count_torsions(template.encode()))
        rotation, shift, _, _ = cross.alignment(meta[src], meta[dst])
        refs = np.asarray(meta[src]["mapped"]) @ rotation + shift
        start = probe / f"{src}_{dst}_deposited.pdbqt"
        relaxed = probe / f"{src}_{dst}_relaxed.pdbqt"
        start.write_text(crystal_template(meta[src], template, refs[0]))
        args = [
            str(binary),
            "--cpu",
            "1",
            "--seed",
            "43103",
            "--receptor",
            str(w / dst / "receptor.pdbqt"),
        ]
        for a, v in zip("xyz", meta[dst]["center"]):
            args.extend(["--center_" + a, str(v), "--size_" + a, "24"])
        timed_process(
            [*args, "--ligand", str(start), "--local_only", "--out", str(relaxed)],
            w,
            relaxed.with_suffix(".log"),
        )
        for stage, p in (("deposited", start), ("relaxed", relaxed)):
            if stage == "relaxed":
                if not p.is_file():
                    rows.append(
                        dict(
                            query=src,
                            receptor=dst,
                            stage=stage,
                            seed="",
                            status="missing",
                            rmsd="",
                            inter="",
                            intra="",
                            chpi_raw="",
                            gap_native="",
                            gap_chpi="",
                        )
                    )
                    continue
            body = p.read_text()
            xyz = {
                int(l[6:11]): cross.coords(l)
                for l in body.splitlines()
                if l.startswith(("ATOM  ", "HETATM"))
            }
            err = np.sqrt(
                (
                    (
                        refs
                        - np.asarray([xyz[i] for i in meta[src]["protonated_serial"]])
                    )
                    ** 2
                )
                .sum(axis=2)
                .mean(axis=1)
            ).min()
            single = (
                "\n".join(
                    l
                    for l in body.splitlines()
                    if not l.startswith(("MODEL", "ENDMDL"))
                )
                + "\n"
            )
            clean = probe / f"{src}_{dst}_{stage}_single.pdbqt"
            clean.write_text(single)
            log = clean.with_suffix(".log")
            timed_process([*args, "--ligand", str(clean), "--score_only"], w, log)
            report = log.read_text()
            inter = float(
                re.search(r"Final Intermolecular Energy\s*:\s*([-0-9.]+)", report)[1]
            )
            intra = float(
                re.search(r"Final Total Internal Energy\s*:\s*([-0-9.]+)", report)[1]
            )
            corr = sum(model.components(single))
            for seed in cross.SEEDS:
                blocks = (
                    (w / "runs" / f"{src}_{dst}_protonated_pro_{seed}.pdbqt")
                    .read_text()
                    .split("ENDMDL")
                )
                raw = []
                adj = []
                for block in (b for b in blocks if b.strip()):
                    value = float(
                        next(
                            l.split()[-1]
                            for l in block.splitlines()
                            if l.startswith("REMARK INTER + INTRA:")
                        )
                    )
                    raw.append(value)
                    adj.append(value + WEIGHT * sum(model.components(block)))
                rows.append(
                    dict(
                        query=src,
                        receptor=dst,
                        stage=stage,
                        seed=seed,
                        status="ok",
                        rmsd=err,
                        inter=inter,
                        intra=intra,
                        chpi_raw=corr,
                        gap_native=scale * (inter + intra - min(raw)),
                        gap_chpi=scale * (inter + intra + WEIGHT * corr - min(adj)),
                    )
                )
    write_csv(out / "crystal_diagnostic.csv", rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("binary", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    analyze(a.work.resolve(), a.binary.resolve(), a.output.resolve())
