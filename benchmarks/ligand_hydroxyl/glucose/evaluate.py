"""Compare bidirectional alcohol score arms on shared glucose pools."""

import argparse
import csv
import json
from pathlib import Path
import runpy
import sys
import time

import numpy as np
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from acceptor import lp_gate
from model import (
    ACCEPTORS,
    Receptor,
    atoms,
    hydrogenated_positions,
    magnitude,
    marginal,
    quality,
    ring,
    score,
    topology,
    torsional_prior,
)

count = runpy.run_path(
    str(HERE.parents[1] / "terminal_rotor_score/independent/verify.py")
)["independently_count"]


def alcohol_acceptance(rec, query, block):
    """Additional penalties relative to the previously frozen donor marginal."""
    heavy, serial, sites = topology(query)
    xyz = atoms(block)
    all_xyz = None
    delta_scope = delta_geom = 0.0
    active = reverse_only = 0
    inverse = {s: j for j, s in serial.items()}
    for sid, site in sites.items():
        o = xyz[sid]
        if (
            o.kind not in ACCEPTORS
            or heavy.GetAtomWithIdx(inverse[site.carbon]).GetHybridization()
            != Chem.HybridizationType.SP3
        ):
            continue
        contacts = []
        for ix in rec.tree.query_ball_point(o.xyz, 3.6):
            j = rec.ids[ix]
            partner = rec.atoms[j]
            weight = magnitude(o, partner)
            hs = rec.hydrogens[j]
            if weight and (hs or partner.kind in ACCEPTORS):
                contacts.append((partner, weight, [rec.atoms[h].xyz for h in hs]))
        if not any(hs for _, _, hs in contacts):
            continue
        circle, around, axis = ring(site, xyz, 72)
        carbon = xyz[site.carbon].xyz
        old = np.zeros(len(circle))
        new = np.zeros(len(circle))
        for a, w, hs in contacts:
            forward = quality(o.xyz, circle, a.xyz) if a.kind in ACCEPTORS else 0.0
            reverse = [quality(a.xyz, h, o.xyz) for h in hs]
            rev_plain = max(reverse, default=0.0)
            rev_lp = np.zeros(len(circle))
            for g, h in zip(reverse, hs):
                rev_lp = np.maximum(rev_lp, g * lp_gate(o.xyz, carbon, circle, h))
            old += w * (1 - np.maximum(forward, rev_plain))
            new += w * (1 - np.maximum(forward, rev_lp))
            if a.kind not in ACCEPTORS:
                reverse_only += 1
                delta_scope += w * (1 - rev_plain)
        if any(i not in serial for i, _ in site.terms) and all_xyz is None:
            all_xyz = hydrogenated_positions(heavy, serial, xyz)
        prior = torsional_prior(site, around, axis, xyz, serial, all_xyz)
        delta_geom += marginal(new, prior) - marginal(old, prior)
        active += 1
    return delta_scope, delta_geom, active, reverse_only


def write(path, rows):
    with path.open("w", newline="") as stream:
        out = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def evaluate(w, out):
    out.mkdir(exist_ok=True, parents=True)
    meta = json.loads((w / "metadata.json").read_text())
    recs = {}
    poses, cells = [], []
    for p in sorted((w / "runs").glob("*.pdbqt")):
        case, cap, seed, engine = p.stem.split("_")
        m = meta[case]
        root = w / m["pdb"]
        query = (root / (m["comp"] + ".pdbqt")).read_text()
        scale = 1 / (1 + 0.05846 * count(query.encode()))
        if engine == "pro" and root not in recs:
            recs[root] = Receptor((root / "receptor.pdbqt").read_text())
        rows = []
        common = dict(case=case, cap=cap, seed=seed, engine=engine)
        for rank, block in enumerate(
            (b for b in p.read_text().split("ENDMDL") if b.strip()), 1
        ):
            xyz = atoms(block)
            ref = np.asarray(m["mapped"])
            err = ref - np.array([xyz[i].xyz for i in m["serial"]])
            rmsd = float(np.sqrt(np.square(err).sum(axis=2).mean(axis=1)).min())
            native = float(
                next(
                    l for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")
                ).split()[3]
            )
            row = dict(**common, rank=rank, rmsd=rmsd, scale=scale, native=native)
            fields = (
                "donor",
                "bidirectional",
                "lone_pairs",
                "scope_delta_raw",
                "geometry_delta_raw",
                "sites",
                "reverse_only",
                "donor_s",
                "extension_s",
            )
            row.update(dict.fromkeys(fields, ""))
            if engine == "pro":
                start = time.perf_counter()
                base = score(recs[root], query, block)["bonded"]
                intermediate = time.perf_counter()
                scope, lp, n, rev = alcohol_acceptance(recs[root], query, block)
                end = time.perf_counter()
                row.update(
                    donor=native + scale * base,
                    bidirectional=native + scale * (base + scope),
                    lone_pairs=native + scale * (base + scope + lp),
                    scope_delta_raw=scope,
                    geometry_delta_raw=lp,
                    sites=n,
                    reverse_only=rev,
                    donor_s=intermediate - start,
                    extension_s=end - intermediate,
                )
            rows.append(row)
        if not rows:
            raise ValueError(f"Empty pose file {p}")
        poses.extend(rows)
        for arm in (
            ("native", "donor", "bidirectional", "lone_pairs")
            if engine == "pro"
            else ("native",)
        ):
            order = sorted(rows, key=lambda a: (a[arm], a["rank"]))
            cells.append(
                dict(
                    **common,
                    arm=arm,
                    n=len(rows),
                    winner=order[0]["rank"],
                    rmsd1=order[0]["rmsd"],
                    rmsd5=min(a["rmsd"] for a in order[:5]),
                    oracle=min(a["rmsd"] for a in rows),
                )
            )
    write(out / "poses.csv", poses)
    write(out / "cells.csv", cells)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("work", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    evaluate(a.work, a.out)
