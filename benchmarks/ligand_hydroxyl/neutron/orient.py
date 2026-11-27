"""Predict ligand alcohol orientations from heavy atoms and fixed external donors.

Ligand D observations are isolated from prediction in read_structure.  This
small controlled calculation uses the deposited asymmetric unit.
"""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import softmax

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "frozen_core/hydroxyl_partition"))
from scorer import RT, quality

sys.path.insert(0, str(HERE.parent))
from acceptor import lp_gate

CHECKSUM = {
    "6EXY": "6aa696bfc840b63085899c6762505d74873a91fff7c15a122fee8830169cb702",
    "6EYM": "ec293ca0766a0a552fe3e55b0a0ad76941d3bbd2a4cfce4344d2ea358999b03b",
}
LIGANDS = {"GOL", "BGC", "GAL"}


@dataclass(eq=False)
class Atom:
    idx: int
    group: str
    site: str
    name: str
    element: str
    xyz: np.ndarray
    occupancy: float
    b_factor: float

    @property
    def label(self):
        return f"{self.group}:{self.site}:{self.name}"


def read_structure(pdb):
    data = gzip.decompress((HERE / f"{pdb}.pdb.gz").read_bytes())
    if hashlib.sha256(data).hexdigest() != CHECKSUM[pdb]:
        raise ValueError(f"Coordinate checksum mismatch for {pdb}")
    ligand, other, observed = {}, [], {}
    graph = defaultdict(set)
    for line in data.decode().splitlines():
        if line.startswith("CONECT"):
            ids = [
                int(line[i : i + 5]) for i in range(6, 31, 5) if line[i : i + 5].strip()
            ]
            for j in ids[1:]:
                graph[ids[0]].add(j)
                graph[j].add(ids[0])
        if not line.startswith(("ATOM  ", "HETATM")) or line[16] not in " A":
            continue
        a = Atom(
            int(line[6:11]),
            line[17:20].strip(),
            line[21:27].strip(),
            line[12:16].strip(),
            line[76:78].strip(),
            np.array([float(line[j : j + 8]) for j in (30, 38, 46)]),
            float(line[54:60]),
            float(line[60:66]),
        )
        if a.group not in LIGANDS:
            other.append(a)
        elif a.element not in {"H", "D"}:
            ligand[a.idx] = a
        elif a.element == "D" and a.name.startswith("DO"):
            observed[a.group, a.site, a.name[1:]] = a
    heavy_graph = {i: graph[i] & ligand.keys() for i in ligand}
    return ligand, heavy_graph, other, observed


def unit(x):
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def circle(oxygen, carbon, frame, nodes):
    axis = unit(oxygen - carbon)
    ref = frame - carbon
    u = unit(ref - (ref @ axis) * axis)
    v = np.cross(axis, u)
    phi = np.arange(nodes) * 2 * np.pi / nodes
    oh = axis / 3 + np.sqrt(8 / 9) * (
        np.cos(phi)[:, None] * u + np.sin(phi)[:, None] * v
    )
    return oxygen + 0.97 * oh, phi, u, v


def environments(oxygen, other, ligand, graph):
    near = [a for a in other if np.linalg.norm(a.xyz - oxygen.xyz) < 5.0]
    groups = defaultdict(list)
    for a in near:
        if a.element in {"H", "D"}:
            groups[a.group, a.site].append(a)
    for a in near:
        if a.element not in {"N", "O"}:
            continue
        hs = [
            h for h in groups[a.group, a.site] if np.linalg.norm(h.xyz - a.xyz) < 1.28
        ]
        # Deposited D occupancy resolves histidine tautomer for this diagnostic.
        acceptor = a.element == "O" or (
            a.group == "HIS" and a.name in {"ND1", "NE2"} and not hs
        )
        yield "water" if a.group == "HOH" else "protein", a, acceptor, hs
    shell = {oxygen.idx}
    for _ in range(3):
        shell |= {j for i in shell for j in graph[i]}
    for i, a in ligand.items():
        if i not in shell and a.element == "O":
            yield "within_ligand", a, True, []


def evaluate_site(oxygen, carbon, frame, other, ligand, graph, nodes):
    h, phi, u, v = circle(oxygen.xyz, carbon.xyz, frame.xyz, nodes)
    # Declare primary arms; use each donor-acceptor pair once.
    names = ("flat", "donor", "bidirectional", "lone_pairs")
    sums = {
        (env, name): np.zeros(nodes)
        for env in ("dry", "deposited_water", "intraligand", "water_and_intraligand")
        for name in names
    }
    contacts = []
    for kind, partner, acceptor, hs in environments(oxygen, other, ligand, graph):
        r = np.linalg.norm(partner.xyz - oxygen.xyz)
        radius = 1.8 if partner.element == "N" else 1.7
        a = 0.587439 * np.clip((1.7 + radius - r) / 0.7, 0, 1)
        if not a or not (acceptor or hs):
            continue
        forward = quality(oxygen.xyz, h, partner.xyz) if acceptor else np.zeros(nodes)
        reverse = np.zeros(nodes)
        lp_reverse = np.zeros(nodes)
        for d in hs:
            donation = quality(partner.xyz, d.xyz, oxygen.xyz)
            reverse = np.maximum(reverse, donation)
            lp_reverse = np.maximum(
                lp_reverse, donation * lp_gate(oxygen.xyz, carbon.xyz, h, d.xyz)
            )
        terms = {
            "donor": forward,
            "bidirectional": np.maximum(forward, reverse),
            "lone_pairs": np.maximum(forward, lp_reverse),
        }
        arms = {
            "protein": (
                "dry",
                "deposited_water",
                "intraligand",
                "water_and_intraligand",
            ),
            "water": ("deposited_water", "water_and_intraligand"),
            "within_ligand": ("intraligand", "water_and_intraligand"),
        }
        for env in arms[kind]:
            for name, q in terms.items():
                sums[env, name] -= a * q
        contacts.append(
            dict(
                kind=kind,
                partner=partner.label,
                distance_A=r,
                occupancy=partner.occupancy,
                accepts=int(acceptor),
                donor_H=" ".join(d.name for d in hs),
                reverse_quality=reverse[0],
                radial_magnitude=a,
            )
        )
    return phi, u, v, sums, contacts


def stats(phi, energy, reference):
    diff = np.abs(np.angle(np.exp(1j * (phi - reference))))
    near = diff <= np.pi / 6 + 1.0e-12
    minima = np.isclose(energy, energy.min(), atol=1.0e-9, rtol=0)
    informative = not np.all(minima)
    return dict(
        span=float(np.ptp(energy)),
        posterior_pm30=float(softmax(-energy / RT)[near].sum()),
        uniform_pm30=float(near.mean()),
        minimum_deg=float(np.degrees(diff[minima].min())) if informative else "",
        maximum_tied_deg=float(np.degrees(diff[minima].max())) if informative else "",
    )


def run(nodes):
    rows, contacts = [], []
    for pdb in CHECKSUM:
        ligand, graph, other, observed = read_structure(pdb)
        sites = [
            o for o in ligand.values() if o.element == "O" and len(graph[o.idx]) == 1
        ]
        if {(o.group, o.site, o.name) for o in sites} != set(observed):
            raise ValueError("Missing D observation, or ambiguous hydroxyl assignment")
        for oxygen in sites:
            carbon = ligand[next(iter(graph[oxygen.idx]))]
            frame = ligand[min(graph[carbon.idx] - {oxygen.idx})]
            phi, u, v, energy, near = evaluate_site(
                oxygen, carbon, frame, other, ligand, graph, nodes
            )
            # Observations enter only here, after computing every angular field.
            d = observed[oxygen.group, oxygen.site, oxygen.name]
            od = d.xyz - oxygen.xyz
            observed_phi = np.arctan2(od @ v, od @ u)
            common = dict(pdb=pdb, site=oxygen.label)
            contacts.extend({**common, **item} for item in near)
            for (env, arm), e in energy.items():
                rows.append(
                    dict(
                        **common,
                        env=env,
                        arm=arm,
                        d_occupancy=d.occupancy,
                        d_B=d.b_factor,
                        OD_length_A=np.linalg.norm(od),
                        observed_deg=np.degrees(observed_phi),
                        **stats(phi, e, observed_phi),
                    )
                )
    return rows, contacts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--nodes", type=int, default=720)
    args = parser.parse_args()
    if args.nodes < 3:
        parser.error("At least three nodes required")
    args.output.mkdir(parents=True, exist_ok=True)
    for name, rows in zip(("sites", "contacts"), run(args.nodes), strict=True):
        with (args.output / (name + ".csv")).open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
