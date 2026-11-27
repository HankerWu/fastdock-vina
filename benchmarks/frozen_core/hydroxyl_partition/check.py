"""Separate scalar radial-angle construction and fine quadrature for two contrasting rank reversals."""
import argparse
import csv
import math
from pathlib import Path
import tempfile

import numpy as np
from scipy.spatial.transform import Rotation

from evaluate import archive_cases, blocks, write
from scorer import AXES, Receptor, atoms, ligand_sites, score_pose

CASES = {
    ("FABP4", "3FR2", "2NNQ"): (1, 3),
    ("FKBP1A", "1FKG", "1FKH"): (1, 2),
}
RT = 1.98720425864083e-3 * 298.15


def radial(d, a):
    radius = {"O": 1.7, "N": 1.8}
    distance = math.dist(d.xyz, a.xyz)
    return .587439 * min(1., max(0., (radius[d.kind[0]] + radius[a.kind[0]] - distance) / .7))


def contact_factor(oh, da):
    # dot(OH, HA) = dot(OH, OA)-|OH|^2. Scalar atom distances supply |HA|.
    length2 = (oh*oh).sum(axis=-1)
    dot = oh @ da
    ha2 = np.maximum(0, da@da + length2 - 2*dot)
    return np.maximum(0, dot-length2)**2 / np.maximum(1e-30, length2 * ha2)


def pmf(prior, cost):
    weights = np.exp(-(prior-prior.min())/RT)
    return -RT * math.log(np.dot(weights,np.exp(-cost/RT))/weights.sum())


def detailed(rec, lig_text, block, n=4096):
    ligand, owners = ligand_sites(lig_text)
    pose = atoms(block)
    rec_atoms = rec.atoms
    donors = rec.hydrogens
    active, contacts, sites = {}, [], []
    nonrotor = 0.
    for i, li in pose.items():
        if li.kind not in ("O","OA","N","NA"):
            continue
        for j,rj in rec_atoms.items():
            if rj.kind not in ("O","OA","N","NA"):
                continue
            e = radial(rj, li)
            if not e:
                continue
            hs1 = owners[i] if rj.kind in ("OA", "NA") else []
            hs2 = donors[j] if li.kind in ("OA", "NA") else []
            if not(hs1 or hs2):
                continue
            g1 = max((contact_factor(pose[h].xyz-li.xyz,rj.xyz-li.xyz) for h in hs1),default=0.)
            g2 = max((contact_factor(rec_atoms[h].xyz-rj.xyz,li.xyz-rj.xyz) for h in hs2),default=0.)
            name = AXES.get((rj.residue[:3],rj.name)) if hs2 else None
            fixed = e*(1-max(g1,g2))
            contacts.append(dict(ligand_atom=i,ligand_type=li.kind,receptor_residue=rj.residue.strip(),
                                 receptor_name=rj.name,receptor_serial=j,radial=-e,donor_ligand=g1,donor_protein=g2,
                                 fixed=fixed,moving=int(bool(name))))
            if name:
                active.setdefault(j, []).append((e, g1, li))
            else:
                nonrotor += fixed
    for j, terms in active.items():
        o = rec_atoms[j]
        (h,) = donors[j]
        h0 = rec_atoms[h].xyz-o.xyz
        (parent,) = [a for a in rec_atoms.values()
                     if a.residue==o.residue and a.name==AXES[o.residue[:3],o.name]]
        axis = (o.xyz-parent.xyz) / np.linalg.norm(o.xyz-parent.xyz)
        # Independent general SO(3) rotation operator, not the production circle formula.
        angles = np.arange(n)*2*math.pi/n
        oh = Rotation.from_rotvec(angles[:,None]*axis).apply(h0)
        prior = np.zeros(n)
        for a in rec_atoms.values():
            if a.residue==o.residue or a.kind not in ("NA","OA"):
                continue
            e = radial(o,a)
            if e:
                prior -= e * contact_factor(oh, a.xyz-o.xyz)
        single_costs = [e*(1-np.maximum(g1,contact_factor(oh, li.xyz-o.xyz))) for e,g1,li in terms]
        combined = sum(single_costs)
        ensemble = pmf(prior,combined)
        minimum = combined.min()
        sites.append(dict(oxygen=o.residue.strip()+" "+o.name,acceptors=len(terms),
                          independent_pair_minimum=sum(k.min() for k in single_costs),
                          joint_minimum=minimum,independent_pair_marginal=sum(pmf(prior,k) for k in single_costs),
                          joint_marginal=ensemble,entropy_and_prior_cost=ensemble-minimum,
                          fraction_within_RT_of_minimum=float(np.mean(combined<=minimum+RT)),
                          prior_range=float(np.ptp(prior))))
    predicted = score_pose(rec, lig_text, block)
    assert abs(sum(c["fixed"] for c in contacts)-predicted["raw_fixed"]) < 1e-10
    assert abs(nonrotor + sum(r["joint_marginal"] for r in sites) - predicted["raw_ensemble"]) < .0001
    return contacts,sites


def toy_checks():
    # Exact two-state examples: unused, constant penalty, and a mirror pair.
    flat = np.zeros(2)
    assert abs(pmf(np.array([-11.,3.]), np.zeros(2))) < 1e-12
    assert abs(pmf(np.array([5.,-3.]), np.full(2, .3))-.3) < 1e-12
    costs = np.array([0.,.587439])
    assert abs(pmf(flat, costs+costs[::-1])-.587439)<1e-12
    assert pmf(flat,costs+costs[::-1]) > 2*pmf(flat,costs)
    assert abs(pmf(flat,costs) + RT*math.log((1+math.exp(-.587439/RT))/2)) < 1e-12


def main(out):
    toy_checks()
    contacts, sites_out = [], []
    with tempfile.TemporaryDirectory() as temp:
        for case, sites, poses, table in archive_cases(Path(temp)):
            key = case["target"],case["query"],case["receiver"]
            if key not in CASES or case["design"]!="per_complex" or case["seed"]!="93529":
                continue
            receptor = Receptor((sites/case["receiver"]/"receptor.pdbqt").read_text())
            ligand = (sites/case["query"]/"ligand.pdbqt").read_text()
            for rank in CASES[key]:
                base=dict(target=case["target"],query=case["query"],receiver=case["receiver"],seed=case["seed"],rank=rank)
                c,s=detailed(receptor,ligand,blocks(poses)[rank-1])
                contacts.extend({**base,**row} for row in c)
                sites_out.extend({**base,**row} for row in s)
    assert len(sites_out)==4
    write(out/"case_contacts.csv",contacts)
    write(out/"case_rotors.csv",sites_out)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    main(parser.parse_args().out)
