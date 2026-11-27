"""Score a ligand's alcohol and phenol donor orientation on frozen Vina poses."""
import argparse
import csv
from functools import lru_cache
from pathlib import Path
import sys
from typing import NamedTuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH / "frozen_core/hydroxyl_partition"))
from scorer import ACCEPTORS, RADII, Receptor, atoms, ligand_sites, magnitude, marginal, quality


class Hydroxyl(NamedTuple):
    oxygen: int              # PDBQT serials
    carbon: int
    hydrogen: int
    reference: int           # heavy neighbor of C that sets a common azimuth
    terms: tuple             # (all-H graph atom index for X, MMFF V1,V2,V3)


@lru_cache(maxsize=128)
def topology(ligand_text):
    parsed, hydrogens = ligand_sites(ligand_text)
    smiles = next(line[14:] for line in ligand_text.splitlines()
                  if line.startswith("REMARK SMILES ") and not line.startswith("REMARK SMILES IDX "))
    heavy = Chem.MolFromSmiles(smiles)
    if heavy is None:
        raise ValueError("Invalid Meeko SMILES graph")
    serial = {}
    for line in ligand_text.splitlines():
        if line.startswith("REMARK SMILES IDX "):
            values = list(map(int, line.split()[3:]))
            serial.update((i-1, s) for i, s in zip(values[::2], values[1::2], strict=True))
    if set(serial) != set(range(heavy.GetNumAtoms())) or len(set(serial.values())) != len(serial):
        raise ValueError("Meeko SMILES IDX must map the whole heavy graph bijectively")
    full = Chem.AddHs(heavy)
    props = AllChem.MMFFGetMoleculeProperties(full, mmffVariant="MMFF94s")
    sites = {}
    for oxygen in heavy.GetAtoms():
        if (oxygen.GetAtomicNum() != 8 or oxygen.GetFormalCharge() != 0 or
                oxygen.GetDegree() != 1 or oxygen.GetTotalNumHs() != 1):
            continue
        c = oxygen.GetNeighbors()[0]
        if c.GetAtomicNum() != 6 or not (c.GetIsAromatic() or c.GetHybridization() == Chem.HybridizationType.SP3):
            continue
        o, center = oxygen.GetIdx(), c.GetIdx()
        reference = sorted(n.GetIdx() for n in c.GetNeighbors() if n.GetIdx() != o)
        if not reference:
            continue
        hs = hydrogens[serial[o]]
        if len(hs) != 1 or parsed[serial[o]].kind not in {"O", "OA"}:
            raise ValueError("SMILES hydroxyl disagrees with PDBQT oxygen or H PARENT")
        h = next(n.GetIdx() for n in full.GetAtomWithIdx(o).GetNeighbors() if n.GetAtomicNum() == 1)
        torsions = []
        for x in full.GetAtomWithIdx(center).GetNeighbors():
            i = x.GetIdx()
            if i == o:
                continue
            p = props.GetMMFFTorsionParams(full, i, center, o, h) if props is not None else None
            torsions.append((i, None if p is None else p[1:]))
        sites[serial[o]] = Hydroxyl(serial[o], serial[center], hs[0], serial[reference[0]], tuple(torsions))
    return heavy, serial, sites


def unit(vec):
    norm = np.linalg.norm(vec)
    if norm < 1.e-9:
        raise ValueError("Degenerate hydroxyl frame")
    return vec / norm


def ring(site, xyz, nodes, phase=0.):
    oxygen, carbon, hydrogen, ref = (xyz[i].xyz for i in site[:4])
    axis = unit(oxygen - carbon)
    oh = hydrogen - oxygen
    axial = np.dot(oh, axis) * axis
    radius = np.linalg.norm(oh-axial)
    refvec = ref-carbon
    u = unit(refvec - np.dot(refvec, axis)*axis)
    v = np.cross(axis, u)
    phi = (np.arange(nodes) + phase)*(2*np.pi/nodes)
    around = np.cos(phi)[:,None] * u + np.sin(phi)[:,None] * v
    return oxygen + axial + radius*around, around, axis


def hydrogenated_positions(heavy, serial, xyz):
    """RDKit reconstructs only the unrepresented aliphatic C–H reference directions."""
    mol = Chem.Mol(heavy)
    conf = Chem.Conformer(heavy.GetNumAtoms())
    for i, s in serial.items():
        conf.SetAtomPosition(i, tuple(xyz[s].xyz))
    mol.AddConformer(conf)
    return Chem.AddHs(mol, addCoords=True).GetConformer().GetPositions()


def torsional_prior(site, around, axis, xyz, serial, full_xyz):
    carbon = xyz[site.carbon].xyz
    prior = np.zeros(len(around))
    for i, params in site.terms:
        if params is None:
            raise ValueError(f"MMFF94s torsion unavailable at ligand O {site.oxygen}")
        x = xyz[serial[i]].xyz if i in serial else full_xyz[i]
        arm = x - carbon
        cos = np.clip(around @ unit(arm - np.dot(arm, axis)*axis), -1., 1.)
        v1, v2, v3 = params
        prior += .5*(v1*(1+cos) + 2*v2*(1-cos*cos) + v3*(1+4*cos**3-3*cos))
    return prior


ARMS = ("fixed", "turned", "minimum", "uniform", "bonded",
        "full_fixed", "full_turned", "full_uniform", "full_bonded")


def score(rec, ligand_text, block, n=72, phase=0.):
    if n < 3:
        raise ValueError("Provide at least three angular nodes")
    parsed, owners = ligand_sites(ligand_text)
    mol, serial, sites = topology(ligand_text)
    xyz = atoms(block)
    if xyz.keys() != parsed.keys() or any(a.kind != parsed[i].kind for i,a in xyz.items()):
        raise ValueError("Different topology in pose and query")
    full_xyz = None
    out = dict.fromkeys(ARMS, 0.)
    other = 0.
    out.update(sites=len(sites), active_sites=0, contacts=0)
    for i, ligand in xyz.items():
        if ligand.kind not in RADII:
            continue
        active = []
        for offset in rec.tree.query_ball_point(ligand.xyz, 3.6):
            j = rec.ids[offset]
            protein = rec.atoms[j]
            lh = owners[i] if protein.kind in ACCEPTORS else ()
            rh = rec.hydrogens[j] if ligand.kind in ACCEPTORS else ()
            if not (lh or rh):
                continue
            a = magnitude(ligand, protein)
            if not a:
                continue
            reverse = max((quality(protein.xyz, rec.atoms[h].xyz, ligand.xyz) for h in rh), default=0.)
            if i in sites and lh:
                active.append((protein.xyz, a, reverse))
            else:
                forward = max((quality(ligand.xyz, xyz[h].xyz, protein.xyz) for h in lh), default=0.)
                other += a*(1-max(forward, reverse))
        if not active:
            continue
        out["active_sites"] += 1
        out["contacts"] += len(active)
        site = sites[i]
        circle, around, axis = ring(site, xyz, n, phase)
        oh = xyz[site.hydrogen].xyz - ligand.xyz
        flipped_h = ligand.xyz + 2*np.dot(oh,axis)*axis - oh
        points = np.vstack((xyz[site.hydrogen].xyz, flipped_h, circle))
        p = sum(a*(1-np.maximum(reverse,quality(ligand.xyz,points,acceptor)))
                for acceptor,a,reverse in active)
        out["fixed"] += p[0]
        out["turned"] += p[1]
        penalty = p[2:]
        out["minimum"] += float(penalty.min())
        out["uniform"] += marginal(penalty, np.zeros(n))
        if any(x not in serial for x,_ in site.terms) and full_xyz is None:
            full_xyz = hydrogenated_positions(mol, serial, xyz)
        prior = torsional_prior(site, around, axis, xyz, serial, full_xyz)
        out["bonded"] += marginal(penalty, prior)
    for arm in ("fixed","turned","uniform","bonded"):
        out["full_"+arm] = other + out[arm]
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receptor",type=Path)
    parser.add_argument("ligand",type=Path)
    parser.add_argument("pool",type=Path)
    parser.add_argument("--scale",required=True,type=float)
    args = parser.parse_args()
    if not 0 < args.scale <= 1:
        parser.error("Native pose normalization must be in (0,1]")
    rec = Receptor(args.receptor.read_text())
    query = args.ligand.read_text()
    writer = None
    for rank,block in enumerate((b for b in args.pool.read_text().split("ENDMDL") if b.strip()),1):
        native = float(next(l for l in block.splitlines() if l.startswith("REMARK VINA RESULT:")).split()[3])
        terms = score(rec, query, block)
        row = dict(rank=rank,native=native,**{arm:native+args.scale*terms[arm] for arm in ARMS},
                   **{k:v for k,v in terms.items() if k not in ARMS})
        if writer is None:
            writer = csv.DictWriter(sys.stdout,fieldnames=row,lineterminator="\n")
            writer.writeheader()
        writer.writerow(row)
    if writer is None:
        parser.error("Empty pose pool")
