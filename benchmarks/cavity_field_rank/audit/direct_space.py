"""Post hoc only. Audit geometry, and remove the spatial lattice from the fixed-pose receptor sum."""
import argparse
import csv
import json
import subprocess
import sys
import tarfile
import time
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from scipy.stats import spearmanr

BENCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(BENCH / "cubic_grid/human_panel"))
from verify import reference

CODES = dict(zip(("6YQW","6Y7L","3FMZ","1RBP","4GHI","5TBM"), ("82I","QMG","2T1","RTL","0X3","79A"), strict=True))
PDB = {"82I":("6YQW","6Y7L"), "3FMZ-A":("3FMZ","1RBP"), "3FMZ-HA":("3FMZ","1RBP"),
       "4GHI":("4GHI","5TBM"), "5TBM":("5TBM","4GHI")}
RADII = np.array([1.9,1.9,1.8,1.8,1.8,1.8,1.7,1.7,1.7,1.7,2,2.1,1.5,1.8,2,2.2,2.2,2.3,1.2])

def table(path, data):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)

def deposited_frames(work):
    for archive, group in [
        ("torsion_proposals/fixtures/inputs.tar.xz", ("6YQW","6Y7L")),
        ("cubic_grid/human_panel/depositions.tar.xz", ("1RBP","3FMZ","4GHI","5TBM")),
    ]:
        with tarfile.open(BENCH / archive) as tar:
            for pdb in group:
                for filename in (pdb+".pdb", CODES[pdb]+".cif"):
                    tar.extract("raw/"+filename, work, filter="data")
    return {p:reference(work, p, c) for p,c in CODES.items()}

def equivalent_targets(name, ligfile, frames):
    q,r = PDB[name]
    qa,qc,crystal = frames[q]
    ra,rc,_ = frames[r]
    if r == "6Y7L":
        ra = {(n-116,code,aa):v for (n,code,aa),v in ra.items()}
    keys = sorted(k for k in qa.keys() & ra.keys() if np.linalg.norm(qa[k]-qc)<14 and np.linalg.norm(ra[k]-rc)<14)
    assert len(keys) >= 12
    x,y = [np.array([ca[k] for k in keys]) for ca in (qa,ra)]
    rot,_ = Rotation.align_vectors(y-y.mean(0), x-x.mean(0))
    geom = rot.apply(crystal.GetConformer().GetPositions()-x.mean(0))+y.mean(0)
    idx = {}
    for line in ligfile.read_text().splitlines():
        if line.startswith("REMARK SMILES IDX "):
            v = list(map(int, line.split()[3:]))
            idx.update(zip(v[::2],v[1::2],strict=True))
        elif line.startswith("REMARK SMILES "):
            mol = Chem.MolFromSmiles(line.split()[2])
    serial = [idx[i+1] for i in range(mol.GetNumAtoms())]
    oxygen = None
    if name=="3FMZ-A":
        ((_,a,b),) = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])[O-]"))
        oxygen = (a,b)
        mol.GetAtomWithIdx(b).SetFormalCharge(0)
        mol.GetAtomWithIdx(b).SetNoImplicit(False)
        Chem.SanitizeMol(mol)
    mappings = crystal.GetSubstructMatches(mol, useChirality=True, uniquify=False, maxMatches=10001)
    assert 0 < len(mappings) < 10001
    refs = geom[np.array(mappings)]
    if oxygen:
        a,b = oxygen
        swapped = refs.copy()
        swapped[:,[a,b]] = refs[:,[b,a]]
        refs = np.concatenate([refs,swapped])
    return serial, refs

def distances(x, ref):
    return np.sqrt(np.min(np.mean(np.sum((x[:,None]-ref[None])**2, axis=-1),axis=-1),axis=-1))

def direct_scores(xyz, rec, rtypes, ltypes, center, size):
    # Same 8 A cutoff, midpoint sampling of squared distance, per-atom curl and native wall.
    # Here each pair uses its Cartesian separation. No map nodes or trilinear interpolation.
    n, heavy = xyz.shape[:2]
    half = .375 * np.ceil(size/.375) / 2
    bounded = np.clip(xyz,center-half,center+half)
    wall = 1e6*np.abs(xyz-bounded).sum(axis=(1,2))
    flat = bounded.reshape(-1,3)
    tree = cKDTree(rec)
    neighbours = cKDTree(flat).sparse_distance_matrix(tree, 8., output_type="coo_matrix")
    i, j, r = neighbours.row, neighbours.col, neighbours.data
    ligand = ltypes[i % heavy]
    receptor = rtypes[j]
    radius = RADII[ligand]+RADII[receptor]
    hydrophobic = [0,12,13,14,15]
    hp = np.isin(ligand,hydrophobic)&np.isin(receptor,hydrophobic)
    donor, acceptor = [3,5,7,9,18],[4,5,8,9]
    hb = ((np.isin(ligand,donor)&np.isin(receptor,acceptor))|
          (np.isin(ligand,acceptor)&np.isin(receptor,donor)))
    binned = np.floor(r*r*32)
    energy = np.zeros(len(r))
    for offset in (0,1):
        gap = np.sqrt((binned+offset)/32)-radius
        pair = (-.035579*np.exp(-(gap/.5)**2) -.005156*np.exp(-((gap-3)/2)**2)
                +.840245*np.minimum(gap,0)**2 -.035069*hp*np.clip(1.5-gap,0,1)
                -.587439*hb*np.clip(-gap/.7,0,1))
        energy += .5*pair*(binned+offset < 8**2*32)
    atom = np.bincount(i, weights=energy, minlength=len(flat))
    atom /= 1+np.maximum(atom,0)/10
    closest = tree.query(xyz)[0]
    return atom.reshape(n,heavy).sum(1)+wall, closest.min(1), (closest<2.).sum(1)

def run(pose, auditwork, binary, out):
    start, cpu = time.perf_counter(),time.process_time()
    auditwork.mkdir(parents=True,exist_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    systems = json.loads((pose/"inputs.json").read_text())["cells"]
    frames = deposited_frames(auditwork)
    rows, obs, details, bounds = [],[],{},[]
    checked, largest_pair_error = 0,0.
    for f in sorted((pose/"banks").glob("*.source.npz")):
        case = f.name.split(".")[0]
        drug,seed = case.split("_")
        c = systems[drug]
        inputs = [pose/c[k] for k in ("receptor","ligand","bonds")]
        record = subprocess.check_output([binary,*inputs],text=True)
        receptor = np.array([r.split()[1:] for r in record.splitlines() if r.startswith("R ")], float)
        ligand = {int(r.split()[1]):int(r.split()[2]) for r in record.splitlines() if r.startswith("L ")}
        serial, refs = equivalent_targets(drug,inputs[1],frames)
        with np.load(f) as src, np.load(pose/"banks"/(case+".npz")) as lib:
            order = src["serials"].tolist()
            assert len(ligand)==len(order) and sorted(order)==sorted(ligand)
            types = np.array([ligand[s] for s in order])
            reverse_map = {s:ligand[t] for s,t in zip(sorted(order),np.roll(sorted(order),-(len(order)//2)),strict=True)}
            shuffled = np.array([reverse_map[s] for s in order])
            coords = (src["xyz"]-src["xyz"].mean(0)) @ lib["matched_rotation"].transpose(0,2,1)+lib["center"][:,None]
            src_xyz = src["xyz"][[order.index(s) for s in serial]]
            err = distances(coords[:,[order.index(s) for s in serial]],refs)
            shape_winner = int(lib["matched_best"])
        # Label-assisted lower bound for exactly this internal shape. No ligand
        # superposition enters the reported prediction RMSDs.
        local = src_xyz - src_xyz.mean(0)
        shape_limits = []
        for target in refs:
            centered = target - target.mean(0)
            rotation, _ = Rotation.align_vectors(centered, local)
            shape_limits.append(np.linalg.norm(centered-rotation.apply(local))/np.sqrt(len(local)))
        centers = coords.mean(1)
        centroid_error = np.linalg.norm(centers-refs[0].mean(0), axis=1)
        # All chemistry-preserving mappings share a centroid. Recombining any
        # retained orientation with any retained centre gives the following exact minimum.
        centered_loss = np.maximum(0,err**2-centroid_error**2)
        crossed = np.sqrt(centroid_error.min()**2+centered_loss.min())
        # The reference enters only here. For a random bijection each ordered
        # pair occurs with probability 1/N, regardless of dependence across rows.
        n_pairs = int(np.searchsorted(np.sort(centered_loss),4-centroid_error**2,side="right").sum())
        bounds.append(dict(case=case,ideal_rigid_A=min(shape_limits),center_min_A=centroid_error.min(),
                           orientation_min_A=np.sqrt(centered_loss.min()),any_center_orientation_A=crossed,
                           matched_min_A=err.min(),present_le2=int(np.count_nonzero(err<=2)),
                           le2_pairs_in_full_cross=n_pairs,uniform_reassignment_expected_hits=n_pairs/len(err)))
        near = err<=2
        for arm,atom_types in (("true",types),("perm",shuffled)):
            exact, closest, clashes = direct_scores(coords,receptor[:,1:],receptor[:,0].astype(int),atom_types,
                                                    np.array(c["center"]),c["size"])
            grid = np.loadtxt(pose/"raw"/f"{case}_{arm}.csv",delimiter=",",skiprows=1)[:,1]
            a,b = int(grid.argmin()),int(exact.argmin())
            subset = sorted({0, a, b, shape_winner, int(err.argmin()), 1023})
            half = .375*np.ceil(c["size"]/.375)/2
            text = []
            for k in subset:
                bound = np.clip(coords[k],np.array(c["center"])-half,np.array(c["center"])+half)
                text.append(f"{len(bound)} {1e6*np.abs(coords[k]-bound).sum():.17g}")
                for t, point in zip(atom_types,bound,strict=True):
                    text.append(f"{t} "+" ".join(f"{v:.17g}" for v in point))
            native = subprocess.check_output([binary,*inputs,"--pairs"],input="\n".join(text)+"\n",text=True)
            delta = np.max(np.abs(np.fromstring(native,sep=" ")-exact[subset]))
            if delta > 2e-7:
                raise ValueError(f"Native Cartesian table mismatch {case} {arm}: {delta}")
            checked += len(subset)
            largest_pair_error = max(delta,largest_pair_error)
            ranked = np.argsort(exact,kind="stable")
            ranks = np.argsort(grid,kind="stable")
            rows.append(dict(case=case,arm=arm,lib_min_A=err.min(),n_near=near.sum(),shape_winner=shape_winner,
                grid_winner=a,direct_winner=b,grid_winner_A=err[a],direct_winner_A=err[b],
                grid_first_near=next((i+1 for i,k in enumerate(ranks) if near[k]),""),
                direct_first_near=next((i+1 for i,k in enumerate(ranked) if near[k]),""),
                median_abs_grid_minus_direct_kcal_mol=np.median(abs(grid-exact)),
                spearman_grid_direct=spearmanr(grid,exact).statistic,grid_winner_min_contact_A=closest[a],
                median_near_min_contact_A=np.median(closest[near]) if near.any() else "",
                grid_winner_atoms_under_2A=clashes[a],near_median_atoms_under_2A=np.median(clashes[near]) if near.any() else ""))
            details[case+"_"+arm] = np.c_[err,grid,exact,closest,clashes]
        for arm in ("official","official-high","pro","pro-high","shape","true","perm"):
            blocks=(pose/"raw"/f"{case}_{arm}.pdbqt").read_text().split("ENDMDL")
            for k,block in enumerate(blocks[:-1],1):
                atoms={int(l[6:11]):[float(l[j:j+8]) for j in (30,38,46)] for l in block.splitlines() if l.startswith(("ATOM  ","HETATM"))}
                rms=distances(np.array([[atoms[s] for s in serial]]),refs)[0]
                obs.append(dict(case=case,arm=arm,rank=k,rmsd_A=rms))
    table(out/"direct.csv",rows)
    table(out/"representation.csv",bounds)
    table(out/"independent_modes.csv", obs)
    np.savez_compressed(out/"direct.npz",**details)
    (out/"time.json").write_text(json.dumps(dict(wall_s=time.perf_counter()-start,cpu_self_s=time.process_time()-cpu,
        native_direct_spotchecks=checked,max_native_direct_error_kcal_mol=largest_pair_error,
        columns=["RMSD_A","grid_kcal_per_mol","direct_kcal_per_mol","min_contact_A","heavy_atoms_under_2A"]),indent=2)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("pose","work","types_binary","out"):
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    run(a.pose,a.work,a.types_binary,a.out)
