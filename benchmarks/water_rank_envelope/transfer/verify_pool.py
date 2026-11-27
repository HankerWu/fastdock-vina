"""Independent protein fit, chemical symmetry, and water-state dominance checks."""
import argparse
import itertools
import json
import lzma
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.spatial.transform import Rotation

HERE=Path(__file__).resolve().parent
BENCH=HERE.parents[1]
sys.path.insert(0,str(BENCH/"cubic_grid/human_panel"))
from verify import reference


def independent_geometry(table, raw):
    qc, qc0, truth = reference(raw,"4EY5","HUP")
    rc, rc0, _ = reference(raw,"4EY6","GNT")
    shared=sorted(key for key in qc.keys() & rc if
                  np.linalg.norm(qc[key]-qc0)<14 and np.linalg.norm(rc[key]-rc0)<14)
    a,b=(np.array([points[k] for k in shared]) for points in (qc,rc))
    R,_=Rotation.align_vectors(b-b.mean(axis=0),a-a.mean(axis=0))
    ref=R.apply(truth.GetConformer().GetPositions()-a.mean(axis=0))+b.mean(axis=0)
    distances=[]
    for block in table.read_text().split("ENDMDL")[:-1]:
        lines=block.splitlines()
        smiles=next(l.split()[2] for l in lines if l.startswith("REMARK SMILES "))
        mol=Chem.MolFromSmiles(smiles)
        indices=[int(i) for l in lines if l.startswith("REMARK SMILES IDX ") for i in l.split()[3:]]
        pairs=sorted(zip(indices[::2],indices[1::2],strict=True))
        serials={int(l[6:11]):np.array([float(l[k:k+8]) for k in (30,38,46)]) for l in lines if l.startswith("ATOM")}
        xyz=np.array([serials[j] for _,j in pairs])
        conf=Chem.Conformer(len(xyz))
        for i,x in enumerate(xyz):
            conf.SetAtomPosition(i,x)
        mol.AddConformer(conf)
        initial=Chem.MolToSmiles(mol)
        Chem.AssignStereochemistryFrom3D(mol,replaceExistingTags=True)
        if Chem.MolToSmiles(mol)!=initial:
            raise ValueError("Changed stereo in published mode")
        maps=truth.GetSubstructMatches(mol,uniquify=False,useChirality=True)
        distances.append(min(np.linalg.norm(xyz-ref[list(mapping)],"fro")/np.sqrt(len(xyz)) for mapping in maps))
    return np.array(distances)


def verify(directory):
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        raw=root/"raw"
        raw.mkdir()
        (raw/"4EY5.pdb").write_bytes(lzma.decompress((HERE/"4EY5.pdb.xz").read_bytes()))
        (raw/"HUP.cif").write_bytes((HERE/"HUP.cif").read_bytes())
        with tarfile.open(BENCH/"receptor_amide_states/transfer/transfer_evidence.tar.xz") as ar:
            for name in ("4EY6.pdb","GNT.cif"):
                (raw/name).write_bytes(ar.extractfile("novel/raw/"+name).read())
        distances=independent_geometry(directory/"dry.pdbqt",root)
    mode=pd.read_csv(directory/"modes.csv",float_precision="round_trip")
    factors=pd.read_csv(directory/"factors.csv",float_precision="round_trip")
    results=[]
    RT=0.00198720425864083*298.15
    for order,subset in factors.groupby("order"):
        data=mode[mode.order==order]
        error=max(abs(distances-data.rmsd.to_numpy()))
        if error>1e-8:
            raise ValueError(f"RMSD differs by {error}")
        avail=subset.pivot(index="mode",columns="site",values="available").to_numpy()
        fixed=data.dry.to_numpy()-RT*np.log(.25+.75*avail).sum(axis=1)
        assert np.allclose(data.hydrated,fixed,rtol=0,atol=1e-12)
        variable=avail[:,np.ptp(avail,axis=0)>0]
        n=variable.shape[1]
        if n>24:
            raise ValueError("Exact independent check too large")
        logv=np.log(np.where(variable>0,variable,1.))
        witness=[(i,int(row.rival)-1) for i,row in enumerate(data.itertuples(index=False))
                 if row.status=="single_rival"]
        minimum=np.full(len(witness),np.inf)
        for start in range(0,2**n,2048):
            ints=np.arange(start,min(start+2048,2**n),dtype=np.uint32)
            state=((ints[:,None]>>np.arange(n))&1).astype(float)
            e=data.dry.to_numpy()[None,:]-RT*(state @ logv.T)
            e[(state @ (variable==0).T)>0]=np.inf
            for k,(i,j) in enumerate(witness):
                ok=np.isfinite(e[:,i])
                if ok.any():
                    minimum[k]=min(minimum[k],(e[ok,i]-e[ok,j]).min())
        reported=data[data.status=="single_rival"].margin.to_numpy()
        assert np.allclose(minimum,reported,rtol=0,atol=1e-10) and np.all(minimum>.001)
        results.append(dict(order=int(order),distinct_binary_dimensions=n,states=2**n,
                            verified_rivals=len(witness),rmsd_max_error_A=error,
                            minimum_gap_error=float(max(abs(minimum-reported),default=0))))
    return results


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory",type=Path)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    a.output.write_text(json.dumps(verify(a.directory),indent=2)+"\n")
