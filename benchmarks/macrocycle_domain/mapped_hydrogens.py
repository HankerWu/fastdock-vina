"""Use a rigid copy of the official human BACE1 ligand with only a mapped HD outside."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull

from boundary import call


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("solution", "output", "before", "after"):
        p.add_argument(key, type=Path)
    a = p.parse_args()
    a.output.mkdir(exist_ok=True, parents=True)
    prefix = a.solution / "BACE_1_receptor"
    lines = (a.solution/"BACE_1_ligand.pdbqt").read_text().splitlines()
    atom_lines = [i for i,l in enumerate(lines) if l.startswith(("ATOM  ","HETATM"))]
    types = np.array([lines[i].split()[-1] for i in atom_lines])
    charge = np.array([float(lines[i][70:76]) for i in atom_lines])
    xyz = np.array([[float(lines[i][j:j+8]) for j in (30,38,46)] for i in atom_lines])
    with prefix.with_suffix(".HD.map").open() as f:
        head = [next(f).split() for _ in range(6)]
    size = float(head[3][1]) * np.array(head[4][1:],dtype=int)
    center = np.array(head[5][1:],dtype=float)
    top = center + size/2
    bottom = center - size/2
    hd = types == "HD"
    # Orient a supporting face of the convex hull of all non-HD sites, including G0.
    # Its normal is chosen only for a positive HD margin, without looking at scores.
    planes = ConvexHull(xyz[~hd]).equations
    normal = planes[(xyz[hd]@planes[:,:3].T+planes[:,3]).max(axis=0).argmax(), :3]
    axis = np.eye(3)[abs(normal).argmin()]
    tangent = np.cross(normal, axis)
    tangent /= np.linalg.norm(tangent)
    rotation = np.stack((normal,tangent,np.cross(normal,tangent)), axis=1)
    assert np.allclose(np.linalg.det(rotation), 1)
    rotated = xyz @ rotation
    gap = rotated[hd,0].max() - rotated[~hd,0].max()
    shifted = rotated + center - (rotated.max(0)+rotated.min(0))/2
    shifted[:,0] += top[0] - (shifted[hd,0].max()+shifted[~hd,0].max())/2
    shifted = np.round(shifted, 3)
    assert np.all((shifted[~hd] >= bottom) & (shifted[~hd] <= top))
    assert np.all(shifted >= bottom) and np.all(shifted[:,1:] <= top[1:])
    outside = np.flatnonzero(shifted[:,0] > top[0])
    assert len(outside) > 0 and np.all(hd[outside])
    # Model is printed to 3 decimals; choose an exact decimal just inside the face.
    clamped = shifted.copy()
    clamped[outside,0] = np.floor(top[0]*1000)/1000
    assert abs(clamped[outside[0],0] - top[0]) < 1e-8, "Need an exactly printable wall for the subtraction"
    wall = float(np.sum((shifted[outside,0] - top[0]) * 1e6 * (1 + charge[outside] + np.abs(charge[outside]))))
    report = []
    for pose_name, coords in (("inside",xyz),("only_HD_outside",shifted),("clamp_HD_to_face",clamped)):
        out = a.output / (pose_name+".pdbqt")
        doc = lines.copy()
        for i,fields in zip(atom_lines,coords):
            doc[i] = doc[i][:30] + ''.join(f'{v:8.3f}' for v in fields) + doc[i][54:]
        out.write_text('\n'.join(doc)+'\n')
        for sf in ("ad4","vina"):
            args = ["--scoring",sf,"--ligand",out,"--score_only"]
            if sf=="ad4":
                args += ["--maps",prefix]
            else:
                args += ["--receptor",prefix.with_suffix(".pdbqt")]
                for d,v,s in zip("xyz",center,size):
                    args += ["--center_"+d,str(v),"--size_"+d,str(s)]
            for arm in ("before","after"):
                rc,score = call(getattr(a,arm),args)
                report.append(dict(pose=pose_name,force_field=sf,operation="score_only",arm=arm,status=rc,score=score))
    values = {(r['pose'],r['force_field'],r['arm']):(r['status'],r['score']) for r in report}
    cl = "clamp_HD_to_face"
    assert values["only_HD_outside","ad4","after"][0] and values["only_HD_outside","ad4","before"][0] == 0
    measured = values["only_HD_outside","ad4","before"][1] - values[cl,"ad4","before"][1]
    assert abs(measured - wall) < 0.003
    for sf in ("ad4","vina"):
        for pose in ("inside",cl):
            assert values[pose,sf,"before"] == values[pose,sf,"after"] and values[pose,sf,"after"][0] == 0
    for arm in ("before","after"):
        assert values["only_HD_outside","vina",arm] == values[cl,"vina",arm]
    for arm in ("before", "after"):
        local_args = ["--scoring", "ad4", "--ligand", a.output/"only_HD_outside.pdbqt",
                      "--maps", prefix, "--local_only", "--out", a.output/(arm+"-local.pdbqt")]
        rc, score = call(getattr(a, arm), local_args)
        assert (rc == 0) == (arm == "before")
        report.append(dict(pose="only_HD_outside", force_field="ad4", operation="local_only",
                           arm=arm, status=rc, score=score))
    with (a.output/"mapped_hydrogens.csv").open("w", newline="") as f:
        w = csv.DictWriter(f,report[0],lineterminator="\n")
        w.writeheader()
        w.writerows(report)
    (a.output/"wall.json").write_text(json.dumps(dict(
        serials=[int(lines[atom_lines[i]][6:11]) for i in outside],
        outward_A=(shifted[outside,0]-top[0]).tolist(), partial_charge=charge[outside].tolist(),
        original_hull_gap_A=gap, expected_wall_kcal_mol=wall, measured_wall_kcal_mol=measured),indent=2)+'\n')


if __name__ == "__main__":
    main()
