"""Recheck the field increment on all supplied ABL tutorial poses and maps."""
import argparse
from pathlib import Path
import tempfile

from fixed import atoms, correction, score, write_csv


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ("solution","old","new","output"):
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    rigid=a.solution/"1fpu_receptor_rigid.pdbqt"
    flex=a.solution/"1fpu_receptor_flex.pdbqt"
    frozen=a.solution/"1fpu_ligand_flex_ad4_out.pdbqt"
    rigid_ids={line[12:26] for line in atoms(rigid.read_text())}
    if rigid_ids & {line[12:26] for line in atoms(flex.read_text())}:
        raise ValueError("Tutorial split duplicates an atom")
    rows=[]
    with tempfile.TemporaryDirectory() as tmp:
        lig,side=[Path(tmp)/name for name in ("ligand.pdbqt","flex.pdbqt")]
        for rank,block in enumerate(frozen.read_text().split("ENDMDL")[:-1],1):
            first,rest=block.split("\nROOT\n",1)[1].split("BEGIN_RES",1)
            ligand,sidechain="ROOT\n"+first,"BEGIN_RES"+rest
            lig.write_text(ligand)
            side.write_text(sidechain)
            pred=correction(ligand,sidechain)
            old,new=[score(binary.resolve(),a.solution.resolve()/"1fpu_receptor_rigid",lig,side)
                     for binary in (a.old,a.new)]
            err=max(abs(new[k]-old[k]-pred) for k in ("side","score"))
            if err>.00101 or any(new[k]!=old[k] for k in ("grid","internal","torsion")):
                raise ValueError(f"Tutorial model {rank} mismatch {err}")
            rows.append(dict(rank=rank,old_score=old["score"],new_score=new["score"],old_side=old["side"],
                             new_side=new["side"],analytic_delta=pred,rounding_error=err))
    write_csv(a.output,rows)


if __name__=="__main__":
    main()
