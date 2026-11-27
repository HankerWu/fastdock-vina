"""A rigid ligand wider than the map diagonal cannot yield a legal released pose."""
import argparse
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("binary",type=Path)
    binary=p.parse_args().binary.resolve()
    with TemporaryDirectory() as tmp:
        root=Path(tmp)
        prefix=root/"receptor"
        hdr="GRID_PARAMETER_FILE x\nGRID_DATA_FILE x\nMACROMOLECULE x\nSPACING 1\nNELEMENTS 2 2 2\nCENTER 0 0 0\n"
        for t in ("C","e","d","C_H"):
            (root/f"receptor.{t}.map").write_text(hdr+"0\n"*27)
        atom="ATOM      1   C1 UNL     1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
        other="ATOM      2   C2 UNL     1       6.000   0.000   0.000  1.00  0.00     0.000 C\n"
        for sf in ("ad4","vina","vinardo"):
            for feasible,lines in ((False,atom+other),(True,atom)):
                ligand=root/"ligand.pdbqt";out=root/"output.pdbqt"
                ligand.write_text("ROOT\n"+lines+"ENDROOT\nTORSDOF 0\n")
                out.write_text("stale result")
                result=subprocess.run([binary,"--scoring",sf,"--maps",prefix,"--ligand",ligand,
                    "--out",out,"--cpu","1","--exhaustiveness","2","--max_evals","500","--seed","331"],
                    check=True,capture_output=True,text=True)
                if feasible:
                    assert "MODEL 1" in out.read_text(), (sf,result.stderr)
                else:
                    assert out.read_bytes()==b"", (sf,out.read_text())
                    assert "No saved poses lie inside the affinity map domain" in result.stderr
    print("All three map-only force fields release feasible poses and replace impossible results with zero models.")

if __name__=="__main__": main()
