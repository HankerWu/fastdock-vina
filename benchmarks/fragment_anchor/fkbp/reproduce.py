"""Clean preparation and one-pass replay of each candidate arm, keeping production and audit clocks separate."""
import argparse
from pathlib import Path
import sys

from assay_runner import LIVE, sha, table, time_one


def replay(source,fresh,engine,out):
    dest=Path(__file__).resolve().parent
    fresh.mkdir(parents=True,exist_ok=True)
    prep_wall,prep_cpu=time_one([sys.executable,str(dest/"assay_runner.py"),"prepare",str(fresh)],fresh,fresh/"prepare.log")
    rows=[dict(action="input_archive_to_prepared",files=0,wall_s=prep_wall,cpu_s=prep_cpu,pool_sha256="")]
    for arm in LIVE:
        original=source/"pools"/arm/"1FKG_8"
        stem=fresh/"pools"/arm/original.name
        stem.parent.mkdir(parents=True,exist_ok=True)
        stem.with_suffix(".conf").write_bytes(original.with_suffix(".conf").read_bytes())
        cmd=[sys.executable,str(dest/"assay_runner.py"),"live",str(fresh),"1FKG","8",arm,str(engine)]
        wall,cpu=time_one(cmd,fresh,stem.with_suffix(".log"))
        suffixes=(".q",".source",".npz",".rank",".transform",".stages",".history",".xyz",".pdbqt")
        for suffix in suffixes: assert stem.with_suffix(suffix).read_bytes()==original.with_suffix(suffix).read_bytes()
        rows.append(dict(action=arm,files=len(suffixes),wall_s=wall,cpu_s=cpu,pool_sha256=sha(stem.with_suffix(".pdbqt"))))
    table(out,rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for s in ("source","fresh","engine","out"):p.add_argument(s,type=Path)
    a=p.parse_args()
    replay(a.source.resolve(),a.fresh.resolve(),a.engine.resolve(),a.out.resolve())
