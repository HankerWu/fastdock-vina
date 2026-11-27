"""Repeat 1FKG compound 8 from separately prepared inputs for each online arm."""
import argparse
import json
from pathlib import Path
import sys

from trial import HERE, LIVE, native


def run(work, fresh, engine, output):
    fresh.mkdir(parents=True,exist_ok=True)
    script=HERE/"trial.py"
    wall,cpu=native.time_one([sys.executable,str(script),"prepare",str(fresh)],fresh,fresh/"prepare.log")
    assert native.sha(fresh/"inputs.json")==native.sha(work/"inputs.json")
    rows=[dict(action="input_preparation",same_files=1,pool_sha256="",wall_s=wall,cpu_s=cpu)]
    for arm in LIVE:
        src=work/"pools"/arm/"1FKG_8"
        target=fresh/"pools"/arm/"1FKG_8"
        target.parent.mkdir(parents=True,exist_ok=True)
        target.with_suffix(".conf").write_bytes(src.with_suffix(".conf").read_bytes())
        wall,cpu=native.time_one([sys.executable,str(script),"live",str(fresh),"1FKG","8",arm,str(engine)],
                                 fresh,target.with_suffix(".log"))
        for suffix in (".source",".q",".npz",".transform",".stages",".history",".xyz",".pdbqt"):
            assert native.sha(src.with_suffix(suffix))==native.sha(target.with_suffix(suffix))
        def without_clocks(p):
            data=json.loads(p.with_suffix(".rank").read_text())
            return {k:v for k,v in data.items() if not k.endswith(("_cpu_s","_wall_s")) and k!="radius_A2"}
        assert without_clocks(src)==without_clocks(target)
        rows.append(dict(action=arm,same_files=8,pool_sha256=native.sha(target.with_suffix(".pdbqt")),wall_s=wall,cpu_s=cpu))
    native.table(output,rows)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("work","fresh","engine","output"): p.add_argument(key,type=Path)
    a=p.parse_args()
    run(a.work.resolve(),a.fresh.resolve(),a.engine.resolve(),a.output.resolve())
