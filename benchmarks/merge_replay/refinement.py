"""Run or verify the single preregistered slot-substitution in the normal CLI."""
import argparse
import csv
import json
import lzma
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stage_local_hunt"))
from coverage import specs, run_panel

KEY = "reciprocal-1RBP-3FMZ-85237-candidate"
ORDER = ("control", "substitution", "substitution-repeat", "control-repeat")


def parse(d, job):
    _, _, q, _, ref, _, _ = job
    base = (HERE / "raw" / (KEY+".tsv.xz")).read_bytes()
    original = lzma.decompress(base).decode().splitlines()
    merge = (d/"merge.tsv").read_text().splitlines()
    assert merge[:-1] == original and merge[-1].startswith("U\t")
    selection = list(map(int,merge[-1].split("\t")[1:]))
    states = [r.split("\t")[4] for r in original if r.startswith("C\t")]
    native = list(map(int,next(r for r in original if r.startswith("I\t")).split("\t")[1:]))
    expected = native.copy()
    if d.name.startswith("substitution"):
        expected[0] = 20
    assert native[0] == 18 and expected == selection
    records = [r.split("\t") for r in (d/"refine.tsv").read_text().splitlines()]
    inputs = {int(r[1]):r[3] for r in records if r[0]=="I"}
    assert [inputs[i] for i in range(len(inputs))] == [states[j] for j in expected]
    refine = {r[2]:(expected[int(r[1])],r) for r in records if r[0]=="O"}
    final = {r[3]:int(r[1])+1 for r in records if r[0]=="F"}
    out = d/(KEY+".pdbqt")
    observed = run_panel.replay.reciprocal.measure(out,q["serial"],ref)
    rounded = [float(l.split()[3]) for l in out.read_text().splitlines() if l.startswith("REMARK VINA RESULT:")]
    rows = []
    for r in records:
        if r[0]!="R":
            continue
        ident, post = refine[r[2]]
        rank = final.get(r[2])
        published_rank = rank if rank and rank <= len(observed) else None
        row = dict(run=d.name, input_id=ident, refined_sort_rank=int(r[1])+1, public_rank=published_rank,
                   public_crystal_A=observed[rank-1] if published_rank else None,
                   public_score=rounded[rank-1] if published_rank else None, evals=int(post[3]),
                   refined_grid=float(post[4]), refined_inter_pairs=float(post[5]),
                   refined_intra_pairs=float(post[6]), refined_lig_intra=float(post[7]),
                   native_intra_reference=float(r[3]), control_intra_reference=float(r[4]),
                   control_reference_score=float(r[5]))
        names = ("native_score","lig_grid","inter_pairs","flex_grid","intra_pairs","lig_intra","conf_independent","unbound")
        assert len(r[6:]) == len(names)
        row.update({k:float(v) for k,v in zip(names,r[6:])})
        if published_rank:
            assert abs(row["native_score"] - row["public_score"]) <= .000501
        rows.append(row)
    assert len(rows)==9 and {r["input_id"] for r in rows}==set(expected)
    return rows


def drive(tmp, out, binary=None):
    job = next(job for _,job in specs(tmp) if (job[2]["pdb"], job[3]["pdb"], job[5])==("1RBP","3FMZ",85237))
    _, inputs, q, r, _, _, _ = job
    anchor = None
    receipts = []
    result = []
    src = json.loads((HERE / "raw/costs.json").read_text())
    previous = next(i for i in src if i["key"]==KEY)
    if binary:
        out.mkdir(parents=True,exist_ok=False)
    for tag in ORDER:
        d = out/tag
        if binary:
            d.mkdir()
            os.environ["FD_MERGE_TRACE"]=str(d/"merge.tsv")
            os.environ["FD_BFGS_TRACE"]=str(d/"refine.tsv")
            if anchor is not None:
                os.environ["FD_COMMON_REFERENCE"]=format(anchor,".17g")
            else:
                os.environ.pop("FD_COMMON_REFERENCE",None)
            if tag.startswith("substitution"):
                os.environ["FD_SUBSTITUTE_20"]="yes"
            else:
                os.environ.pop("FD_SUBSTITUTE_20",None)
            commands = []
            execute = run_panel.execute
            def timed(command,*args):
                commands.append(command)
                return execute(command,*args)
            run_panel.execute = timed
            try:
                row,_ = run_panel.dock(("candidate",binary,d,job))
            finally:
                run_panel.execute = execute
            receipts.append(dict(binary_sha256=run_panel.digest(binary),tag=tag,
                                 ligand_sha256=run_panel.digest(inputs/q["pdb"]/"ligand.pdbqt"),
                                 receptor_sha256=run_panel.digest(inputs/r["pdb"]/"receptor.pdbqt"),
                                 args=commands[0],**row))
            (out/"receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")
        else:
            row = next(r for r in json.loads((out/"receipts.json").read_text()) if r["tag"]==tag)
        assert row["status"]=="ok"
        assert run_panel.digest(d/(KEY+".pdbqt"))==row["pool_sha256"]
        if tag.startswith("control"):
            assert row["pool_sha256"] == previous["pool_sha256"]
        assert all(row[k]==previous[k] for k in ("hunt_evals","second_evals","explicit_calls"))
        parsed = parse(d,job)
        assert sum(r["evals"] for r in parsed)==row["explicit_evals"]
        if anchor is None:
            anchor = parsed[0]["native_intra_reference"]
        assert all(anchor == r["control_intra_reference"] for r in parsed)
        result.extend(parsed)
    def stripped(rows):
        return [{k:v for k,v in r.items() if k!="run"} for r in rows]
    for mode in ("control","substitution"):
        orig = [r for r in result if r["run"]==mode]
        rep = [r for r in result if r["run"]==mode+"-repeat"]
        assert stripped(orig)==stripped(rep)
        for f in ("merge.tsv","refine.tsv",KEY+".pdbqt"):
            assert (out/mode/f).read_bytes()==(out/(mode+"-repeat")/f).read_bytes()
    control = {r["input_id"]:r for r in result if r["run"]=="control"}
    subs = {r["input_id"]:r for r in result if r["run"]=="substitution"}
    for i in control.keys() & subs.keys():
        assert all(control[i][k]==subs[i][k] for k in ("control_reference_score","evals","refined_grid","refined_lig_intra"))
    path = out/"scores.csv"
    if binary:
        run_panel.csvfile(path,result)
    else:
        assert list(csv.DictReader(path.open()))==[{k:("" if v is None else str(v)) for k,v in r.items()} for r in result]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    parser.add_argument("--binary",type=Path)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        drive(Path(tmp),args.out.resolve(),args.binary.resolve() if args.binary else None)
