"""Create deterministic per-arm raw archives or verify and restore into a clean prepared directory."""
import argparse
import json
from pathlib import Path
import tarfile
from assay_runner import ARMS, sha


def create(work, dest, *, arms=ARMS, replay_files=()):
    dest.mkdir(parents=True,exist_ok=True)
    parts={arm:sorted((work/"pools"/arm).glob("*")) for arm in arms}
    # Production replay PDBQTs are byte-identical to the corresponding pools. Keep logs, hashes and clocks.
    parts["replay"]=[*sorted((work/"ordinary_pro").glob("*.log")),work/"unmodified_pro.csv",
                     *sorted((work/"failures").glob("*")), *replay_files]
    manifest={}
    for arm,files in parts.items():
        fname=dest/(arm+".tar.xz")
        entries={}
        with tarfile.open(fname,"w:xz",preset=4) as tar:
            for f in files:
                name=f.relative_to(work).as_posix()
                entries[name]=sha(f)
                info=tar.gettarinfo(str(f),arcname=name)
                info.uid=info.gid=info.mtime=0
                info.uname=info.gname=""
                with f.open("rb") as src: tar.addfile(info,src)
        manifest[fname.name]=dict(sha256=sha(fname),files=entries)
    (dest/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")


def unpack(work,src):
    for name,row in json.loads((src/"manifest.json").read_text()).items():
        f=src/name
        assert sha(f)==row["sha256"]
        with tarfile.open(f) as tar: tar.extractall(work,filter="data")
        for p,h in row["files"].items(): assert sha(work/p)==h


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("pack","unpack"))
    p.add_argument("work",type=Path)
    p.add_argument("raw",type=Path)
    a=p.parse_args()
    (create if a.action=="pack" else unpack)(a.work.resolve(),a.raw.resolve())
