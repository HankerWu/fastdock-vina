"""Archive raw predictions, timing logs, graphs, and one failed cost bracket with exact-byte deduplication."""
import argparse
import hashlib
from pathlib import Path
import tarfile
import run_panel as panel


def pack(p,dest):
    entries=[]
    for name in ("dump","raw","timed","inputs","official_graph","failure"):
        entries.extend(path for path in (p/name).rglob("*") if path.is_file() or path.is_symlink())
    for pattern in ("*.csv","*.json","*.sha256","evaluation_time.txt"):
        entries.extend(p.glob(pattern))
    archive=dest/"raw.tar.xz"
    hashes=[]
    content={}
    with tarfile.open(archive,"w:xz",preset=6) as tar:
        for f in sorted(set(entries)):
            name=f.relative_to(p).as_posix()
            sha=panel.digest(f)
            hashes.append(f"{sha}  {name}\n")
            info=tar.gettarinfo(f,arcname=name)
            if not f.is_symlink() and sha in content:
                info.type=tarfile.LNKTYPE
                info.linkname=content[sha]
                info.size=0
                tar.addfile(info)
            elif f.is_symlink():
                tar.addfile(info)
            else:
                content[sha]=name
                with f.open("rb") as stream: tar.addfile(info,stream)
    (dest/"results/raw_files.sha256").write_text("".join(hashes))
    (dest/"ARCHIVE.sha256").write_text(f"{panel.digest(archive)}  raw.tar.xz\n")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("work",type=lambda s:Path(s).resolve())
    a=cli.parse_args()
    pack(a.work,panel.HERE)
