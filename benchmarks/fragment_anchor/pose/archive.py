"""Pack small per-input archives and a hash of every replay file from a completed run."""
import argparse
import hashlib
import os
from pathlib import Path
import tarfile

import panel


def pack(work, dest):
    files = {}
    def add(key, paths):
        files[key] = sorted({p for p in paths if p.is_file() or p.is_symlink()})
    add("setup", [work/"inputs.json", *sorted((work/"inputs").glob("*")), *work.glob("*.csv"),
                  *work.glob("*.sha256"), work/"control_binaries.json", work/"observer.sha256",
                  work/"binary_context.json",work/"evaluation_time.txt", *sorted((work/"official_graph").glob("*.csv"))])
    for name in panel.CASES:
        raw = (work/"raw").glob(name+"_*")
        logs = (work/"timed").glob(name+"_*")
        dumps = (work/"dump").glob(name+"_*/*")
        add(name, [*raw,*logs,*dumps])
    allhashes = []
    archives = []
    output = dest/"archives"
    output.mkdir(exist_ok=True)
    for key, group in files.items():
        path = output/(key+".tar.xz")
        with tarfile.open(path,"w:xz",preset=3) as tar:
            for p in group:
                rel = p.relative_to(work).as_posix()
                tar.add(p,arcname=rel,recursive=False)
                allhashes.append((rel,panel.digest(p)))
        archives.append((path.relative_to(dest).as_posix(),panel.digest(path)))
    for f, entries in ((dest/"results/raw_files.sha256",allhashes),(dest/"ARCHIVE.sha256",archives)):
        f.write_text("".join(f"{h}  {p}\n" for p,h in sorted(entries)))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path)
    args=parser.parse_args()
    pack(args.work.resolve(), Path(__file__).resolve().parent)
