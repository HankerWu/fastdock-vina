"""Preserve every observed coordinate, log, phase counter and child-process receipt."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

def sha(data):
    return hashlib.sha256(data).hexdigest()

def pack(work, raw):
    raw.mkdir()
    members = {}
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream,mode="w:xz",preset=4) as archive:
        for p in sorted(work.rglob("*")):
            key = p.relative_to(work).as_posix()
            if not p.is_file() or key.startswith("fixtures/"):
                continue
            data = p.read_bytes()
            members[key] = sha(data)
            info = tarfile.TarInfo(key)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    data = stream.getvalue()
    chunks = {}
    for n,offset in enumerate(range(0,len(data),720000)):
        part = data[offset:offset+720000]
        name = f"stages.tar.xz.{n:02}"
        (raw/name).write_bytes(part)
        chunks[name]=sha(part)
    (raw/"manifest.json").write_text(json.dumps(dict(members=members,chunks=chunks,sha256=sha(data)),indent=2)+"\n")

def restore(raw, work):
    manifest=json.loads((raw/"manifest.json").read_text())
    chunks=[]
    for name,expected in manifest["chunks"].items():
        part=(raw/name).read_bytes()
        assert sha(part)==expected
        chunks.append(part)
    data=b"".join(chunks)
    assert sha(data)==manifest["sha256"]
    work.mkdir()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:xz") as archive:
        assert sorted(archive.getnames())==sorted(manifest["members"])
        archive.extractall(work,filter="data")
    for name,expected in manifest["members"].items():
        assert sha((work/name).read_bytes())==expected

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("method",choices=("pack","restore"))
    p.add_argument("source",type=Path)
    p.add_argument("dest",type=Path)
    a=p.parse_args()
    {"pack":pack,"restore":restore}[a.method](a.source,a.dest)
