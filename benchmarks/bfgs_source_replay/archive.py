"""Pack diagnostic logs losslessly, restore reconstructed controls and verify original hashes."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile

from run_replay import ARMS, CASES, digest, run_panel


def pack(work, dest):
    dest.mkdir(exist_ok=True, parents=True)
    members = {}
    packed = io.BytesIO()
    with tarfile.open(fileobj=packed, mode="w:xz", preset=3) as archive:
        for query, receiver, seed in sorted(CASES):
            path = work / f"{query}-{receiver}-{seed}"
            for source in sorted(path.iterdir()):
                if source.suffix == ".inputs" or source.name in {f"{a}-{a}.shadow" for a in ARMS}:
                    continue  # Recover these exactly from native proposal/hunt records.
                key = source.relative_to(work).as_posix()
                info = tarfile.TarInfo(key)
                info.size = source.stat().st_size
                with source.open("rb") as stream:
                    archive.addfile(info, stream)
                members[key] = digest(source)
    data = packed.getvalue()
    parts = []
    for i in range(0, len(data), 750000):
        blob = dest / f"trace.tar.xz.{len(parts):02d}"
        blob.write_bytes(data[i:i+750000])
        parts.append(dict(name=blob.name, sha256=digest(blob)))
    for name in ("costs.json", "candidate_costs.json"):
        shutil.copyfile(work/name, dest/name)
    manifest = dict(parts=parts, archive_sha256=hashlib.sha256(data).hexdigest(), members=members)
    (dest/"files.json").write_text(json.dumps(manifest, indent=2)+"\n")


def restore(src, dest):
    dest.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((src/"files.json").read_text())
    data = io.BytesIO()
    for part in manifest["parts"]:
        blob = src/part["name"]
        assert digest(blob) == part["sha256"]
        data.write(blob.read_bytes())
    assert hashlib.sha256(data.getvalue()).hexdigest() == manifest["archive_sha256"]
    data.seek(0)
    with tarfile.open(fileobj=data, mode="r:xz") as archive:
        assert sorted(archive.getnames()) == sorted(manifest["members"])
        archive.extractall(dest, filter="data")
    assert all(digest(dest/name) == checksum for name,checksum in manifest["members"].items())
    for name in ("costs.json","candidate_costs.json"):
        shutil.copyfile(src/name,dest/name)
    for path in dest.iterdir():
        if not path.is_dir():
            continue
        for arm in ARMS:
            lines=[l.split("\t") for l in (path/(arm+".native")).read_text().splitlines()]
            props=[f[3] for f in lines if f[0]=="P"]
            native=[f for f in lines if f[0]=="H"]
            (path/(arm+".inputs")).write_text("\n".join(props)+"\n")
            (path/f"{arm}-{arm}.shadow").write_text("".join("\t".join([r[0],str(i),*r[2:]])+"\n"
                                                         for i,r in enumerate(native)))
    for name in ("costs.json","candidate_costs.json"):
        for row in json.loads((dest/name).read_text()):
            key=(row["driver"]+".native" if row["kind"]=="native" else
                 row["driver"]+"-"+row["evaluator"]+".shadow")
            assert digest(dest/row["case"]/key)==row["trace_sha256"]
    # Use the same immutable panel extraction/checks as a full run.
    list(run_panel.panel(dest/"fixtures"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=("pack","restore"))
    p.add_argument("source",type=Path)
    p.add_argument("dest",type=Path)
    a=p.parse_args()
    {"pack":pack,"restore":restore}[a.mode](a.source,a.dest)
