"""Publish the additions to the existing hash-locked candidate/pool archives."""
import argparse
import json
from pathlib import Path
import shutil
import tarfile

from frozen_assay import ARMS, sha


def archive(work, out):
    out.mkdir(parents=True, exist_ok=True)
    reports=out/"results"
    reports.mkdir(exist_ok=True)
    for path in (work/"results").iterdir():
        shutil.copyfile(path,reports/path.name)
    for name in ("rank_cost.csv","source_hashes.json","official_high_replay.csv","merge_cost.csv","end_to_end.csv","ordinary-typed.json","ordinary-perm.json"):
        shutil.copyfile(work/name, reports/name)
    selected=[work/"stage_audit.csv",work/"candidate_audit.csv",work/"rank_cost.csv",
              work/"official_high_replay.csv",work/"merge_cost.csv",work/"end_to_end.csv",work/"source_hashes.json"]
    for arm in (*ARMS,"official_high"):
        selected.extend((work/"pools"/arm).iterdir())
        selected.append(work/f"runs-{arm}.csv")
    for subdir in ("merge","ranks","failed_parse"):
        selected.extend((work/subdir).iterdir())
    hashes={str(p.relative_to(work)):sha(p) for p in sorted(selected)}
    (reports/"raw_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    path=out/"field_forks.tar.xz"
    with tarfile.open(path,"w:xz",preset=3) as tar:
        for name in hashes:
            tar.add(work/name,arcname=name,recursive=False)
    (out/"ARCHIVE.sha256").write_text(f"{sha(path)}  field_forks.tar.xz\n")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("work",type=Path)
    p.add_argument("out",type=Path)
    a=p.parse_args()
    archive(a.work.resolve(), a.out.resolve())
