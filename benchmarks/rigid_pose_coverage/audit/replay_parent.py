"""Verify archives and repeat all four specified parent analyses in a clean directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import time

D = Path(__file__).resolve().parents[2] / 'fragment_anchor'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(manifest, root):
    lines = manifest.read_text().splitlines()
    for line in lines:
        sha, path = line.split('  ', 1)
        if digest(root/path) != sha:
            raise ValueError(path)
    return len(lines)


def main(work, observer, receipt):
    pose, out = work/'parent', work/'out'
    pose.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=True)
    n_archive = verify(D/'pose/ARCHIVE.sha256', D/'pose')
    for path in sorted((D/'pose/archives').glob('*.tar.xz')):
        with tarfile.open(path) as archive:
            archive.extractall(pose, filter='data')
    n_files = verify(D/'pose/results/raw_files.sha256', pose)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    start = time.perf_counter()
    costs = {}
    for script, args, label in [
            ('edge_audit.py', ['pose',work/'frames',out,'--pose',pose], 'pose'),
            ('edge_audit.py', ['at9283',work/'at',out], 'at9283'),
            ('omitted.py', [pose,out], 'omitted'),
            ('forks.py', [pose,work/'forks',observer,out], 'forks')]:
        t = time.perf_counter()
        subprocess.run([sys.executable,D/'audit'/script,*args],env=env,check=True)
        costs[label] = time.perf_counter()-t
    checks = {}
    for name in ('gates_pose','gates_at9283','omitted','near_forks'):
        path = out/(name+'.csv')
        if path.read_bytes() != (D/'audit/results'/path.name).read_bytes():
            raise ValueError(f'{name} differs from parent archive')
        checks[name] = dict(sha256=digest(path),records=len(path.read_text().splitlines())-1)
    receipt.write_text(json.dumps(dict(archives=n_archive,raw_files=n_files,observer_sha256=digest(observer),
        matched=checks,analysis_wall_s=time.perf_counter()-start,stage_wall_s=costs),indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('work','observer','receipt'):
        parser.add_argument(name,type=lambda x:Path(x).resolve())
    args = parser.parse_args()
    main(args.work,args.observer,args.receipt)
