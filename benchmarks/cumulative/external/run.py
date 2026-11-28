"""Evaluate a frozen external set using the same four-version pose runner."""
import argparse
import json
from pathlib import Path
import sys
import tarfile

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'geometry_comparison/upstream/four_versions'))
import run_panel as four_pose


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('arm', choices=four_pose.ARMS)
    p.add_argument('binary', type=Path)
    p.add_argument('destination', type=Path)
    a = p.parse_args()
    work, binary = a.destination.resolve(), a.binary.resolve(strict=True)
    raw, fixtures = work / 'raw', work / 'fixtures'
    raw.mkdir(parents=True, exist_ok=False)
    with tarfile.open(HERE / 'inputs.tar.xz') as archive:
        archive.extractall(fixtures, filter='data')
    frozen = (HERE / 'preparation.json').read_bytes()
    if (fixtures / 'manifest.json').read_bytes() != frozen:
        raise ValueError('changed input manifest')
    cells, ranks = [], []
    for entry in json.loads(frozen):
        case = fixtures / entry['pdb']
        for name, sha in entry['sha256'].items():
            if four_pose.digest(case / name) != sha:
                raise ValueError(f'changed fixture {case.name}/{name}')
        if entry['status'] != 'ready':
            continue
        meta = json.loads((case / 'metadata.json').read_text())
        for seed in four_pose.replay.SEEDS:
            job = ('external', fixtures, meta, meta, np.asarray(meta['mapped']), seed, 12000)
            cell, modes = four_pose.dock((a.arm, binary, raw, job))
            cells.append(cell)
            ranks.extend(modes)
    four_pose.csvfile(work / 'cells.csv', cells)
    if ranks:
        four_pose.csvfile(work / 'ranks.csv', ranks)
    (work / 'binary.sha256').write_text(four_pose.digest(binary) + '\n')


if __name__ == '__main__':
    main()
