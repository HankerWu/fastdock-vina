"""Check separation at three native selection stages in the archived 66-run panel."""

import csv
import io
import itertools
from pathlib import Path
import sys
import tarfile

import numpy as np


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / 'accepted_retention/coverage/results'
STAGES = ((6, 'mc_exit', 1.0), (7, 'merged', 2.0), (8, 'explicit_input', 1.0))
# Every stored Cartesian component was rounded to 0.001 Å. For a pair of
# structures the norm of the per-atom displacement error is at most sqrt(3)*0.001 Å.
DISTANCE_ERROR = np.sqrt(3.0) * 0.001


def inspect_panel():
    for family in ('DYRK1A', 'EPAS1', 'FABP4', 'FKBP1A', 'JAK2', 'RBP4'):
        with tarfile.open(SOURCE / f'{family}.tar.xz') as archive:
            for entry in archive:
                if not entry.name.endswith('.npz'):
                    continue
                with np.load(io.BytesIO(archive.extractfile(entry).read())) as data:
                    stage_ids = data['stage']
                    xyz = data['xyz_milli'].astype(float) / 1000
                for stage_id, name, radius in STAGES:
                    poses = xyz[stage_ids == stage_id]
                    pairs = [(np.sqrt(np.square(a - b).sum(axis=1).mean()), i + 1, j + 1)
                             for (i, a), (j, b) in itertools.combinations(enumerate(poses), 2)]
                    minimum = min(pairs) if pairs else None
                    yield (entry.name[:-4], name, len(poses), radius,
                           sum(d + DISTANCE_ERROR < radius for d, _, _ in pairs),
                           sum(abs(d - radius) <= DISTANCE_ERROR for d, _, _ in pairs),
                           '' if minimum is None else f'{minimum[0]:.6f}',
                           '' if minimum is None else minimum[1], '' if minimum is None else minimum[2])


if __name__ == '__main__':
    output = csv.writer(sys.stdout, lineterminator='\n')
    output.writerow(('trajectory', 'stage', 'poses', 'radius_A', 'definite_conflicts',
                     'rounding_ambiguous_pairs', 'closest_rmsd_A', 'closest_rank_a', 'closest_rank_b'))
    output.writerows(inspect_panel())
