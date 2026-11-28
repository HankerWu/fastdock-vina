"""Combine four external result folders using the frozen panel's summaries."""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'geometry_comparison/upstream/four_versions'))
from evaluate import collect, summaries


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    collect(args.work, args.out)
    summaries(args.out)
