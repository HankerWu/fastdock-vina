"""Convert existing cavity transforms into native ligand-root translations and quaternions."""

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation


def write_bank(path, x, q, centers, rotations):
    """Preserve candidate order. Vina quaternions use the scalar-first convention."""
    if rotations.shape != (1024, 3, 3) or centers.shape != (1024, 3):
        raise ValueError("expected the frozen bank of 1024 rigid transforms")
    # The original frozen banks composed one orientation at a time. Batching
    # rigid operations can change a translation's last bit and an entire MC path.
    poses = np.empty((len(centers), 7))
    centered_root = q[:3] - x.mean(axis=0)
    source_orientation = Rotation.from_quat(q[[4, 5, 6, 3]])
    for i, (matrix, center) in enumerate(zip(rotations, centers, strict=True)):
        rotation = Rotation.from_matrix(matrix)
        poses[i, :3] = rotation.apply(centered_root) + center
        poses[i, 3:] = (rotation * source_orientation).as_quat()[[3, 0, 1, 2]]
    if not np.isfinite(poses).all():
        raise ValueError("nonfinite rigid transform")
    np.savetxt(path, poses, fmt="%.17g")
    return poses


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("candidates", type=Path)
    parser.add_argument("bank", type=Path)
    parser.add_argument("--geometry", choices=["matched", "sham"], default="matched")
    args = parser.parse_args()
    with np.load(args.source) as source, np.load(args.candidates) as candidates:
        write_bank(args.bank, source["xyz"], source["q"], candidates["center"],
                   candidates[args.geometry + "_rotation"])
