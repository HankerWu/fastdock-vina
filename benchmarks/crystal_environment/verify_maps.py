"""Reconstruct deposited maps with explicit reciprocal symmetry and a NumPy FFT."""

import argparse
from pathlib import Path

import gemmi
import numpy as np


def verify(path):
    block = gemmi.as_refln_blocks(gemmi.cif.read(str(path)))[0]
    reference = block.transform_f_phi_to_map("pdbx_FWT", "pdbx_PHWT", sample_rate=3)
    reference.normalize()
    miller = block.make_miller_array()
    amplitude = block.make_float_array("pdbx_FWT")
    phase = block.make_float_array("pdbx_PHWT")
    present = np.isfinite(amplitude) & np.isfinite(phase)
    hkl = miller[present]
    factor = amplitude[present] * np.exp(1j * np.deg2rad(phase[present]))
    coefficients = np.zeros(reference.shape, dtype=np.complex64)
    for op in block.spacegroup.operations():
        rotation = np.array(op.rot) / op.DEN
        translation = np.array(op.tran) / op.DEN
        related = np.rint(hkl @ np.linalg.inv(rotation)).astype(int)
        value = factor * np.exp(2j * np.pi * (related @ translation))
        coefficients[tuple(related.T)] = value
        coefficients[tuple(-related.T)] = value.conjugate()
    grid = np.fft.fftn(coefficients).real
    grid -= grid.mean()
    grid /= grid.std()
    diff = grid - reference.array
    maximum = np.max(np.abs(diff))
    if maximum > 3e-5:
        raise ValueError((path, maximum))
    print(path.name, reference.shape, "max", maximum, "RMS", np.sqrt(np.mean(diff**2)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    paths = sorted(
        parser.parse_args().directory.glob("*_validation_2fo-fc_map_coef.cif.gz")
    )
    if len(paths) != 11:
        raise ValueError(f"Expected eleven experimental map files, found {len(paths)}")
    for path in paths:
        verify(path)
