"""Physical nulls on actual frozen ligand coordinates with production topology."""

import argparse
import csv
import io
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from run import nodes


def compute(binary, rec, lig, pool, grid):
    text = subprocess.check_output([binary, rec, lig, pool, grid, "298.15"], text=True)
    return {k: float(v) for k, v in list(csv.DictReader(io.StringIO(text)))[-1].items()}


def run(args):
    p = args.run / "galectin_affinity"
    rec = p / "prepared/1GZW.pdbqt"
    lig = p / "prepared/10.pdbqt"
    pool = p / "pools/1GZW_10_141421_native.pdbqt"
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        out = root / "grid"
        nodes(out, 256)
        raw = np.loadtxt(out)
        assert np.max(np.linalg.norm(raw[:, :3], axis=1)) < 1
        rotation = raw[:, 3:12].reshape(-1, 3, 3)
        np.testing.assert_allclose(
            rotation @ rotation.transpose(0, 2, 1),
            np.broadcast_to(np.eye(3), rotation.shape),
            atol=1e-12,
        )
        assert np.min(np.linalg.det(rotation)) > 0.999999999
        assert (
            np.min(rotation.trace(axis1=1, axis2=2))
            >= 1 + 2 * np.cos(np.deg2rad(15)) - 1e-12
        )
        np.testing.assert_equal(
            np.bincount((16 * raw[:, 12]).astype(int)), np.full(16, 16)
        )
        # At infinite protein separation, all internal deformation energies
        # survive in both integrals. Matched reference must cancel them exactly.
        atom = next(l for l in rec.read_text().splitlines() if l.startswith("ATOM"))
        far = root / "far.pdbqt"
        far.write_text(atom[:30] + " 800.000" * 3 + atom[54:] + "\n")
        a = compute(args.integrator, far, lig, pool, out)
        assert a["hinge_moment"] > 0 and abs(a["free"]) > 1e-6
        np.testing.assert_allclose(
            [a["inter"], a["rigid"], a["joint"], a["hinge"]], 0, atol=1e-9
        )
        raw[:, :3] = 0
        raw[:, 3:12] = np.eye(3).ravel()
        np.savetxt(out, raw, fmt="%.16g")
        b = compute(args.integrator, rec, lig, pool, out)
        np.testing.assert_allclose(b["rigid"], 0, atol=1e-9)
        np.testing.assert_allclose(b["joint"], b["hinge"], atol=1e-9)
        assert abs(b["hinge"]) > 1e-3
        # A macrocycle put through an open-cut torsion would change its closure
        # distance. The present bounded model must give the rigid limit there.
        k = args.run / "torsion_kinase_affinity"
        nodes(out, 256)
        c = compute(
            args.integrator,
            k / "prepared/2E9P.pdbqt",
            k / "prepared/153999.pdbqt",
            k / "pools/2E9P_153999_60913_native.pdbqt",
            out,
        )
        np.testing.assert_allclose(
            [c["hinge_moment"], c["hinge"], c["free"]], 0, atol=1e-9
        )
        np.testing.assert_allclose(c["joint"], c["rigid"], atol=1e-9)
        print(
            "SO(3), matched free reference, zero rigid perturbation, macrocycle fallback: passed"
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("integrator", type=Path)
    p.add_argument("run", type=Path)
    run(p.parse_args())
