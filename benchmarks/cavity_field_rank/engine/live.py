"""Generate and screen the bank in one live, single-CPU XS docking run (POSIX)."""

import argparse
import errno
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


def open_bank(stream, child):
    # The observer closes both source snapshots before opening its read end.
    # Check for the child stopping before the selected step, including errors.
    while child.poll() is None:
        try:
            fd = os.open(stream, os.O_WRONLY | os.O_NONBLOCK)
        except OSError as exc:
            if exc.errno != errno.ENXIO:
                raise
            time.sleep(0.002)
            continue
        os.set_blocking(fd, True)
        return os.fdopen(fd, "w")
    raise RuntimeError(f"engine exited before the requested rank event: {child.returncode}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", required=True)
    parser.add_argument("--arm", required=True, choices=["true","perm"])
    parser.add_argument("--trajectory", type=int, default=1)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("native", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if min(args.trajectory, args.step) < 1:
        parser.error("trajectory and step start at one")
    native = args.native[1:] if args.native[:1] == ["--"] else args.native
    options = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    options.add_argument("--receptor", required=True, type=Path)
    options.add_argument("--seed", required=True, type=int)
    options.add_argument("--cpu", required=True, choices=[1], type=int)
    options.add_argument("--out", required=True, type=Path)
    for axis in "xyz":
        options.add_argument("--size_"+axis, required=True, type=float)
        options.add_argument("--center_"+axis, required=True, type=float)
    cfg, _ = options.parse_known_args(native)
    with tempfile.TemporaryDirectory(dir=cfg.out.parent) as tmp:
        work = Path(tmp).resolve()
        pipe, log, spec = work/"bank.fifo", work/"selection.log", work/"event"
        q, xyz, pool = work/"start.q", work/"start.xyz", work/"pool.pdbqt"
        spec.write_text(f"{args.trajectory} {args.step} 0 0 0 1 0 0 0\n")
        os.mkfifo(pipe, mode=0o600)
        env = {k:v for k,v in os.environ.items() if not k.startswith("FD_")}
        env.update(FD_FIELD_BANK=str(pipe), FD_FIELD_ARM=args.arm, FD_RIGID_LOG=str(log),
                   FD_RIGID_INTERVENTION=str(spec), FD_SOURCE_Q=str(q), FD_SOURCE_XYZ=str(xyz),
                   FD_SOURCE_STEP=str(args.step), FD_SOURCE_TRAJECTORY=str(args.trajectory))
        forwarded = []
        skip = False
        for word in native:
            if skip:
                skip = False
            elif word == "--out":
                skip = True
            elif not word.startswith("--out="):
                forwarded.append(word)
        child = subprocess.Popen([args.engine, *forwarded, "--out", str(pool)], env=env)
        try:
            import numpy as np
            from prepare_bank import write_bank
            sys.path.insert(0, str(Path(__file__).resolve().parents[3]/"tools"))
            from pocket_triangle import Pocket
            with open_bank(pipe, child) as dest:
                header, coords = xyz.read_text().splitlines()
                serials = np.array(header.split()[1:],int)
                if len(np.unique(serials)) != len(serials):
                    raise ValueError("duplicate ligand atom serial")
                x = np.array(coords.split()[3:],float).reshape(len(serials),3)[np.argsort(serials)]
                seed_conf = np.loadtxt(q)[2:]
                center = [getattr(cfg,"center_"+a) for a in "xyz"]
                size = [getattr(cfg,"size_"+a) for a in "xyz"]
                pocket = Pocket.from_pdbqt(cfg.receptor, center, size)
                bank = pocket.place(x,cfg.seed,arms=("matched",))
                write_bank(dest, x, seed_conf, bank["center"], bank["matched_rotation"])
            code = child.wait()
            if code:
                raise subprocess.CalledProcessError(code, child.args)
            summary = [row for row in log.read_text().splitlines() if row.startswith("field_rank ")]
            if len(summary) != 1:
                raise ValueError("did not observe exactly one field rank event")
            print(summary[0], flush=True)
            print(f"site_nearest {pocket.n_grid} burial {pocket.n_burial} "
                  f"shape_queries {len(x)*len(bank['center'])}",flush=True)
            os.replace(pool, cfg.out)
        finally:
            if child.poll() is None:
                child.terminate()
            child.wait()


if __name__ == "__main__":
    main()
