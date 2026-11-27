"""One-pass live placement. Forward ordinary explicit Vina arguments after -- (single CPU)."""

import argparse
import json
import os
from pathlib import Path
import socket
import subprocess
import tempfile
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", required=True)
    parser.add_argument("--arm", choices=["matched", "sham"], required=True)
    parser.add_argument("--trajectory", type=int, default=1)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("vina_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    native = args.vina_args[1:] if args.vina_args[:1] == ["--"] else args.vina_args
    needed = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    needed.add_argument("--receptor", required=True, type=Path)
    needed.add_argument("--cpu", required=True, type=int, choices=[1])
    needed.add_argument("--seed", required=True, type=int)
    needed.add_argument("--out", required=True, type=Path)
    for a in "xyz":
        needed.add_argument("--center_" + a, required=True, type=float)
        needed.add_argument("--size_" + a, required=True, type=float)
    cfg, _ = needed.parse_known_args(native)
    if args.trajectory < 1 or args.step < 1:
        parser.error("trajectory and step are counted from 1")
    center = [getattr(cfg, "center_" + a) for a in "xyz"]
    size = [getattr(cfg, "size_" + a) for a in "xyz"]
    # Vina may report a complete unmodified pool when the selected event never
    # occurs. Publish only after the handshake and an ordinary successful exit.
    cleaned = []
    index = 0
    while index < len(native):
        token = native[index]
        if token == "--out":
            index += 2
        elif token.startswith("--out="):
            index += 1
        else:
            cleaned.append(token)
            index += 1
    parent, child = socket.socketpair()
    with parent, child, tempfile.TemporaryDirectory(dir=cfg.out.parent) as tmp:
        pending = Path(tmp) / "pool.pdbqt"
        native = [*cleaned, "--out", str(pending)]
        fd = child.fileno()
        env = os.environ | {"FD_TRIPLET_SOCKET": str(fd), "FD_TRIPLET_STEP": str(args.step),
                            "FD_TRIPLET_TRAJECTORY": str(args.trajectory)}
        engine = subprocess.Popen([args.engine, *native], pass_fds=(fd,), env=env)
        child.close()
        try:
            with parent.makefile("rwb") as stream:
                # EOF means the engine exited without requesting the selected event.
                line = stream.readline()
                if not line:
                    raise RuntimeError(f"engine exited before the live event (exit={engine.wait()})")
                import numpy as np
                from scipy.spatial.transform import Rotation
                sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tools"))
                from pocket_triangle import Pocket
                q = np.fromstring(line.decode(), sep=" ")
                n = int(stream.readline())
                x = np.array([list(map(float, stream.readline().split())) for _ in range(n)])
                if q.shape != (7,) or x.shape != (n,3):
                    raise ValueError("incomplete native coordinate stream")
                pocket = Pocket.from_pdbqt(cfg.receptor, center, size)
                placement = pocket.place(x,cfg.seed,arms=(args.arm,))
                i = placement[args.arm + "_best"]
                rotation = Rotation.from_matrix(placement[args.arm + "_rotation"][i])
                origin = rotation.apply(q[:3] - x.mean(axis=0)) + placement["center"][i]
                quat = rotation * Rotation.from_quat(q[[4,5,6,3]])
                out = np.r_[origin, quat.as_quat()[[3,0,1,2]]]
                stream.write((" ".join(format(v,".17g") for v in out) + "\n").encode())
                stream.flush()
                print("Live cavity placement " + json.dumps(dict(
                    arm=args.arm, selected=i, candidates=len(placement["center"]),
                    atom_queries=len(placement["center"])*len(x),
                    site_nearest_queries=pocket.n_grid, site_count_queries=pocket.n_burial,
                    trajectory=args.trajectory,step=args.step)), flush=True)
            # Closing before wait also lets a duplicate request fail cleanly.
            parent.close()
            status = engine.wait()
            if status:
                raise subprocess.CalledProcessError(status,[args.engine,*native])
            os.replace(pending, cfg.out)
        except BaseException:
            parent.close()
            if engine.poll() is None:
                engine.terminate()
            engine.wait()
            raise


if __name__ == "__main__":
    main()
