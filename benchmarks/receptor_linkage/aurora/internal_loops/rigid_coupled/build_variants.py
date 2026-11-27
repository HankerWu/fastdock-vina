"""Compile the patched mutation unit against a clean build of the frozen source."""

import argparse
import os
import shlex
import subprocess
from pathlib import Path


def compile_variants(log, destination):
    # All compiler options, ABI and library order are identical to tools/build.py.
    commands = [
        shlex.split(line) for line in log.read_text().splitlines()
        if " -o " in line and not line.startswith("Built ")
    ]
    compile_cmd, = [cmd for cmd in commands if "-c" in cmd and cmd[-1].endswith("/mutate.o")]
    link_cmd, = [cmd for cmd in commands if "-c" not in cmd and cmd[-1].endswith("/vina")]
    source = Path(compile_cmd[compile_cmd.index("-c") + 1]).resolve()
    src = destination / "source"
    src.mkdir(parents=True, exist_ok=False)
    old = src / "src/lib/mutate.cpp"
    old.parent.mkdir(parents=True)
    old.write_bytes(source.read_bytes())
    patch = Path(__file__).with_name("implementation.patch").resolve()
    subprocess.run(["git", "apply", patch], cwd=src, check=True,
                   env=dict(os.environ, GIT_CEILING_DIRECTORIES=str(src.parent)))
    for arm in ("pivot", "sham"):
        out = destination / arm
        out.mkdir(exist_ok=True)
        obj = str(out / "mutate.o")
        cmd = [s if s != str(source) else str(old) for s in compile_cmd]
        cmd[-1] = obj
        if arm == "sham":
            cmd.insert(1, "-DFD_PIVOT_SHAM")
        subprocess.run(cmd, check=True)
        link = [obj if item == compile_cmd[-1] else item for item in link_cmd]
        link[-1] = str(out / "vina")
        subprocess.run(link, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_log", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    compile_variants(args.base_log.resolve(), args.destination.resolve())
