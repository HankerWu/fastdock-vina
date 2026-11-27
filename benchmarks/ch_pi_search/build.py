"""Build the isolated search prototype; optionally replay the numerical check."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--boost-prefix", type=Path)
    parser.add_argument("--jobs", default=4, type=int)
    parser.add_argument("--check-ligand", type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    with TemporaryDirectory() as tmp:
        base = Path(tmp)
        shutil.copytree(ROOT / "src", base / "src")
        (base / "tools").mkdir()
        build = base / "tools/build.py"
        shutil.copyfile(ROOT / "tools/build.py", build)
        shutil.copyfile(HERE / "search_grid.h", base / "src/lib/search_grid.h")
        subprocess.run(
            ["patch", "-p1", "-i", str(HERE / "search.patch")], cwd=base, check=True
        )
        cmd = [sys.executable, str(build), "--jobs", str(args.jobs)]
        if args.boost_prefix:
            cmd += ["--boost-prefix", str(args.boost_prefix.resolve())]
        subprocess.run(cmd + ["--output-dir", str(output)], check=True)
        if args.check_ligand:
            shutil.copyfile(HERE / "check.cpp", base / "src/main/main.cpp")
            check_dir = base / "check"
            subprocess.run(cmd + ["--output-dir", str(check_dir)], check=True)
            shutil.copyfile(check_dir / "vina", output / "check")
            (output / "check").chmod(0o755)
            subprocess.run(
                [str(output / "check"), str(args.check_ligand.resolve())], check=True
            )


if __name__ == "__main__":
    main()
