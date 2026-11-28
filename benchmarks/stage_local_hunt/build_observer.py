"""Build the source-pinned hunt observer with the original optimizer call."""

import argparse
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bfgs_source_replay"))
from instrument_candidate import edit, instrument

REVISION = "0bcae6ef09857f97656291bb4821d20b3b619bf5"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("boost", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    with subprocess.Popen(
        ["git", "-C", HERE.parents[1], "archive", REVISION, "src", "tools/build.py"],
        stdout=subprocess.PIPE,
    ) as command:
        with tarfile.open(fileobj=command.stdout, mode="r|") as archive:
            archive.extractall(out, filter="data")
        if command.wait():
            raise RuntimeError("candidate source export failed")
    instrument(out)
    # Use exactly the optional solver arguments of the public MC hunt entry.
    edit(
        out / "src/lib/bfgs_trace.h",
        "local(m, p, ig, pose, cap, cost);",
        "local(m, p, ig, pose, cap, cost, nullptr, 0, true);",
    )
    subprocess.run(
        [
            sys.executable,
            out / "tools/build.py",
            "--jobs",
            "3",
            "--boost-prefix",
            args.boost.resolve(),
            "--output-dir",
            out / "build",
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
