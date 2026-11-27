"""Build the published research core, verify its numerical test, and rerun fixed first hunts."""
import argparse
import os
import shlex
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "pair_gradient"))
from build_arms import export_src, git

CORES = {"tangent": "4ba83b5087452639de42502bc0194b6915363c26",
         "placement": "6e2513f213e23db4b0da41ac5aa89945f882ac6b"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("arm", choices=CORES)
    ap.add_argument("work", type=Path, help="fresh output directory")
    ap.add_argument("--boost-prefix", type=Path)
    ap.add_argument("--jobs", type=int, default=2)
    args = ap.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True)
    source = work / "source"
    source.mkdir()
    rev = CORES[args.arm]
    local = HERE if args.arm == "tangent" else HERE / "placement"
    export_src(REPO, rev, source)
    (source / "tools").mkdir()
    builder = source / "tools/build.py"
    builder.write_bytes(git(REPO, "show", rev + ":tools/build.py"))
    subprocess.run(["git", "apply", str(local / "observe.patch")], cwd=source, check=True,
                   env={**os.environ, "GIT_CEILING_DIRECTORIES": str(work)})
    # Keep the actual compile/link flags from the normal builder for the numerical test.
    build = work / "build"
    command = [sys.executable, builder, "--output-dir", build, "--jobs", str(args.jobs)]
    if args.boost_prefix:
        command += ["--boost-prefix", args.boost_prefix.resolve()]
    log = subprocess.check_output(command, text=True)
    (work / "build.log").write_text(log)
    link = next(shlex.split(line) for line in log.splitlines() if line.endswith(" -o " + str(build / "vina")))
    unit = HERE / ("check_grid.cpp" if args.arm == "tangent" else "placement/check_placement.cpp")
    exclude = {"main.o"}
    if args.arm == "placement":
        exclude.add("monte_carlo.o")
    unit_link = [arg for arg in link[:-2] if Path(arg).name not in exclude]
    compiler = unit_link.pop(0)
    check = work / "check"
    subprocess.run([compiler, str(unit), *unit_link, "-o", check], check=True)
    # Restore provenance-qualified inputs from the older structural panel.
    fixtures = work / "inputs"
    reference_panel = HERE.parent / "search_backtracking/pose"
    with tarfile.open(reference_panel / "inputs.tar.xz") as archive:
        archive.extractall(fixtures, filter="data")
    with tarfile.open(reference_panel / "runs.tar.xz") as archive:
        for member in archive.getmembers():
            if member.name.startswith("crystal_diagnostic/") and member.name.endswith("-best_projected-input.pdbqt"):
                archive.extract(member, fixtures, filter="data")
    argv = [check]
    if args.arm == "placement":
        argv.append(fixtures / "4UYE/ligand.pdbqt")
    subprocess.run(argv, check=True)
    results = work / "results"
    subprocess.run([sys.executable, local / "first_hunt.py", fixtures, build / "vina", results], check=True)
    if (results / "rows.csv").read_bytes() != (local / "first_hunts.csv").read_bytes():
        raise RuntimeError("Numerical results differ from the archived table; inspect rows.csv and platform versions")
    print(args.arm, "matched the published per-case table")


if __name__ == "__main__":
    main()
