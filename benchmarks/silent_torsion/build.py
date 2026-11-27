"""Build one research MC object against a clean parent build from tools/build.py."""

import argparse
import shlex
import subprocess
import tempfile
from pathlib import Path

here = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "checkout", type=Path, help="source checkout of the protocol's parent"
)
parser.add_argument(
    "build_log", type=Path, help="stdout of tools/build.py for this exact checkout"
)
parser.add_argument("mode", choices=["always", "adaptive"])
parser.add_argument("output", type=Path)
parser.add_argument("--observer", action="store_true")
args = parser.parse_args()
source = args.checkout.resolve()
commands = [
    shlex.split(s) for s in args.build_log.read_text().splitlines() if s.strip()
]
compile_cmd = next(
    c
    for c in commands
    if "-c" in c and Path(c[c.index("-c") + 1]).name == "monte_carlo.cpp"
)
link_cmd = next(
    c for c in commands if "-o" in c and Path(c[c.index("-o") + 1]).name == "vina"
)
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    directory = work / "src/lib"
    directory.mkdir(parents=True)
    cpp = directory / "monte_carlo.cpp"
    cpp.write_bytes((source / "src/lib" / cpp.name).read_bytes())
    subprocess.run(
        ["patch", "-p1", "-i", str(here / "proposal.patch")], cwd=work, check=True
    )
    obj = work / "research.o"
    old = compile_cmd[compile_cmd.index("-o") + 1]
    compile_cmd[compile_cmd.index("-c") + 1] = str(cpp)
    compile_cmd[compile_cmd.index("-o") + 1] = str(obj)
    compile_cmd.extend(
        [
            "-I" + str(source / "src/lib"),
            "-DSILENT_TORSION_ALWAYS=" + str(int(args.mode == "always")),
        ]
    )
    if args.observer:
        compile_cmd.append("-DSILENT_TORSION_OBSERVER")
    subprocess.run(compile_cmd, check=True)
    link_cmd[link_cmd.index(old)] = str(obj)
    link_cmd[link_cmd.index("-o") + 1] = str(args.output.resolve())
    subprocess.run(link_cmd, check=True)
