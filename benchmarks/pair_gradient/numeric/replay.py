"""Check table forces and full-model derivatives with each prebuilt docking arm."""
import argparse
import csv
import gzip
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
FIXTURES = HERE.parents[1] / "torsion_proposals" / "fixtures"


def read_probe(args):
    output = subprocess.check_output(args, text=True)
    return list(csv.DictReader(io.StringIO(output)))


def run(args):
    pair_rows, model_rows = [], []
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for filename in ("inputs.tar.xz", "crystal_diagnostic.tar.xz"):
            with tarfile.open(FIXTURES / filename) as archive:
                archive.extractall(work, filter="data")
        for name in ("native", "linear", "hermite", "shifted"):
            root, build = getattr(args, name)
            obj = build / "obj"
            objects = sorted(p for p in obj.glob("*.o") if p.stem not in {"main", "split"})
            if not objects:
                raise RuntimeError(f"No compiled objects in {obj}")
            flags = ["-O3", "-DNDEBUG", "-std=c++14", "-pthread", f"-I{root / 'src/lib'}"]
            if name == "native":
                flags.append("-DPAIR_AUDIT_NATIVE")
            elif name == "hermite":
                flags.append("-DPAIR_AUDIT_MONOTONE")
            else:
                flags.append("-DFASTDOCK_EXACT_PAIR_GRADIENT")
            link = [f"-lboost_{s}" for s in ("thread", "serialization", "filesystem", "program_options")]
            if args.boost_prefix:
                flags.append(f"-I{args.boost_prefix / 'include'}")
                lib = args.boost_prefix / "lib"
                link[:0] = [f"-L{lib}", f"-Wl,-rpath,{lib}"]
            probes = ("shift", "model") if name == "shifted" else ("pair", "model")
            for probe in probes:
                subprocess.run([args.cxx, *flags, HERE / f"{probe}.cpp", *objects, *link,
                                "-o", work / probe], check=True)
            if name == "shifted":
                with (args.output / "shift.csv").open("w") as out:
                    subprocess.run([work / "shift"], stdout=out, check=True)
            else:
                pair_rows.extend(dict(arm=name, **row) for row in read_probe([work / "pair"]))
            if name == "native":
                trace = read_probe([work / "model", work / "2YWP/receptor.pdbqt",
                                    work / "crystal_diagnostic/2YWP-2YWP-crystal-input.pdbqt",
                                    "--trace-y-torque"])
                with (args.output / "2ywp_torque.csv.gz").open("wb") as raw:
                    with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
                        with io.TextIOWrapper(gz, encoding="utf-8", newline="") as f:
                            writer = csv.DictWriter(f, list(trace[0]))
                            writer.writeheader()
                            writer.writerows(trace)
            for ligand in sorted((work / "crystal_diagnostic").glob("*-crystal-input.pdbqt")):
                source, receptor, *_ = ligand.name.split("-")
                for row in read_probe([work / "model", work / receptor / "receptor.pdbqt", ligand]):
                    model_rows.append(dict(arm=name, ligand=source, receptor=receptor, **row))
    with (args.output / "pair.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)
    path = args.output / "model.csv.gz"
    with path.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        with io.TextIOWrapper(gz, encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, list(model_rows[0]))
            writer.writeheader()
            writer.writerows(model_rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    for arm in ("native", "linear", "hermite", "shifted"):
        p.add_argument("--" + arm, nargs=2, type=Path, required=True,
                       metavar=("CHECKOUT", "BUILD"), help="Matching headers and terminal build/obj directory")
    p.add_argument("--boost-prefix", type=Path)
    p.add_argument("--output", required=True, type=Path)
    run(p.parse_args())
