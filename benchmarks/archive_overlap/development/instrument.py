"""Add passive archive stream recording to a disposable source checkout."""

import argparse
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent


def instrument(source):
    lib = source / "src/lib"
    shutil.copyfile(HERE / "observer.h", lib / "observer.h")
    for name, stage in (
        ("monte_carlo", "save"),
        ("parallel_mc", "merge"),
        ("vina", "prefine"),
    ):
        path = lib / (name + ".cpp")
        body = path.read_text()
        token = "add_to_output_container("
        if body.count(token) != 1:
            raise ValueError(f"Unexpected candidate insertion sites in {path}")
        body = body.replace(token, f'observed_insert("{stage}", ')
        if name == "monte_carlo":
            body = body.replace(
                "int evalcount = 0;", "int evalcount = 0;\n    trace_atom_order(m);"
            )
        path.write_text('#include "observer.h"\n' + body)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkout", type=Path)
    instrument(parser.parse_args().checkout)
