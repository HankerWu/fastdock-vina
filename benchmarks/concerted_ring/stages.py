"""Observe pre/post-hunt docking states and ring activation without changing normal binaries."""

import argparse
import csv
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
MACRO = HERE.parent / "macrocycle_proposal"
sys.path.insert(0, str(MACRO))
sys.path.insert(0, str(MACRO / "external"))
from panel import command, counts, sha, write
from external_score import main as score


def observe(source, objects, output, boost):
    output.mkdir(parents=True, exist_ok=True)
    lib = source / "src/lib"
    mc = (lib / "monte_carlo.cpp").read_text()
    mc = mc.replace(
        '#include "monte_carlo.h"',
        """#include "monte_carlo.h"
#include <fstream>
static void trace(const model& m, const conf& c, const char* stage, sz step) {
    model scratch = m; // Preserve the real model's cached coordinates.
    scratch.set(c);
    const char* base = std::getenv("RING_OBSERVER");
    std::ofstream out(std::string(base)+"."+stage+".pdbqt",std::ios::app);
    scratch.write_model(out,step+1,"");
}
""",
    )
    for anchor, extra in (
        ("candidate.c = tmp.c;", 'trace(m,candidate.c,"before",step);'),
        (
            "const int before_hunt = evalcount;",
            'std::fprintf(stderr,"Ring step: %zu\\n",sz(step)); trace(m,candidate.c,"proposal",step);',
        ),
        (
            "quasi_newton_par(m, p, ig, candidate, hunt_cap, evalcount);",
            'trace(m,candidate.c,"hunt",step);',
        ),
        ("tmp.e = candidate.e;", 'trace(m,tmp.c,"accepted",step);'),
        (
            "saved.coords = m.get_heavy_atom_movable_coords();",
            'trace(m,saved.c,"second",step);',
        ),
    ):
        assert mc.count(anchor) == 1
        mc = mc.replace(anchor, anchor + "\n" + extra)
    (output / "monte_carlo.cpp").write_text(mc)
    original = (lib / "mutate.cpp").read_text()
    if "concerted_ring::propose" in original:
        anchor = "concerted_ring::propose(c, m, i, which, angle, amplitude);"
        assert original.count(anchor) == 1
        original = original.replace(
            anchor,
            """bool success = concerted_ring::propose(c, m, i, which, angle, amplitude);
                std::fprintf(stderr,"Ring proposal: axis=%zu draw=%.17g applied=%d\\n",which,angle,int(success));""",
        )
    else:
        anchor = "selected_entity))\n\t\t\t\treturn false;"
        assert original.count(anchor) == 1
        original = original.replace(
            anchor,
            """selected_entity)) {
                std::fprintf(stderr,"Ring proposal: axis=%zu draw=%.17g applied=0\\n",which,angle);
                return false;
            }""",
        )
    (output / "mutate.cpp").write_text("#include <cstdio>\n" + original)
    flags = [
        "c++",
        "-O3",
        "-DNDEBUG",
        "-std=c++14",
        "-pthread",
        "-fvisibility=hidden",
        "-I" + str(lib),
    ]
    libs = []
    if boost:
        flags += ["-I" + str(boost / "include")]
        libs += ["-L" + str(boost / "lib"), "-Wl,-rpath," + str(boost / "lib")]
    for unit in ["monte_carlo", "mutate"]:
        subprocess.run(
            [
                *flags,
                "-c",
                str(output / (unit + ".cpp")),
                "-o",
                str(output / (unit + ".o")),
            ],
            check=True,
        )
    unchanged = sorted(
        o
        for o in objects.glob("*.o")
        if o.stem not in {"monte_carlo", "mutate", "split"}
    )
    subprocess.run(
        [
            *flags,
            *unchanged,
            output / "monte_carlo.o",
            output / "mutate.o",
            *libs,
            "-lboost_filesystem",
            "-lboost_serialization",
            "-lboost_program_options",
            "-lboost_thread",
            "-o",
            output / "vina",
        ],
        check=True,
    )


def replay(inputs, panel, output, engines):
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    records, events = [], []
    for row in csv.DictReader((panel / "runs.csv").open()):
        if row["arm"] not in engines or row["cap"] != "3000":
            continue
        key, pdb, arm = (row[n] for n in ("key", "pdb", "arm"))
        directory = inputs / pdb
        case = (
            pdb,
            "",
            directory / "receptor.pdbqt",
            directory / "ligand.pdbqt",
            panel / (pdb + ".box"),
        )
        if list(raw.glob(key + ".*.pdbqt")):
            raise ValueError("Refuse to append duplicate observer output: " + key)
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "result.pdbqt"
            log = subprocess.check_output(
                command(case, row["seed"], row["cap"], engines[arm], out),
                env={**os.environ, "RING_OBSERVER": str((raw / key).resolve())},
                stderr=subprocess.STDOUT,
                text=True,
            )
            if sha(out) != row["sha256"] or any(
                int(row[k]) != v for k, v in counts(log).items()
            ):
                raise ValueError("Observer differs: " + key)
        (raw / (key + ".log")).write_text(log)
        pending = None
        for line in log.splitlines():
            if line.startswith("Ring proposal:"):
                m = re.fullmatch(
                    r"Ring proposal: axis=(\d+) draw=(\S+) applied=(\d+)", line
                )
                if not m or pending:
                    raise ValueError(line)
                pending = m.groups()
            if line.startswith("Ring step:") and pending:
                axis, draw, applied = pending
                events.append(
                    dict(
                        key=key,
                        step=int(line.split()[-1]) + 1,
                        axis=axis,
                        draw_rad=draw,
                        applied=applied,
                    )
                )
                pending = None
        if pending:
            raise ValueError(key)
        for stage in ("before", "proposal", "hunt", "accepted", "second"):
            sk = key + "." + stage
            stages = (raw / (sk + ".pdbqt")).read_text()
            ids = [
                int(s.split()[1]) for s in stages.splitlines() if s.startswith("MODEL ")
            ]
            records.append(
                dict(
                    key=sk,
                    pdb=pdb,
                    arm=arm,
                    seed=row["seed"],
                    cap=row["cap"],
                    stage=stage,
                    modes=len(ids),
                )
            )
    write(output / "runs.csv", records)
    write(output / "events.csv", events)
    score(inputs, output)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    q = sub.add_parser("build")
    for k in ["source", "objects", "output"]:
        q.add_argument(k, type=Path)
    q.add_argument("--boost-prefix", type=Path)
    q = sub.add_parser("replay")
    for k in ["inputs", "panel", "output", "null", "concerted"]:
        q.add_argument(k, type=Path)
    a = p.parse_args()
    if a.action == "build":
        observe(a.source, a.objects, a.output, a.boost_prefix)
    else:
        replay(a.inputs, a.panel, a.output, {"null": a.null, "concerted": a.concerted})
