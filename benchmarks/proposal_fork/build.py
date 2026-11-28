"""Export the two production solver identities and add a benchmark-only shadow observer."""

import argparse
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bfgs_source_replay"))
from instrument_candidate import edit, instrument

SOLVERS = {
    "stable": "0449b1d0bc3adec768dd9a9f243f3ce5a29cef58",
    "candidate": "0bcae6ef09857f97656291bb4821d20b3b619bf5",
}


def prepare(arm, out):
    out.mkdir(parents=True, exist_ok=False)
    with subprocess.Popen(
        [
            "git",
            "-C",
            HERE.parents[1],
            "archive",
            SOLVERS[arm],
            "src",
            "tools/build.py",
        ],
        stdout=subprocess.PIPE,
    ) as p:
        with tarfile.open(fileobj=p.stdout, mode="r|") as archive:
            archive.extractall(out, filter="data")
        if p.wait():
            raise RuntimeError("source export failed")
    instrument(out)
    lib = out / "src/lib"
    trace = lib / "bfgs_trace.h"
    edit(
        trace,
        "namespace bfgs_trace {",
        "namespace fork_probe { void select(sz); }\nnamespace bfgs_trace {",
    )
    edit(
        trace,
        "        reset();",
        "        fork_probe::select(index);\n        reset();",
    )
    if arm == "candidate":
        edit(
            trace,
            "local(m, p, ig, pose, cap, cost);",
            "local(m, p, ig, pose, cap, cost, nullptr, 0, true);",
        )
    qn = lib / "quasi_newton.cpp"
    edit(
        qn,
        '#include "quasi_newton.h"',
        '#include "quasi_newton.h"\n' + (HERE / "trial_geometry.inc").read_text(),
    )
    edit(
        qn,
        "    const auto rigid =",
        """    fork_probe::moving = &m;
    fork_probe::maps = dynamic_cast<const cache*>(&ig);
    const auto rigid =""",
    )
    header = lib / "bfgs.h"
    edit(
        header,
        "#include <cstring>",
        """#include <cstring>
namespace fork_probe {
extern thread_local bool enabled;
void write(const conf&, const conf&, unsigned, fl, fl, fl, fl, fl, fl, bool);
}""",
    )
    edit(
        header,
        "\tVINA_U_FOR(trial, max_trials)",
        """    fl fd[2] = {};
    if (fork_probe::enabled) {
        fl max_component = 0;
        for (fl component : p) max_component = std::max(max_component, std::abs(component));
        VINA_CHECK(max_component > 0);
        for (unsigned j=0; j<2; ++j) {
            const fl h = (j == 0 ? 1e-4 : 1e-5) / max_component;
            flv scratch(g_new.size());
            conf probe = x;
            probe.set_incremented(x, p, h);
            const fl plus = f(probe, scratch);
            probe.set_incremented(x, p, -h);
            const fl minus = f(probe, scratch);
            fd[j] = (plus - minus) / (2*h);
        }
    }
    VINA_U_FOR(trial, max_trials)""",
    )
    edit(
        header,
        "\t\t++evalcount;\n\t\tif (StrictArmijo)",
        """        ++evalcount;
        if (fork_probe::enabled) {
            const bool pass = StrictArmijo
                ? std::isfinite(f1) && f1 < f0 && f1-f0 <= c0*alpha*pg &&
                    std::all_of(g_new.begin(), g_new.end(), [](fl g){return std::isfinite(g);})
                : f1-f0 < c0*alpha*pg;
            fork_probe::write(x, x_new, trial, alpha, pg, fd[0], fd[1], f0, f1, pass);
        }
        if (StrictArmijo)""",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("arm", choices=SOLVERS)
    p.add_argument("out", type=Path)
    p.add_argument("--boost", type=Path, required=True)
    p.add_argument("--jobs", default="3")
    a = p.parse_args()
    out = a.out.resolve()
    prepare(a.arm, out)
    subprocess.run(
        [
            sys.executable,
            out / "tools/build.py",
            "--boost-prefix",
            a.boost.resolve(),
            "--output-dir",
            out / "build",
            "--jobs",
            a.jobs,
        ],
        check=True,
    )
