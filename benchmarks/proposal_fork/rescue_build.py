"""Export and observe the production last-trial solver and its same-state counterfactual."""

import argparse
from pathlib import Path
import subprocess
import sys

from build import SOLVERS, edit, prepare

REV = "1c672a8e2ce12c7a16e20172a3315ea5a19c23d8"


def instrument(out):
    SOLVERS["last_trial"] = REV
    prepare("last_trial", out)
    lib = out / "src/lib"
    header = lib / "bfgs.h"
    qn = lib / "quasi_newton.cpp"
    trace = lib / "bfgs_trace.h"
    mc = lib / "monte_carlo.cpp"
    edit(mc, "quasi_newton_par, path);", "quasi_newton_par, hunt_grid_spacing, path);")
    edit(
        trace,
        "quasi_newton& local, const char* path",
        "quasi_newton& local, fl spacing, const char* path",
    )
    edit(
        trace,
        "local(m, p, ig, pose, cap, cost);",
        "local(m, p, ig, pose, cap, cost, nullptr, 0, spacing);",
    )
    edit(
        header,
        "extern thread_local bool enabled;",
        """extern thread_local bool enabled;
extern thread_local bool counterfactual;
void rescue(fl, fl, fl, fl, fl, fl, bool, bool);
""",
    )
    edit(
        qn,
        "thread_local bool enabled = false;",
        """thread_local bool enabled = false;
thread_local bool counterfactual = std::getenv("FD_FORK_COUNTERFACTUAL") != nullptr;
void rescue(fl pg, fl f0, fl old_a, fl new_a, fl old_e, fl new_e, bool old_pass, bool new_pass) {
    static std::ofstream out(std::getenv("FD_FORK_RESCUES"));
    VINA_CHECK(out.good());
    out << std::setprecision(17) << input << ',' << pg << ',' << f0 << ',' << old_a << ',' << new_a << ','
        << old_e << ',' << new_e << ',' << old_pass << ',' << new_pass << '\\n';
}
""",
    )
    edit(
        header,
        """\t\tif (trial + 1 == max_trials && last_trial)
\t\t\talpha = last_trial(x, p, alpha);""",
        """        fl alternate_alpha = alpha, alternate_energy = 0;
        bool shadow = trial + 1 == max_trials && last_trial && fork_probe::counterfactual;
        if (trial + 1 == max_trials && last_trial) {
            alternate_alpha = last_trial(x, p, alpha);
            if (shadow) {
                VINA_CHECK(!StrictArmijo);
                conf alternate = x;
                alternate.set_incremented(x, p, alternate_alpha);
                flv alternate_gradient(g_new.size());
                alternate_energy = f(alternate, alternate_gradient);
            } else alpha = alternate_alpha;
        }""",
    )
    edit(
        header,
        "        if (fork_probe::enabled) {\n            const bool pass",
        """        if (shadow)
            fork_probe::rescue(pg, f0, alpha, alternate_alpha, f1, alternate_energy,
                f1-f0 < c0*alpha*pg, alternate_energy-f0 < c0*alternate_alpha*pg);
        if (fork_probe::enabled) {
            const bool pass""",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument("--boost", type=Path, required=True)
    p.add_argument("--jobs", default="3")
    a = p.parse_args()
    out = a.out.resolve()
    instrument(out)
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
