"""Observe MC's configured optimizer on frozen proposals using the production object."""
import argparse
from pathlib import Path
from instrument import edit


def instrument(source):
    lib = source / "src/lib"
    path = lib / "bfgs_trace.h"
    path.write_text(Path(__file__).with_name("trace.h").read_text())
    edit(path, "unsigned steps, const char* path", "quasi_newton& local, const char* path")
    edit(path, "    quasi_newton local(size, steps);", "")
    edit(path, "extern thread_local fl first_angle;", """extern thread_local fl first_angle;
extern thread_local fl f_start, pg_start, fail_pg, fail_last_angle;""")
    edit(path, "first_angle = -1; }", "first_angle = -1; f_start = pg_start = fail_pg = fail_last_angle = NAN; }")
    edit(path, "    out << '\\n';", """    out << '\\t' << f_start << '\\t' << pg_start << '\\t' << fail_pg << '\\t' << fail_last_angle << '\\n';""")
    mc = lib / "monte_carlo.cpp"
    edit(mc, '#include "trajectory_archive.h"', '''#include "trajectory_archive.h"
#include "bfgs_trace.h"
namespace bfgs_trace {
thread_local int failures = 0, non_descent = 0;
thread_local fl first_angle = -1, f_start = NAN, pg_start = NAN, fail_pg = NAN, fail_last_angle = NAN;
}''')
    edit(mc, "VINA_U_FOR(step, global_steps)", '''if (const char* path = std::getenv("FD_BFGS_SHADOW"))
        bfgs_trace::shadow(m, p, ig, hunt_cap, quasi_newton_par, path);
    VINA_U_FOR(step, global_steps)''')
    bfgs = lib / "bfgs.h"
    edit(bfgs, '#include <cstring>', '''#include <cstring>
namespace bfgs_trace {
extern thread_local int failures, non_descent;
extern thread_local fl first_angle, f_start, pg_start, fail_pg, fail_last_angle;
}''')
    edit(bfgs, 'if (!std::isfinite(pg) || pg >= 0)\n\t\treturn 0;', '''if (std::isnan(bfgs_trace::f_start)) {
        bfgs_trace::f_start = f0;
        bfgs_trace::pg_start = pg;
    }
    if (!std::isfinite(pg) || pg >= 0) {
        ++bfgs_trace::non_descent;
        bfgs_trace::fail_pg = pg;
        return 0;
    }''')
    edit(bfgs, 'alpha *= multiplier;', '''if (trial + 1 == max_trials) {
            ++bfgs_trace::failures;
            fl angle = std::hypot(std::hypot(p[3], p[4]), p[5]);
            for (sz i = 6; i < p.size(); ++i) angle = std::max(angle, std::abs(p[i]));
            bfgs_trace::first_angle = angle;
            bfgs_trace::fail_last_angle = alpha * angle;
            bfgs_trace::fail_pg = pg;
        }
        alpha *= multiplier;''')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    instrument(parser.parse_args().source)
