"""Apply a read-only observer to an exported historical tree; build with its tools/build.py."""
import argparse
from pathlib import Path
import shutil


def edit(path, before, after):
    text = path.read_text()
    if text.count(before) != 1:
        raise ValueError(f"{path}: expected exactly one {before!r}")
    path.write_text(text.replace(before, after))


def instrument(source):
    lib = source / "src/lib"
    mc, bfgs, vina = (lib / p for p in ("monte_carlo.cpp", "bfgs.h", "vina.cpp"))
    shutil.copyfile(Path(__file__).with_name("trace.h"), lib / "bfgs_trace.h")
    edit(mc, '#include "trajectory_archive.h"', '''#include "trajectory_archive.h"
#include "bfgs_trace.h"
namespace bfgs_trace {
thread_local int failures = 0, non_descent = 0;
thread_local fl first_angle = -1;
}''')
    edit(bfgs, '#include <cstring>', '''#include <cstring>
namespace bfgs_trace {
extern thread_local int failures, non_descent;
extern thread_local fl first_angle;
}''')
    edit(bfgs, 'return 0;\n\tconst fl c0', '{ ++bfgs_trace::non_descent; return 0; }\n\tconst fl c0')
    edit(bfgs, '\t\talpha *= multiplier;', '''        if (trial + 1 == max_trials) {
            ++bfgs_trace::failures;
            if (bfgs_trace::first_angle < 0 && p.size() >= 6) {
                fl angle = std::sqrt(p[3]*p[3] + p[4]*p[4] + p[5]*p[5]);
                for (sz i = 6; i < p.size(); ++i) angle = std::max(angle, std::abs(p[i]));
                bfgs_trace::first_angle = angle;
            }
        }
        alpha *= multiplier;''')
    edit(mc, 'VINA_U_FOR(step, global_steps)', '''if (const char* path = std::getenv("FD_BFGS_SHADOW"))
        bfgs_trace::shadow(m, p, ig, hunt_cap, local_steps, path);
    bfgs_trace::stream() << "T\\n";
    VINA_U_FOR(step, global_steps)''')
    edit(mc, 'candidate.c = tmp.c;', '''bfgs_trace::pose('S', step, tmp);
        candidate.c = tmp.c;''')
    edit(mc, 'const int before_hunt = evalcount;', '''bfgs_trace::pose('P', step, candidate);
        bfgs_trace::reset();
        const int before_hunt = evalcount;''')
    edit(mc, 'if (counters) {\n\t\t\t++counters->hunt_calls;', '''bfgs_trace::result('H', step, candidate, m, evalcount - before_hunt);
        if (counters) {
            ++counters->hunt_calls;''')
    edit(mc, 'std::swap(tmp.c, candidate.c);', '''bfgs_trace::stream() << "A\\t" << step << '\\n';
            std::swap(tmp.c, candidate.c);''')
    edit(mc, 'quasi_newton_par(m, p, ig, saved, authentic_v, evalcount);', '''bfgs_trace::reset();
                quasi_newton_par(m, p, ig, saved, authentic_v, evalcount);
                bfgs_trace::result('D', step, saved, m, evalcount - before_second);''')
    edit(mc, 'reference_from_save = !separate_hunt_chain;', '''for (const auto& retained : out) bfgs_trace::pose('M', step, retained);
                reference_from_save = !separate_hunt_chain;''')
    edit(vina, '#include "vina.h"', '#include "vina.h"\n#include "bfgs_trace.h"')
    edit(vina, 'VINA_FOR(p, 5){', '''bfgs_trace::pose('I', i, poses[i]);
                    VINA_FOR(p, 5){''')
    edit(vina, 'poses[i].coords = late_selection ?', '''bfgs_trace::result('O', i, poses[i], m_model, evalcount);
                    poses[i].coords = late_selection ?''')
    edit(vina, '++reported;', '''bfgs_trace::result('F', reported, poses[i], m_model, 0);
                ++reported;''')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    instrument(parser.parse_args().source)
