"""Instrument an exported 0bb827e tree. The two stage flags only select its original scalar policy."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bfgs_source_replay"))
from instrument_candidate import instrument as candidate_observer, edit


def instrument(source):
    candidate_observer(source)
    lib = source / "src/lib"
    trace, bfgs, mc, vina = [lib / f for f in ("bfgs_trace.h", "bfgs.h", "monte_carlo.cpp", "vina.cpp")]
    decl = "pg_start, fail_pg, fail_last_angle;"
    for f in (bfgs, trace):
        edit(f, decl, decl + """
extern thread_local int capped;
extern thread_local fl shift_before, shift_after;
extern thread_local long long elapsed_ns;""")
    edit(trace, '#include <cstdlib>', '#include <cstdlib>\n#include <time.h>')
    edit(trace, 'inline void reset()', '''inline long long ns() {
    timespec t; VINA_CHECK(clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t) == 0);
    return 1000000000LL*t.tv_sec + t.tv_nsec;
}
inline void reset()''')
    edit(trace, 'fail_last_angle = NAN; }','fail_last_angle = NAN; capped = 0; shift_before = shift_after = 0; elapsed_ns = 0; }')
    edit(trace, " << fail_last_angle << '\\n';", " << fail_last_angle << '\\t' << capped << '\\t' << shift_before << '\\t' << shift_after << '\\t' << elapsed_ns << '\\n';")
    edit(trace, 'const vec& cap, quasi_newton& local, const char* path', 'const vec& hunt, const vec& save, quasi_newton& local, const char* path')
    edit(trace,'decode(fields, pose.c);','''char phase; fields >> phase;
        VINA_CHECK(phase == 'H' || phase == 'D');
        decode(fields, pose.c);''')
    edit(trace,'local(m, p, ig, pose, cap, cost);', '''const auto begin = ns();
        local(m, p, ig, pose, phase == 'H' ? hunt : save, cost);
        elapsed_ns = ns() - begin;''')
    edit(trace, "result('H', index++", "result(phase, index++")
    edit(bfgs, 'fl alpha = limit_angles ? principal_angle_step(x, p) : 1;', '''fl alpha = limit_angles ? principal_angle_step(x, p) : 1;
    if (limit_angles && alpha < 1) {
        const fl norm = std::hypot(std::hypot(p[0], p[1]), p[2]);
        if (norm > 0) {
            ++bfgs_trace::capped;
            bfgs_trace::shift_before += norm;
            bfgs_trace::shift_after += alpha * norm;
        }
    }''')
    edit(mc, 'thread_local int failures = 0, non_descent = 0;', '''thread_local int failures = 0, non_descent = 0, capped = 0;
thread_local fl shift_before = 0, shift_after = 0;
thread_local long long elapsed_ns = 0;''')
    edit(mc,'quasi_newton_par.limit_angles = true;', '''const std::string mode = std::getenv("FD_STAGE_MODE");
    VINA_CHECK(mode == "00" || mode == "10" || mode == "01" || mode == "11");
    const bool hunt_limit = mode[0] == '1', save_limit = mode[1] == '1';
    quasi_newton_par.limit_angles = hunt_limit;''')
    edit(mc, 'ig, hunt_cap, quasi_newton_par, path);', 'ig, hunt_cap, authentic_v, quasi_newton_par, path);')
    edit(mc, 'VINA_U_FOR(step, global_steps)', '''bfgs_trace::stream() << "T\\n";
    VINA_U_FOR(step, global_steps)''')
    edit(mc,'candidate.c = tmp.c;', '''bfgs_trace::pose('S', step, tmp);
        candidate.c = tmp.c;''')
    edit(mc,'const int before_hunt = evalcount;', '''bfgs_trace::pose('P', step, candidate);
        quasi_newton_par.limit_angles = hunt_limit;
        bfgs_trace::reset();
        const int before_hunt = evalcount;
        const auto hunt_clock = bfgs_trace::ns();''')
    edit(mc,'quasi_newton_par(m, p, ig, candidate, hunt_cap, evalcount);', '''quasi_newton_par(m, p, ig, candidate, hunt_cap, evalcount);
        bfgs_trace::elapsed_ns = bfgs_trace::ns() - hunt_clock;
        bfgs_trace::result('H', step, candidate, m, evalcount - before_hunt);''')
    edit(mc, 'std::swap(tmp.c, candidate.c);', '''bfgs_trace::stream() << "A\\t" << step << '\\n';
            std::swap(tmp.c, candidate.c);''')
    edit(mc, 'quasi_newton_par(m, p, ig, saved, authentic_v, evalcount);', '''quasi_newton_par.limit_angles = save_limit;
                bfgs_trace::reset();
                const auto save_clock = bfgs_trace::ns();
                quasi_newton_par(m, p, ig, saved, authentic_v, evalcount);
                bfgs_trace::elapsed_ns = bfgs_trace::ns() - save_clock;
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
