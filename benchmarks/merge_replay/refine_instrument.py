"""Extend the clean candidate merge observer with the registered single-slot substitution."""
from pathlib import Path
import sys

lib = Path(sys.argv[1]) / "src/lib"


def edit(path, before, after):
    old = path.read_text()
    assert old.count(before) == 1, (path, before)
    path.write_text(old.replace(before, after))


edit(lib / "parallel_mc.cpp", 'members("G", sorted);', '''members("G", sorted);
    VINA_CHECK(out.size() == 9 && all.size() == 27);
    VINA_CHECK(identifiers.at(bfgs_trace::state(out[0].c)) == 18);
    if (std::getenv("FD_SUBSTITUTE_20"))
        out[0] = all[20];
    members("U", out);''')
path = lib / "vina.cpp"
edit(path, '#include "vina.h"', '#include "vina.h"\n#include "bfgs_trace.h"')
edit(path, 'VINA_FOR(p, 5){', '''const int before_refine = evalcount;
                    bfgs_trace::pose('I', i, poses[i]);
                    VINA_FOR(p, 5){''')
edit(path, 'poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra;', '''poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra;
                    bfgs_trace::stream() << "O\\t" << i << '\\t' << bfgs_trace::state(poses[i].c)
                        << '\\t' << evalcount-before_refine << '\\t' << all_grids << '\\t' << inter_pairs
                        << '\\t' << intra_pairs << '\\t' << lig_intra << '\\n';''')
edit(path, 'std::vector<double> energies = score_current_conf(intramolecular_energy);', '''std::vector<double> energies = score_current_conf(intramolecular_energy);
            const char* fixed = std::getenv("FD_COMMON_REFERENCE");
            const double anchor = fixed ? std::stod(fixed) : intramolecular_energy;
            const std::vector<double> fixed_score = score_current_conf(anchor);
            auto& log = bfgs_trace::stream();
            log << "R\\t" << i << '\\t' << bfgs_trace::state(poses[i].c) << '\\t' << intramolecular_energy
                << '\\t' << anchor << '\\t' << fixed_score[0];
            for (const double term : energies) log << '\\t' << term;
            log << '\\n';''')
edit(path, '++reported;', '''bfgs_trace::pose('F', reported, poses[i]);
                ++reported;''')
