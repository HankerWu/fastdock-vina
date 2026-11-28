// Benchmark observer for the historical single-ligand, single-worker CLI.
#pragma once
#include "quasi_newton.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace bfgs_trace {
extern thread_local int failures, non_descent;
extern thread_local fl first_angle;

inline std::ostream& stream() {
    static std::ofstream sink(std::getenv("FD_BFGS_TRACE"));
    VINA_CHECK(sink.good());
    return sink << std::setprecision(17);
}
inline void reset() { failures = non_descent = 0; first_angle = -1; }
inline std::string state(const conf& c) {
    VINA_CHECK(c.ligands.size() == 1 && c.flex.empty());
    const auto& ligand = c.ligands.front();
    const auto& r = ligand.rigid;
    const auto& q = r.orientation;
    std::ostringstream s;
    s << std::setprecision(17) << r.position[0] << ' ' << r.position[1] << ' ' << r.position[2]
      << ' ' << q.R_component_1() << ' ' << q.R_component_2()
      << ' ' << q.R_component_3() << ' ' << q.R_component_4();
    for (fl t : ligand.torsions) s << ' ' << t;
    return s.str();
}
inline void decode(std::istream& s, conf& c) {
    auto& ligand = c.ligands.front();
    auto& r = ligand.rigid;
    fl a, b, d, e;
    s >> r.position[0] >> r.position[1] >> r.position[2] >> a >> b >> d >> e;
    r.orientation = qt(a, b, d, e);
    for (fl& t : ligand.torsions) s >> t;
    VINA_CHECK(!s.fail());
}
inline void pose(char event, sz step, const output_type& p) {
    stream() << event << '\t' << step << '\t' << p.e << '\t' << state(p.c) << '\n';
}
inline void result(char event, sz step, const output_type& p, const model& m, int evals) {
    auto& out = stream();
    out << event << '\t' << step << '\t' << p.e << '\t' << state(p.c)
        << '\t' << evals << '\t' << failures << '\t' << non_descent << '\t' << first_angle << '\t';
    // Chemical atoms in ORIGINAL ligand PDBQT record order. Millangstroms for
    // Cartesian comparison, exact binary64 conformation round-trips above.
    for (const auto& atom : m.chemical_heavy_atom_elements()) {
        const auto& xyz = m.movable_coords(atom.first);
        for (sz i=0; i<3; ++i) out << ' ' << std::llround(1000 * xyz[i]);
    }
    out << '\n';
}
inline void shadow(model& m, const precalculate_byatom& p, const igrid& ig,
                   const vec& cap, unsigned steps, const char* path) {
    std::ifstream input(path);
    VINA_CHECK(input.good());
    const auto size = m.get_size();
    quasi_newton local(size, steps);
    std::string line;
    sz index = 0;
    while (std::getline(input, line)) {
        output_type pose(size, 0);
        std::istringstream fields(line);
        decode(fields, pose.c);
        int cost = 0;
        reset();
        local(m, p, ig, pose, cap, cost);
        result('H', index++, pose, m, cost);
    }
    VINA_CHECK(input.eof() && index);
    stream().flush();
    std::exit(EXIT_SUCCESS);
}
}
