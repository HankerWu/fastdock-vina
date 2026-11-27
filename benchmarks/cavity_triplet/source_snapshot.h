// Research observer. A single-CPU process records one native post-mutation
// state on a model copy. No evaluations and no RNG draws occur here.
#pragma once
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>

struct source_snapshot {
    unsigned trajectory = 0, step = 0;
    source_snapshot() {
        if (const char* value = std::getenv("FD_SOURCE_TRAJECTORY")) {
            trajectory = std::stoul(value);
            step = std::stoul(std::getenv("FD_SOURCE_STEP"));
        }
    }
    void capture(unsigned tr, unsigned n, const model& original, const conf& c) const {
        if (tr != trajectory || n != step) return;
        model state(original);
        state.set(c);
        VINA_CHECK(c.ligands.size() == 1 && c.flex.empty());
        std::ofstream g(std::getenv("FD_SOURCE_Q")), x(std::getenv("FD_SOURCE_XYZ"));
        VINA_CHECK(g.good() && x.good());
        const auto& ligand = c.ligands.front();
        const auto& q = ligand.rigid.orientation;
        g << std::setprecision(17) << tr << ' ' << n;
        for (sz j = 0; j < 3; ++j) g << ' ' << ligand.rigid.position[j];
        g << ' ' << q.R_component_1() << ' ' << q.R_component_2()
          << ' ' << q.R_component_3() << ' ' << q.R_component_4();
        for (fl v : ligand.torsions) g << ' ' << v;
        g << '\n';
        std::set<sz> heavy;
        for (const auto& entry : state.chemical_heavy_atom_elements()) heavy.insert(entry.first);
        szv ordered;
        x << std::setprecision(17) << '#';
        for (const auto& record : state.get_ligand(0).cont)
            if (record.second && heavy.count(*record.second)) {
                ordered.push_back(*record.second);
                x << ' ' << std::stoul(record.first.substr(6, 5));
            }
        VINA_CHECK(ordered.size() == heavy.size());
        x << "\nsource " << tr << ' ' << n;
        for (sz index : ordered)
            for (sz j = 0; j < 3; ++j) x << ' ' << state.movable_coords(index)[j];
        x << '\n';
    }
};
