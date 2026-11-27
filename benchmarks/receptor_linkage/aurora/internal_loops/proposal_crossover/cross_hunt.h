// Diagnostic for a single rigid receptor, a single ligand, and --cpu 1.
#pragma once
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>

namespace cross_hunt {
struct Recorder {
    unsigned trajectory = 0;
    std::ofstream proposals, outputs;
    Recorder() {
        if (const char* path = std::getenv("FD_DUMP_PROPOSALS")) {
            proposals.open(path);
            VINA_CHECK(proposals.is_open());
            proposals << std::setprecision(17);
        }
        if (const char* path = std::getenv("FD_QUENCH_OUTPUT")) {
            outputs.open(path);
            VINA_CHECK(outputs.is_open());
            outputs << std::setprecision(17);
        }
    }
    void record(unsigned step, const conf& c) {
        if (!proposals.is_open()) return;
        VINA_CHECK(c.ligands.size() == 1 && c.flex.empty());
        const auto& ligand = c.ligands.front();
        const auto& q = ligand.rigid.orientation;
        proposals << trajectory << ' ' << step;
        for (sz j = 0; j < 3; ++j) proposals << ' ' << ligand.rigid.position[j];
        proposals << ' ' << q.R_component_1() << ' ' << q.R_component_2()
                  << ' ' << q.R_component_3() << ' ' << q.R_component_4();
        for (fl v : ligand.torsions) proposals << ' ' << v;
        proposals << '\n';
    }
    void quench(model m, const precalculate_byatom& p, const igrid& ig,
                const vec& cap, unsigned local_steps) {
        if (!outputs.is_open()) return;
        const char* path = std::getenv("FD_QUENCH_INPUT");
        VINA_CHECK(path && m.num_ligands() == 1 && m.num_flex() == 0);
        std::ifstream in(path);
        VINA_CHECK(in.good());
        std::set<sz> heavy;
        for (const auto& atom : m.chemical_heavy_atom_elements()) heavy.insert(atom.first);
        szv order;
        outputs << "#";
        for (const auto& atom : m.get_ligand(0).cont) {
            if (atom.second && heavy.count(*atom.second)) {
                order.push_back(*atom.second);
                outputs << ' ' << std::stoul(atom.first.substr(6, 5));
            }
        }
        outputs << '\n';
        VINA_CHECK(order.size() == heavy.size());
        const conf_size s = m.get_size();
        quasi_newton hunt(s, local_steps);
        unsigned tr, step;
        while (in >> tr >> step) {
            output_type candidate(s, 0);
            auto& ligand = candidate.c.ligands.front();
            for (sz j = 0; j < 3; ++j) in >> ligand.rigid.position[j];
            fl a, b, c, d;
            in >> a >> b >> c >> d;
            ligand.rigid.orientation = qt(a, b, c, d);
            for (auto& angle : ligand.torsions) in >> angle;
            VINA_CHECK(in.good() && quaternion_is_normalized(ligand.rigid.orientation));
            int evals = 0;
            hunt(m, p, ig, candidate, cap, evals);
            outputs << tr << ' ' << step << ' ' << evals << ' ' << candidate.e;
            for (sz i : order)
                for (sz j = 0; j < 3; ++j) outputs << ' ' << m.movable_coords(i)[j];
            outputs << '\n';
        }
        VINA_CHECK(in.eof() && outputs.good());
    }
};
}
