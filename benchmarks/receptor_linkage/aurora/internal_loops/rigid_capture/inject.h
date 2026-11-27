// Single-CPU diagnostic. Replace one real proposal's rigid transform, keeping
// every torsion and the pre-intervention RNG stream; print the common start.
#pragma once
#include <array>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>

struct rigid_injection {
    unsigned trajectory = 0, step = 0;
    std::array<fl, 7> value{};
    std::ofstream log;
    rigid_injection() {
        const char* path = std::getenv("FD_RIGID_INTERVENTION");
        if (!path) return;
        std::ifstream f(path);
        f >> trajectory >> step;
        for (fl& v : value) f >> v;
        VINA_CHECK(f.good() && trajectory != 0 && step != 0);
        const char* dest = std::getenv("FD_RIGID_LOG");
        VINA_CHECK(dest);
        log.open(dest);
        VINA_CHECK(log.is_open());
        log << std::setprecision(17);
    }
    bool apply(unsigned tr, unsigned n, conf& c) {
        if (tr != trajectory || n != step) return false;
        VINA_CHECK(c.ligands.size() == 1 && c.flex.empty());
        auto& ligand = c.ligands.front();
        auto& rigid = ligand.rigid;
        log << "source";
        for (sz i = 0; i < 3; ++i) log << ' ' << rigid.position[i];
        const qt& q = rigid.orientation;
        log << ' ' << q.R_component_1() << ' ' << q.R_component_2()
            << ' ' << q.R_component_3() << ' ' << q.R_component_4();
        for (fl v : ligand.torsions) log << ' ' << v;
        log << '\n';
        for (sz i = 0; i < 3; ++i) rigid.position[i] = value[i];
        rigid.orientation = qt(value[3], value[4], value[5], value[6]);
        VINA_CHECK(quaternion_is_normalized(rigid.orientation));
        return true;
    }
};

// Snapshot coordinates in PDBQT serial order. Read-only at active states;
// callers sampling stored generalized conformations use a separate model copy.
inline void rigid_stage_xyz(const char* stage, unsigned tr, unsigned step, const model& m) {
    static const char* path = std::getenv("FD_RIGID_COORDS");
    if (!path) return;
    static std::ofstream out(path);
    static szv order;
    if (order.empty()) {
        VINA_CHECK(out.good() && m.num_ligands() == 1 && m.num_flex() == 0);
        std::set<sz> heavy;
        for (const auto& entry : m.chemical_heavy_atom_elements()) heavy.insert(entry.first);
        out << '#' << std::setprecision(17);
        for (const auto& line : m.get_ligand(0).cont) {
            if (line.second && heavy.count(*line.second)) {
                order.push_back(*line.second);
                out << ' ' << std::stoul(line.first.substr(6, 5));
            }
        }
        out << '\n';
        VINA_CHECK(order.size() == heavy.size());
    }
    out << stage << ' ' << tr << ' ' << step;
    for (sz index : order)
        for (sz d = 0; d != 3; ++d) out << ' ' << m.movable_coords(index)[d];
    out << '\n';
}
