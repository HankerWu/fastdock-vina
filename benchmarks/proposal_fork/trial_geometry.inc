// Compiled only in the source-pinned benchmark observer.
#include "cache.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>

namespace fork_probe {
thread_local model* moving = nullptr;
thread_local const cache* maps = nullptr;
thread_local long input = -1;
thread_local int iteration = -1;
thread_local bool enabled = false;

void select(sz n) {
    static const std::set<sz> indices = [] {
        std::set<sz> chosen;
        std::ifstream s(std::getenv("FD_FORK_INDICES"));
        sz i;
        while (s >> i) chosen.insert(i);
        VINA_CHECK(s.eof());
        return chosen;
    }();
    input = n;
    iteration = -1;
    enabled = indices.count(n);
}

fl rms_move(const vecv& a, const vecv& b) {
    fl sum = 0;
    VINA_FOR_IN(i, a) sum += vec_distance_sqr(a[i], b[i]);
    return std::sqrt(sum / a.size());
}

void write(const conf& x, const conf& trial, unsigned attempt, fl alpha, fl pg,
           fl fd4, fl fd5, fl f0, fl f1, bool accepted) {
    static std::ofstream out(std::getenv("FD_FORK_TRIALS"));
    VINA_CHECK(out.good() && moving && maps);
    if (attempt == 0) ++iteration;
    const auto xyz = moving->get_chemical_heavy_atom_movable_coords();
    moving->set(x);
    const auto base = moving->get_chemical_heavy_atom_movable_coords();
    conf component = x;
    auto& ligand = component.ligands.front();
    ligand.rigid.orientation = trial.ligands.front().rigid.orientation;
    moving->set(component);
    const fl rotation = rms_move(base, moving->get_chemical_heavy_atom_movable_coords());
    component = x;
    component.ligands.front().torsions = trial.ligands.front().torsions;
    moving->set(component);
    const fl torsion = rms_move(base, moving->get_chemical_heavy_atom_movable_coords());
    moving->set(trial);
    const fl translation = std::sqrt(vec_distance_sqr(x.ligands.front().rigid.position,
                                         trial.ligands.front().rigid.position));
    const grid_dims grid = maps->gd();
    int crossed = 0, out0 = 0, out1 = 0;
    VINA_FOR_IN(i, xyz) {
        bool face = false, outside_base = false, outside_trial = false;
        VINA_FOR(axis, 3) {
            const auto& g = grid[axis];
            const auto bin = [&](fl coordinate) {
                const fl offset = (coordinate - g.begin) * g.n_voxels / g.span();
                return std::max(fl(0), std::min(std::floor(offset), fl(g.n_voxels-1)));
            };
            face |= bin(base[i][axis]) != bin(xyz[i][axis]);
            outside_base |= base[i][axis] < g.begin || base[i][axis] > g.end;
            outside_trial |= xyz[i][axis] < g.begin || xyz[i][axis] > g.end;
        }
        crossed += face;
        out0 += outside_base;
        out1 += outside_trial;
    }
    out << std::setprecision(17) << input << ',' << iteration << ',' << attempt << ','
        << alpha << ',' << pg << ',' << fd4 << ',' << fd5 << ',' << f0 << ',' << f1 << ','
        << alpha * pg << ',' << f1-f0 << ',' << translation << ',' << rotation << ','
        << torsion << ',' << rms_move(base, xyz) << ',' << crossed << ',' << out0 << ','
        << out1 << ',' << accepted << '\n';
}
}
