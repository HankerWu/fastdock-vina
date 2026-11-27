// Check the production pair list against independent shortest paths in the
// parsed ligand graph. --forces emits actual per-atom intraligand gradients
// keyed by PDBQT serial, including the component from forbidden pairs.
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "parse_pdbqt.h"
#include "precalculate.h"

namespace {
void check(const char* path, bool forces, bool require_valid) {
    model m = parse_ligand_pdbqt_from_file(path, atom_type::XS);
    const sz n = m.num_atoms();
    std::vector<std::vector<sz>> distance(n, std::vector<sz>(n, n + 1));
    for (sz i = 0; i < n; ++i) {
        distance[i][i] = 0;
        for (const bond& b : m.get_atom(i).bonds)
            if (!b.connected_atom_index.in_grid)
                distance[i][b.connected_atom_index.i] = 1;
    }
    // Floyd-Warshall neither uses the engine's traversal order nor truncates paths.
    for (sz k = 0; k < n; ++k)
        for (sz i = 0; i < n; ++i)
            for (sz j = 0; j < n; ++j)
                distance[i][j] = std::min(distance[i][j], distance[i][k] + distance[k][j]);
    const interacting_pairs pairs = m.get_ligand(0).pairs;
    interacting_pairs spurious;
    for (const auto& p : pairs)
        if (distance[p.a][p.b] <= 3) spurious.push_back(p);

    std::map<sz,unsigned> serial;
    for (const auto& kv : ligand_heavy_indices_by_serial(m))
        serial.emplace(kv.second, kv.first);
    if (forces) {
        const ScoringFunction sf(SF_VINA, {-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50,0.05846});
        const precalculate_byatom tables(sf, m);
        vecv xyz, all(n, vec(0,0,0)), erroneous(n, vec(0,0,0));
        for (sz i = 0; i < n; ++i) xyz.push_back(m.get_coords(i));
        const fl energy = eval_interacting_pairs_deriv(tables, 1000, pairs, xyz, all);
        const fl excess = eval_interacting_pairs_deriv(tables, 1000, spurious, xyz, erroneous);
        for (const auto& s : serial) {
            const auto i = s.first;
            std::cout << path << '\t' << s.second << '\t' << energy << '\t' << excess;
            for (int axis = 0; axis < 3; ++axis) std::cout << '\t' << all[i][axis];
            for (int axis = 0; axis < 3; ++axis) std::cout << '\t' << erroneous[i][axis];
            std::cout << '\n';
        }
    } else {
        if (spurious.empty()) std::cout << path << '\t' << pairs.size() << "\t0\t\t\t\n";
        for (const auto& p : spurious)
            std::cout << path << '\t' << pairs.size() << '\t' << spurious.size() << '\t'
                      << serial.at(p.a) << '\t' << serial.at(p.b) << '\t' << distance[p.a][p.b] << '\n';
    }
    if (require_valid && !spurious.empty())
        throw std::runtime_error("nonbonded interaction spans at most three graph edges");
}
}

int main(int argc, char** argv) {
    const bool forces = argc > 1 && std::string(argv[1]) == "--forces";
    const bool require_valid = argc > 1 && std::string(argv[1]) == "--require-valid";
    const int first = forces || require_valid ? 2 : 1;
    if (argc <= first) return 2;
    std::cout << std::setprecision(17) << (forces
        ? "path\tserial\tenergy\texcess\tgx\tgy\tgz\texcess_gx\texcess_gy\texcess_gz\n"
        : "path\tpairs\tspurious\tserial_a\tserial_b\tgraph_distance\n");
    for (int k = first; k < argc; ++k) check(argv[k], forces, require_valid);
}
