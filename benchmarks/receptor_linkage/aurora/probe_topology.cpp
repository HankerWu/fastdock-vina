#include "parse_pdbqt.h"
#include <algorithm>
#include <iostream>
#include <set>
#include <tuple>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) return 1;
    const auto model = argc == 2 ? parse_ligand_pdbqt_from_file(argv[1], atom_type::XS)
                                 : parse_ligand_pdbqt_with_bonds(argv[1], atom_type::XS, argv[2]);
    const auto& lig = model.get_ligand(0);
    std::vector<unsigned> ids(model.num_atoms());
    for (const auto& l : lig.cont)
        if (l.second) ids[*l.second] = std::stoul(l.first.substr(6, 5));
    std::set<std::pair<unsigned,unsigned>> bonds;
    for (sz i = 0; i < model.num_atoms(); ++i) {
        const auto& atom = model.get_atom(i);
        std::cout << "type " << ids[i] << ' ' << atom.xs << '\n';
        for (const auto& b : atom.bonds) {
            const auto j = ids[b.connected_atom_index.i];
            bonds.emplace(std::min(ids[i],j), std::max(ids[i],j));
        }
    }
    for (const auto& b : bonds) std::cout << "bond " << b.first << ' ' << b.second << '\n';
    std::set<std::tuple<unsigned,unsigned,sz>> pairs;
    for (const auto& p : lig.pairs)
        pairs.emplace(std::min(ids[p.a],ids[p.b]), std::max(ids[p.a],ids[p.b]), p.table_index);
    for (const auto& p : pairs)
        std::cout << "pair " << std::get<0>(p) << ' ' << std::get<1>(p) << ' ' << std::get<2>(p) << '\n';
}
