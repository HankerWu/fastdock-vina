// Compare the graph, branch edges and ordered pair tables under the actual AD typing.
#include "parse_pdbqt.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) return 2;
    model m = argc == 2 ? parse_ligand_pdbqt_from_file(argv[1], atom_type::AD)
                        : parse_ligand_pdbqt_with_bonds(argv[1], atom_type::AD, argv[2]);
    const auto ligand = m.get_ligand(0);
    std::vector<unsigned> serial(m.num_atoms());
    for (const auto& record : ligand.cont)
        if (record.second) serial.at(*record.second) = std::stoul(record.first.substr(6,5));
    for (sz i=0; i<m.num_atoms(); ++i) {
        auto atom = m.get_atom(i);
        std::cout << "atom " << serial[i] << ' ' << atom.ad << ' ' << atom.xs << ' '
                  << std::setprecision(17) << atom.charge << '\n';
        std::vector<std::pair<unsigned,bool>> neighbors;
        for (const auto& bond : atom.bonds) {
            if (bond.connected_atom_index.in_grid) throw std::runtime_error("unexpected grid atom");
            neighbors.emplace_back(serial.at(bond.connected_atom_index.i), bond.rotatable);
        }
        std::sort(neighbors.begin(),neighbors.end());
        for (const auto& edge : neighbors)
            std::cout << "edge " << serial[i] << ' ' << edge.first << ' ' << edge.second << '\n';
    }
    for (const auto& pair : ligand.pairs)
        std::cout << "pair " << serial.at(pair.a) << ' ' << serial.at(pair.b) << ' ' << pair.table_index << '\n';
}
