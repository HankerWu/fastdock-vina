// Read the exact prepared ligand models used by the two production CLI paths.
#include "parse_pdbqt.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) return 2;
    model m = argc == 2 ? parse_ligand_pdbqt_from_file(argv[1], atom_type::XS)
                        : parse_ligand_pdbqt_with_bonds(argv[1], atom_type::XS, argv[2]);
    std::ifstream pdbqt(argv[1]);
    std::map<unsigned, vec> deposited;
    std::string line;
    while (std::getline(pdbqt, line)) {
        if (line.compare(0, 6, "ATOM  ") && line.compare(0, 6, "HETATM")) continue;
        const unsigned id = std::stoul(line.substr(6,5));
        const vec xyz(std::stod(line.substr(30,8)),std::stod(line.substr(38,8)),std::stod(line.substr(46,8)));
        if (!deposited.emplace(id,xyz).second) throw std::runtime_error("duplicate prepared serial");
    }
    if (deposited.size() != m.num_movable_atoms()) throw std::runtime_error("atom count mismatch");
    std::vector<unsigned> ids;
    for (sz i=0; i<m.num_movable_atoms(); ++i) {
        const vec& x=m.movable_coords(i);
        std::vector<unsigned> matches;
        for (const auto& entry : deposited)
            if (vec_distance_sqr(entry.second,x) < 1e-12) matches.push_back(entry.first);
        if (matches.size()!=1) throw std::runtime_error("ambiguous prepared coordinate to serial");
        ids.push_back(matches.front());
    }
    std::cout << "atom,atom_xs_type,neighbor_order,neighbor,rotatable\n";
    for (sz i=0;i<m.num_movable_atoms();++i) {
        const atom& a=m.get_atom(atom_index(i,false));
        for (sz n=0;n<a.bonds.size();++n) {
            const bond& b=a.bonds[n];
            if (b.connected_atom_index.in_grid) throw std::runtime_error("unexpected grid atom");
            std::cout << ids[i] << ',' << a.xs << ',' << n << ',' << ids.at(b.connected_atom_index.i)
                      << ',' << b.rotatable << '\n';
        }
    }
}
