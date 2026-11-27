// Exercise declared rotors that fail both the local cutoff and the bead search.
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "parse_pdbqt.h"

void verify(const model& m, fl wanted_length) {
    sz count = 0;
    for (sz i = 0; i < m.num_atoms(); ++i) {
        const auto atom = m.get_atom(i);
        for (const auto& bond : atom.bonds) {
            if (!bond.rotatable) continue;
            if (bond.connected_atom_index.in_grid || std::abs(bond.length - wanted_length) > 1e-8)
                throw std::runtime_error("Incorrect atom or changed rotor length");
            ++count;
        }
    }
    if (count != 2) throw std::runtime_error("Expected exactly one bidirectional rotor bond");
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    for (auto type : {atom_type::XS, atom_type::AD}) {
        for (const char* x : {"   1.500", "   2.100", "  35.000"}) {
            const fl distance = std::stod(x);
            const std::string atom = std::string("ATOM      2  C   UNL     1    ") +
                                     x + "   0.000   0.000  1.00  0.00     0.000 C\n";
            const std::string ligand = "ROOT\n"
                "ATOM      1  C   UNL     1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
                "ENDROOT\nBRANCH 1 2\n" + atom +
                "ATOM      3  C   UNL     1    " + x + "   1.500   0.000  1.00  0.00     0.000 C\n"
                "ENDBRANCH 1 2\nTORSDOF 1\n";
            verify(parse_ligand_pdbqt_from_string(ligand, type), distance);
        }
        verify(parse_receptor_pdbqt("", argv[1], type), 2.1);
    }
    std::cout << "Declared ligand and flex-residue rotor adjacency retained for XS and AD.\n";
}
