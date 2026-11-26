// Read the exact conformation-independent inputs used by this engine.
#include "conf_independent.h"
#include "parse_pdbqt.h"
#include <iostream>

int main(int argc, char** argv) {
    try {
        if (argc < 2) throw std::runtime_error("Supply prepared ligand PDBQT files");
        std::cout << "compound,num_tors,torsdof,num_rotors,heavy_atoms\n";
        for (int i = 1; i < argc; ++i) {
            const std::string path(argv[i]);
            const auto start = path.find_last_of("/\\") + 1;
            const auto name = path.substr(start, path.find_last_of('.') - start);
            const conf_independent_inputs inputs(parse_ligand_pdbqt_from_file(path, atom_type::XS));
            std::cout << name << ',' << inputs.num_tors << ',' << inputs.torsdof
                      << ',' << inputs.num_rotors << ',' << inputs.num_heavy_atoms << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
