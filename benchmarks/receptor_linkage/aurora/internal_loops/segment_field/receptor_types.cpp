// Dump the actual XS types after the production rigid-receptor bond inference.
#include <iomanip>
#include <iostream>

#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: receptor_types receptor.pdbqt\n";
        return 2;
    }
    const model receptor = parse_receptor_pdbqt(argv[1]);
    std::cout << std::setprecision(17);
    for (const auto& atom : receptor.grid_atoms())
        std::cout << atom.coords[0] << ' ' << atom.coords[1] << ' '
                  << atom.coords[2] << ' ' << atom.xs << ' '
                  << (atom.xs < XS_TYPE_SIZE) << '\n';
}
