// Export native XS typing, heavy torsion count and mobile pair identities. No scoring here.
#include "parse_pdbqt.h"
#include <iomanip>
#include <iostream>
#include <map>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::cout << std::setprecision(17);
    const std::string kind(argv[1]);
    if (kind == "receptor") {
        const model m = parse_receptor_pdbqt(argv[2]);
        for (const auto& a : m.grid_atoms())
            if (a.xs < XS_TYPE_SIZE)
                std::cout << "R " << a.xs << ' ' << a.coords[0] << ' ' << a.coords[1] << ' ' << a.coords[2] << '\n';
    } else if (kind == "ligand") {
        const model m = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
        const auto& ligand = m.get_ligand(0);
        std::cout << "T " << ligand.active_heavy_torsions << '\n';
        std::map<sz, unsigned> serial;
        for (const auto& entry : ligand_heavy_indices_by_serial(m)) {
            serial[entry.second] = entry.first;
            std::cout << "L " << entry.first << ' ' << m.get_atom(entry.second).xs << '\n';
        }
        for (const auto& pair : ligand.pairs)
            std::cout << "P " << serial.at(pair.a) << ' ' << serial.at(pair.b) << '\n';
    } else {
        return 2;
    }
}
