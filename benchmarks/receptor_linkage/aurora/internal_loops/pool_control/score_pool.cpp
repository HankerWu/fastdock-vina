// Read a sequence of frozen poses without re-inferring the prepared covalent graph.
#include "vina.h"
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
    if (argc != 6 && argc != 7) {
        std::cerr << "score_pool RECEPTOR LIGAND CX CY CZ [BONDS] < paths.txt\n";
        return 2;
    }
    Vina engine("vina", 1, 1, 0, false);
    engine.set_receptor(argv[1]);
    if (argc == 7)
        engine.set_ligand_from_file_with_bonds(argv[2], argv[6]);
    else
        engine.set_ligand_from_file(argv[2]);
    engine.compute_vina_maps(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]),
                             26., 26., 26., 0.375, false, true);
    std::cout << std::setprecision(17);
    std::string path;
    while (std::getline(std::cin, path)) {
        engine.set_pose_coordinates_for_scoring(path);
        const auto energies = engine.score(0.);
        // Total, explicit intermolecular energy, ligand intramolecular energy.
        std::cout << energies.at(0) << ',' << energies.at(1) << ',' << energies.at(5) << '\n';
    }
}
