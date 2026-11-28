#include "vina.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

// rigid receptor, original PDBQT topology, center xyz, then one single-model
// coordinate file per pose; use the identical public scoring API for each.
int main(int argc, char** argv) {
    if (argc < 7) return 2;
    try {
        Vina engine("vina", 1, 85231, 0);
        engine.set_receptor(argv[1]);
        engine.set_ligand_from_file(argv[2]);
        engine.compute_vina_maps(std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]),
                                 24, 24, 24, 0.375, false, true);
        std::cout << std::setprecision(16);
        int failed = 0;
        for (int i = 6; i < argc; ++i) {
            try {
                engine.set_pose_coordinates_for_scoring(argv[i]);
                auto e = engine.score(0.0);
                std::cout << "POSE " << i-6 << " " << e.at(0) << " " << e.at(1)+e.at(2)
                          << " " << e.at(3)+e.at(4)+e.at(5) << "\n";
            } catch (const std::exception& exc) {
                ++failed;
                std::cerr << "Pose " << i-6 << " failed: " << exc.what() << "\n";
            }
        }
        return failed ? 1 : 0;
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 1;
    }
}
