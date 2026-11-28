// Exercise both public fixed-coordinate score overloads across setup stages.
#include "vina.h"
#include <iostream>
#include <limits>
#include <stdexcept>

template<class F>
void require_rejection(F action, const std::string& expected) {
    try {
        action();
    } catch (const vina_runtime_error& e) {
        if (std::string(e.what()).find(expected) != std::string::npos) return;
        throw;
    }
    throw std::runtime_error("Missing rejection: " + expected);
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    Vina v("vina", 1, 20261128, 0);
    require_rejection([&]{ v.score(0); }, "not initialized");
    v.set_receptor(argv[1]);
    v.set_ligand_from_file(argv[2]);
    require_rejection([&]{ v.score(); }, "Affinity maps");
    require_rejection([&]{ v.score(0); }, "Affinity maps");
    v.compute_vina_maps(15.190, 53.903, 16.917, 30, 30, 30, 0.375, false, true);
    const auto automatic = v.score();
    if (v.score(automatic.at(7)) != automatic)
        throw std::runtime_error("Reusing the native reference changed scoring.");
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                            std::numeric_limits<double>::infinity(),
                            -std::numeric_limits<double>::infinity()})
        require_rejection([&]{ v.score(invalid); }, "must be finite");
    v.compute_vina_maps(55.190, 53.903, 16.917, 30, 30, 30, 0.375, false, true);
    require_rejection([&]{ v.score(); }, "outside the grid box");
    require_rejection([&]{ v.score(0); }, "outside the grid box");
    std::cout << "Public API setup, native reference, nonfinite and box checks passed.\n";
}
