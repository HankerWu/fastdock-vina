#include <iostream>
#include <string>

#include "vina.h"

int main(int argc, char** argv) {
    if (argc != 5) return 2;
    const std::string receptor(argv[1]);
    const std::string ligand(argv[2]);
    const std::string flex_rigid(argv[3]);
    const std::string flex_residue(argv[4]);

    for (const std::string scoring : {"vina", "vinardo"}) {
        for (bool no_refine : {false, true}) {
            for (bool flexible : {false, true}) {
                Vina eager(scoring, 1, 42, 0, no_refine);
                Vina deferred(scoring, 1, 42, 0, no_refine);
                for (Vina* engine : {&eager, &deferred}) {
                    engine->set_receptor(flexible ? flex_rigid : receptor,
                                         flexible ? flex_residue : "");
                    engine->set_ligand_from_file(ligand);
                    engine->compute_vina_maps(15.19, 53.903, 16.917, 20, 20, 20,
                                              0.375, false, engine == &deferred);
                }
                if (eager.score(-3.25) != deferred.score(-3.25)
                        || eager.score() != deferred.score()
                        || eager.optimize(8) != deferred.optimize(8)) {
                    std::cerr << "Deferred map energies differ for " << scoring
                              << ", no_refine=" << no_refine << ", flexible=" << flexible << '\n';
                    return 1;
                }
            }
        }
    }
    std::cout << "Deferred and eager maps agree for scoring and local search\n";
}
