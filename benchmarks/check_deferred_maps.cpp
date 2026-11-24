#include <iostream>
#include <string>

#include "vina.h"

int main(int argc, char** argv) {
    if (argc != 6) return 2;
    const std::string receptor(argv[1]);
    const std::string ligand(argv[2]);
    const std::string flex_rigid(argv[3]);
    const std::string flex_residue(argv[4]);
    const std::string macrocycle(argv[5]);

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

            Vina batch(scoring, 1, 42, 0, no_refine);
            batch.set_receptor(receptor);
            batch.initialize_vina_batch_maps(15.19, 53.903, 16.917, 20, 20, 20, 0.375);
            for (const auto& input : {ligand, macrocycle, ligand}) {
                Vina single(scoring, 1, 42, 0, no_refine);
                single.set_receptor(receptor);
                single.set_ligand_from_file(input);
                single.compute_vina_maps(15.19, 53.903, 16.917, 20, 20, 20, 0.375);

                batch.set_ligand_from_file(input);
                if (batch.score(-3.25) != single.score(-3.25)
                        || batch.score() != single.score()
                        || batch.optimize(8) != single.optimize(8)) {
                    std::cerr << "Batch energies differ for " << scoring << ", no_refine="
                              << no_refine << ", input=" << input << '\n';
                    return 1;
                }
            }
        }
    }
    std::cout << "Deferred, eager and mixed batch maps agree for scoring and local search\n";
}
