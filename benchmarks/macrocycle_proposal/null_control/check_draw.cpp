// Check the actual mutation entry point against its unrestricted outcome and RNG advancement.
#include "fixed_ring.h"
#include "mutate.h"
#include "parse_pdbqt.h"

#include <algorithm>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const conf start = m.get_initial_conf();
    const fixed_ring ring(m);
    if (!ring.active()) return 3;
    const auto& choices = ring.allowed_moves();
    const auto& ligand = start.ligands.at(0);
    const sz n = 2 + ligand.torsions.size();
    unsigned unchanged = 0, legal = 0;
    for (unsigned seed = 1; seed<=1024; ++seed) {
        rng probe(seed), full(seed), restricted(seed);
        sz selected = random_int(0, n-1, probe);
        const bool allow = std::binary_search(choices.begin(),choices.end(),selected);
        conf original = start, trial = start;
        mutate_conf(original,m,2,full,false);
        mutate_conf(trial,m,2,restricted,false,false,false,false,&choices);
        if (full!=restricted) return 4;
        const auto& expected = allow ? original.ligands[0] : ligand;
        const auto& actual = trial.ligands[0];
        if (expected.torsions!=actual.torsions ||
            expected.rigid.orientation!=actual.rigid.orientation ||
            vec_distance_sqr(expected.rigid.position,actual.rigid.position)!=0) return 5;
        if (allow) ++legal; else ++unchanged;
    }
    std::cout << n << ',' << choices.size() << ',' << n-choices.size()
              << ',' << legal << ',' << unchanged << '\n';
}
