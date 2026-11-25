#include <iostream>
#include <stdexcept>

#include "parse_pdbqt.h"

sz anchor_index(const model& receptor) {
    const auto& bonds = receptor.get_atom(atom_index(0, true)).bonds;
    if (bonds.size() != 1 || bonds.front().connected_atom_index.in_grid)
        throw std::runtime_error("Expected one bond from the grid atom to the residue anchor");
    return bonds.front().connected_atom_index.i;
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto receptor = parse_receptor_pdbqt(argv[1], argv[2]);
    const auto ligand = parse_ligand_pdbqt_from_file(argv[3], atom_type::XS);
    const sz original_anchor = anchor_index(receptor);
    auto working = receptor;
    for (sz count = 1; count <= 2; ++count) {
        working.append(ligand);
        if (anchor_index(working) != original_anchor + count * ligand.num_movable_atoms() ||
            anchor_index(receptor) != original_anchor) {
            std::cerr << "Appending a ligand must remap only the working model's receptor bond\n";
            return 1;
        }
        working.verify_bond_lengths();
    }
    std::cout << "Rigid-to-flex receptor bonds survive repeated ligand appends\n";
}
