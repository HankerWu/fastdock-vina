// Export parsed heavy-atom bonds and scored intraligand pairs by PDBQT serial.
#include <algorithm>
#include <iostream>
#include <map>
#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    std::map<sz,unsigned> serial;
    for (const auto& kv : ligand_heavy_indices_by_serial(m)) serial[kv.second] = kv.first;
    const auto write = [&](const char* kind, const sz i, const sz j) {
        if (!serial.count(i) || !serial.count(j)) return;
        const auto ordered = std::minmax(serial.at(i), serial.at(j));
        std::cout << kind << '\t' << ordered.first << '\t' << ordered.second << '\n';
    };
    for (sz i = 0; i < m.num_atoms(); ++i)
        for (const bond& b : m.get_atom(i).bonds)
            if (!b.connected_atom_index.in_grid && i < b.connected_atom_index.i)
                write("bond", i, b.connected_atom_index.i);
    for (const auto& p : m.get_ligand(0).pairs) write("pair", p.a, p.b);
}
