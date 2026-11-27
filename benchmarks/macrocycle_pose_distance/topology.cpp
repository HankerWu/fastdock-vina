// Integer-only signature of the inferred native scoring graph and its ordered tables.
#include "parse_pdbqt.h"
#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const auto lig = m.get_ligand(0);
    std::cout << "torsions " << lig.active_heavy_torsions << ' ' << lig.degrees_of_freedom << '\n';
    for (const auto& line : lig.cont) {
        if (line.second)
            std::cout << "serial " << *line.second << ' ' << line.first.substr(6,5) << '\n';
        else if (line.first.compare(0, 6, "BRANCH") == 0 || line.first.compare(0, 9, "ENDBRANCH") == 0)
            std::cout << line.first << '\n';
    }
    for (sz i = 0; i < m.num_atoms(); ++i) {
        const auto atom = m.get_atom(i);
        std::cout << "atom " << i << ' ' << atom.ad << ' ' << atom.xs << ' ' << atom.charge << '\n';
        std::vector<std::pair<sz, bool>> edges;
        for (const auto& bond : atom.bonds)
            if (!bond.connected_atom_index.in_grid && i < bond.connected_atom_index.i)
                edges.emplace_back(bond.connected_atom_index.i, bond.rotatable);
        std::sort(edges.begin(),edges.end());
        for (const auto& edge : edges)
            std::cout << "bond " << i << ' ' << edge.first << ' ' << edge.second << '\n';
    }
    const auto table = [](const char* label, const interacting_pairs& pairs) {
        for (const auto& p : pairs)
            std::cout << label << ' ' << p.a << ' ' << p.b << ' ' << p.table_index << '\n';
    };
    table("internal", lig.pairs);
    table("glue", m.get_glue_pairs());
}
