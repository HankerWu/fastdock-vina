// Diagnose real nonbonded pairs across the intentionally cut chemical edge.
#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::cout << "ligand,listed_pairs,glue_pairs,serial_a,serial_b,closed_distance\n";
    for (int k = 1; k < argc; ++k) {
        const model m = parse_ligand_pdbqt_from_file(argv[k], atom_type::XS);
        const sz n = m.num_atoms();
        const sz infinity = n + 1;
        std::vector<szv> d(n, szv(n, infinity));
        std::array<szv, 4> ends;
        std::map<sz, unsigned> serial;
        for (const auto& entry : ligand_heavy_indices_by_serial(m)) serial[entry.second] = entry.first;
        for (sz i = 0; i < n; ++i) {
            d[i][i] = 0;
            const auto& a = m.get_atom(i);
            if (a.ad >= AD_TYPE_CG0 && a.ad <= AD_TYPE_CG3) ends[a.ad - AD_TYPE_CG0].push_back(i);
            for (const auto& b : a.bonds)
                if (!b.connected_atom_index.in_grid && !ad_is_closure_dummy(a.ad) &&
                    !ad_is_closure_dummy(m.get_atom(b.connected_atom_index).ad))
                    d[i][b.connected_atom_index.i] = 1;
        }
        for (const auto& pair : ends) {
            if (pair.empty()) continue;
            VINA_CHECK(pair.size() == 2);
            d[pair[0]][pair[1]] = d[pair[1]][pair[0]] = 1;
        }
        for (sz j = 0; j < n; ++j)
            for (sz i = 0; i < n; ++i)
                for (sz b = 0; b < n; ++b)
                    d[i][b] = std::min(d[i][b], d[i][j] + d[j][b]);
        const auto& pairs = m.get_ligand(0).pairs;
        const auto prefix = [&] {
            std::cout << argv[k] << ',' << pairs.size() << ',' << m.get_glue_pairs().size() << ',';
        };
        bool found = false;
        for (const auto& p : pairs) {
            if (d[p.a][p.b] > 3) continue;
            found = true;
            prefix();
            std::cout << serial.at(p.a) << ',' << serial.at(p.b) << ',' << d[p.a][p.b] << '\n';
        }
        if (!found) { prefix(); std::cout << ",,\n"; }
    }
}
