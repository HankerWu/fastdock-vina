// The integer suffix on a PDBQT ring-closure pair identifies the broken bond
// within one molecule. It must not change its chemistry or couple two molecules.
#include <fstream>
#include <iostream>
#include <sstream>

#include "non_cache.h"
#include "parse_pdbqt.h"

namespace {
std::string relabel(std::string pdbqt, char label) {
    size_t pos = 0;
    while ((pos = pdbqt.find("G0", pos)) != std::string::npos) {
        pdbqt[pos + 1] = label;
        pos += 2;
    }
    return pdbqt;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::ifstream input(argv[2]);
    std::ostringstream content;
    content << input.rdbuf();
    require(!content.str().empty(), "Missing cyclododecane fixture");

    for (sz tag = XS_TYPE_C_H_CG0; tag <= XS_TYPE_C_H_CG3; tag += 3) {
        require(xs_radius(tag) == xs_radius(XS_TYPE_C_H) &&
                xs_radius(tag + 1) == xs_radius(XS_TYPE_C_P) && xs_radius(tag + 2) == 0,
                "Incorrect Vina alias radius");
        require(xs_vinardo_radius(tag) == xs_vinardo_radius(XS_TYPE_C_H) &&
                xs_vinardo_radius(tag + 1) == xs_vinardo_radius(XS_TYPE_C_P) &&
                xs_vinardo_radius(tag + 2) == 0, "Incorrect Vinardo alias radius");
    }
    const auto sf = std::make_shared<ScoringFunction>(
        SF_VINA, flv{-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50,0.05846});
    grid_dims box;
    for (auto& d : box) {
        d.begin = -70; d.end = 70; d.n_voxels = 100;
    }
    const model receptor = parse_receptor_pdbqt(argv[1]);
    // Evaluate all suffixes in the same Cartesian poses, including offset
    // between ligands and independent torsion deformations within each ring.
    fl reference_energy = 0;
    flv reference_gradient;
    for (char label = '0'; label <= '3'; ++label) {
        const auto tagged = relabel(content.str(), label);
        for (auto typing : {atom_type::XS, atom_type::AD}) {
            const model first = parse_ligand_pdbqt_from_string(content.str(), typing);
            const model second = parse_ligand_pdbqt_from_string(tagged, typing);
            model combined = first;
            combined.append(second);
            require(first.get_glue_pairs().size() == 2 && second.get_glue_pairs().size() == 2,
                    "Each ring should retain two opposing closure restraints");
            require(combined.get_glue_pairs().size() == 4, "Intermolecular glue pairs");
            for (const auto& bond : combined.get_glue_pairs())
                require(combined.find_ligand(bond.a) == combined.find_ligand(bond.b),
                        "Ring closure crossed ligand boundary");
            require(combined.get_inter_pairs().size() == 12 * 12, "Lost carbon-carbon interactions");
            for (const auto& pair : combined.get_inter_pairs())
                require(combined.get_atom(pair.a).el == EL_TYPE_C && combined.get_atom(pair.b).el == EL_TYPE_C,
                        "Dummy in intermolecular pair list");
        }

        model m = receptor;
        m.append(parse_ligand_pdbqt_from_string(content.str(), atom_type::XS));
        m.append(parse_ligand_pdbqt_from_string(tagged, atom_type::XS));
        for (sz i = 0; i < m.num_atoms(); ++i) {
            const atom a = m.get_atom(i);
            if (a.el != EL_TYPE_C) continue;
            require(xs_radius(a.xs) == xs_radius(XS_TYPE_C_H) &&
                    xs_vinardo_radius(a.xs) == xs_vinardo_radius(XS_TYPE_C_H),
                    "Closure carbon radius depends on suffix");
        }

        auto c = m.get_initial_conf();
        c.ligands[0].torsions.front() += 0.3;
        c.ligands[1].torsions.back() -= 0.4;
        c.ligands[1].rigid.position += vec(9, 0, 0);
        m.set(c);
        non_cache field(m, box, sf, 100);
        field.prepare_for_model(m);
        const precalculate_byatom pairs(*sf, m);
        flv gradient(m.get_size().num_degrees_of_freedom());
        const fl energy = m.eval_deriv(pairs, field, vec(1000,1000,1000), gradient);
        if (label == '0') {
            reference_energy = energy;
            reference_gradient = gradient;
        }
        require(energy == reference_energy && gradient == reference_gradient,
                "Closure labels changed the potential or gradient at fixed coordinates");
    }
    std::cout << "XS/AD intermolecular pair topology, Vina/Vinardo radii, "
                 "and Vina energy/gradient label invariance passed.\n";
}
