// Exercise the extended XS glue objective around the final table bins using
// the actual tabulated pair type and post-curl Cartesian derivative.
#include <iomanip>
#include <iostream>
#include "curl.h"
#include "parse_pdbqt.h"
#include "precalculate.h"

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
    const auto& pairs = m.get_glue_pairs();
    VINA_CHECK(!pairs.empty());
    std::cout << "weight,cap,r,energy,radial_grad,numeric_grad,error\n";
    for (fl w : {0., 37., 50.}) {
        const ScoringFunction sf(SF_VINA, flv{-0.035579,-0.005156,0.840245,-0.035069,-0.587439,w,0.05846});
        const precalculate_byatom p(sf, m);
        for (fl cap : {10.,1000.}) {
            auto evaluate = [&](fl r) {
                auto ed = p.glue_deriv(pairs.front(), r*r);
                fl g = ed.second*r;
                curl(ed.first,g,cap);
                return pr(ed.first,g);
            };
            const fl join = std::sqrt((std::ceil(32*sqr(macrocycle_glue_cutoff))-2)/32);
            for (fl r : {19.5, join-.0001,join,join+.0001,20.,23.,30.}) {
                const auto ed = evaluate(r);
                const fl h = 1e-6;
                const fl numerical = (evaluate(r+h).first-evaluate(r-h).first)/(2*h);
                const fl error = std::abs(ed.second-numerical);
                VINA_CHECK(error < 0.0006 && std::isfinite(ed.first));
                if (w > 0) VINA_CHECK(ed.first > 0 && ed.second > 0);
                std::cout << std::setprecision(12) << w << ',' << cap << ',' << r << ',' << ed.first
                          << ',' << ed.second << ',' << numerical << ',' << error << '\n';
            }
        }
    }
}
