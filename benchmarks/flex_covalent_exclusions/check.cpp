// For flex CG the path CG-CB-CA-C1 has three bonds. The next atom C2
// is four bonds away and must still contribute. The ligand sees both.
#include "cache.h"
#include "non_cache.h"
#include "parse_pdbqt.h"
#include "vina.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

static void near(fl measured, fl expected) {
    if (!std::isfinite(measured) || !std::isfinite(expected) || std::abs(measured - expected) > 1e-10)
        throw std::runtime_error("energy or gradient differs: " + std::to_string(measured)
                                 + " vs " + std::to_string(expected));
}

struct fields {
    model m;
    precalculate_byatom pairs;
    non_cache explicit_field;
    cache maps;

    fields(model input, const grid_dims& box, std::shared_ptr<const ScoringFunction> sf)
        : m(std::move(input)), explicit_field(m, box, sf, 100), maps(box,100) {
        m.append(parse_ligand_pdbqt_from_string(
            "ROOT\nATOM      1    C UNL     1      11.100  50.700  15.300  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n",
            atom_type::XS));
        pairs = precalculate_byatom(*sf,m);
        explicit_field.prepare_for_model(m);
        explicit_field.populate_maps(maps, m, cache::map_types(m.get_movable_atom_types(atom_type::XS)), true);
    }
};

int main(int argc,char** argv) {
    if (argc != 5) return 2;
    auto whole = parse_receptor_pdbqt(argv[1],argv[3]);
    auto far = parse_receptor_pdbqt(argv[2],argv[3]);
    const auto excluded = whole.rigid_bond_exclusions();
    if (excluded.size() != 1 || excluded[0] != szv{0} || !far.rigid_bond_exclusions().empty())
        throw std::runtime_error("incorrect covalent neighborhood");
    auto sf = std::make_shared<ScoringFunction>(SF_VINA,
        flv{-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50,5*0.05846/0.1-1});
    grid_dims box;
    const vec center(11,50,14);
    for (sz k=0; k<3; ++k) {
        box[k].begin=center[k]-5;
        box[k].end=center[k]+5;
        box[k].n_voxels=32;
    }
    fields all(std::move(whole),box,sf), remaining(std::move(far),box,sf);
    cache::cubic_view cubic_all(all.maps), cubic_remaining(remaining.maps);
    const igrid* field1[] = {&all.explicit_field, &all.maps, &cubic_all};
    const igrid* field2[] = {&remaining.explicit_field, &remaining.maps, &cubic_remaining};
    std::cout << std::setprecision(12) << "field,angle,cap,flex_energy,flex_gradient,ligand_contact_change\n";
    for (fl angle : {0.,0.71,-1.29}) {
        conf pose=all.m.get_initial_conf();
        pose.flex.at(0).torsions.at(0)=angle;
        all.m.set(pose);
        remaining.m.set(pose);
        for (fl cap : {0.1,1000.}) {
            for (sz k=0;k<3;++k) {
                const auto& a=*field1[k];
                const auto& b=*field2[k];
                const fl intra=a.eval_intra(all.m,cap);
                near(intra,b.eval_intra(remaining.m,cap));
                if (!std::isfinite(intra) || std::abs(intra)<1e-6) throw std::runtime_error("lost four-bond neighbor");
                const fl full_delta=a.eval(all.m,cap)-b.eval(remaining.m,cap);
                if (!std::isfinite(full_delta) || std::abs(full_delta)<1e-6) throw std::runtime_error("lost ordinary ligand contact");
                flv gradient(all.m.get_size().num_degrees_of_freedom()), reference(gradient.size());
                const vec caps(cap,cap,cap);
                all.m.eval_deriv(all.pairs,a,caps,gradient);
                remaining.m.eval_deriv(remaining.pairs,b,caps,reference);
                near(gradient.back(),reference.back());
                std::cout << k << ',' << angle << ',' << cap << ',' << intra << ','
                          << gradient.back() << ',' << full_delta << '\n';
            }
        }
    }

    // Real public API: defer creation, then append a ligand of a previously
    // absent type. Neither appending a second ligand nor loading a map may
    // attach the flexible atom's private field to a ligand atom.
    const std::string carbon =
        "ROOT\nATOM      1    C UNL     1      11.100  50.700  15.300  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n";
    const std::string oxygen =
        "ROOT\nATOM      1    O UNL     1      10.000  47.000  16.000  1.00  0.00     0.000 OA\nENDROOT\nTORSDOF 0\n";
    {
        ifile source(argv[3]);
        ofile different(std::string(argv[4])+"-different.pdbqt");
        std::string line;
        while (std::getline(source,line)) {
            const auto coord=line.find("51.500");
            if (coord!=std::string::npos) line.replace(coord,6,"51.600");
            different << line << '\n';
        }
    }
    for (const std::string scheme : {"vina","vinardo"}) {
        const std::string prefix = std::string(argv[4]) + "-" + scheme;
        Vina cached(scheme,1,921,0,true), fresh(scheme,1,921,0,true), imported(scheme,1,921,0,true);
        cached.set_receptor(argv[1],argv[3]);
        fresh.set_receptor(argv[1],argv[3]);
        imported.set_receptor("",argv[3]);
        cached.set_ligand_from_string(carbon);
        fresh.set_ligand_from_string(std::vector<std::string>{carbon,oxygen});
        cached.compute_vina_maps(11,50,14,10,10,10,.3125,true,true);
        cached.score();
        cached.set_ligand_from_string(std::vector<std::string>{carbon,oxygen});
        fresh.compute_vina_maps(11,50,14,10,10,10,.3125,true,false);
        const auto expected=fresh.score(), reused=cached.score();
        for (sz i=0;i<expected.size();++i) near(expected[i],reused[i]);
        cached.write_maps(prefix);
        imported.load_maps(prefix);
        imported.set_ligand_from_string(std::vector<std::string>{carbon,oxygen});
        const auto recovered=imported.score();
        for (sz i=0;i<expected.size();++i)
            if (!std::isfinite(recovered[i]) || std::abs(expected[i]-recovered[i])>0.001)
                throw std::runtime_error("map round-trip changed score beyond map text precision");
        bool mismatch_rejected=false;
        try {
            cache mismatch;
            const auto different = parse_receptor_pdbqt("",std::string(argv[4])+"-different.pdbqt");
            mismatch.read(prefix,&different);
        } catch (const std::runtime_error& err) {
            mismatch_rejected=std::string(err.what()).find("Flexible receptor does not match")!=std::string::npos;
        }
        if (!mismatch_rejected) throw std::runtime_error("supplemental fields accepted the wrong flex input");
    }
}
