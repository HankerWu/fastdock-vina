// Fixed-pose rigid configurational integral, using the production Vina atom
// typing, conformational normalization, pair tables, cutoff and per-atom curl.
#include "parse_pdbqt.h"
#include "cache.h"
#include "precalculate.h"
#include "curl.h"
#include "hinge.h"
#include <chrono>
#include <iomanip>
#include <iostream>

struct transform {
    vec translation;
    fl matrix[9];
    fl hinge;
    vec rotate(const vec &v) const {
        return vec(matrix[0]*v[0]+matrix[1]*v[1]+matrix[2]*v[2],
                   matrix[3]*v[0]+matrix[4]*v[1]+matrix[5]*v[2],
                   matrix[6]*v[0]+matrix[7]*v[1]+matrix[8]*v[2]);
    }
};

struct contact {
    vec xyz;
    const precalculate_element *table;
};
struct probe_atom {
    sz index;
    vec arm;
    std::vector<contact> neighbors;
    fl energy(const vec &pos, const precalculate &tables) const {
        fl e = 0;
        for (const contact &p : neighbors) {
            const fl r2 = vec_distance_sqr(p.xyz, pos);
            if (r2 < tables.cutoff_sqr())
                e += p.table->fast_sample(tables.fast_index(r2));
        }
        curl(e, 1000);
        return e;
    }
};

static std::vector<transform> read_nodes(const char *file) {
    std::ifstream in(file);
    std::vector<transform> nodes;
    transform t;
    while (in >> t.translation[0] >> t.translation[1] >> t.translation[2]) {
        for (fl &v : t.matrix)
            if (!(in >> v) || !std::isfinite(v)) throw std::runtime_error("Invalid rotation");
        if (!(in >> t.hinge) || !std::isfinite(t.hinge) || t.hinge<0 || t.hinge>=1) throw std::runtime_error("Invalid hinge node");
        nodes.push_back(t);
    }
    if (!in.eof() || nodes.empty()) throw std::runtime_error("Invalid quadrature nodes");
    return nodes;
}
static std::string top_pose(const char *file) {
    std::ifstream in(file);
    std::string line, block;
    while (std::getline(in, line)) {
        block += line + "\n";
        if (line == "ENDMDL") return block;
    }
    throw std::runtime_error("No complete first pose");
}
struct accumulator {
    fl log_sum = -std::numeric_limits<fl>::infinity(), minimum = 0;
    void add(fl energy, fl rt) {
        minimum = std::min(minimum, energy);
        const fl value = -energy / rt, maximum = std::max(value, log_sum);
        log_sum = maximum + std::log1p(std::exp(std::min(value, log_sum) - maximum));
    }
    fl work(sz n, fl rt) const { return -rt * (log_sum - std::log(fl(n))); }
};
int main(int argc, char **argv) {
    try {
        if (argc != 6 && argc != 7) throw std::runtime_error("usage: integrate receptor ligand pose_pool nodes temperature_K [torsion_bins]");
        const auto began = std::chrono::steady_clock::now();
        const fl rt = 0.00198720425864083 * std::stod(argv[5]);
        if (!std::isfinite(rt) || rt<=0) throw std::runtime_error("Invalid temperature");
        const auto nodes = read_nodes(argv[4]);
        const model receptor = parse_receptor_pdbqt(argv[1]);
        model ligand = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
        replace_fixed_pose_coordinates(ligand, top_pose(argv[3]));
        const ScoringFunction sf(SF_VINA, {-0.035579,-0.005156,0.840245,-0.035069,-0.587439,50,5*0.05846/0.1-1});
        const fl scale = sf.conf_independent(ligand, 1) - sf.conf_independent(ligand, 0);
        const vec center = ligand.ligand_heavy_centroid(0);
        std::bitset<XS_TYPE_SIZE> rtps, ltps;
        for (const atom &a : receptor.grid_atoms())
            if (a.xs < XS_TYPE_SIZE) rtps.set(a.xs);
        for (sz i=0; i<ligand.num_movable_atoms(); ++i) {
            const atom_base &a=ligand.movable_atom(i);
            if (a.ad == AD_TYPE_W)
                throw std::runtime_error("Hydrated-ligand scoring is unsupported");
            if (a.xs < XS_TYPE_SIZE && !is_glue_type(a.xs)) ltps.set(a.xs);
        }
        precalculate tables(sf.get_cutoff());
        tables.populate(sf, rtps, ltps);
        const unsigned bins=argc==7 ? std::stoul(argv[6]) : 16;
        if (!bins || bins>4096) throw std::runtime_error("Invalid torsion bin count");
        const local_hinge hinge(ligand,sf,bins);
        std::vector<probe_atom> atoms;
        for (sz i=0; i<ligand.num_movable_atoms(); ++i) {
            const auto &a = ligand.movable_atom(i);
            if (a.xs >= XS_TYPE_SIZE || is_glue_type(a.xs)) continue;
            const vec xyz = ligand.movable_coords(i);
            probe_atom probe{i,xyz-center, {}};
            fl max_shift = 0;
            for (const transform &t : nodes) {
                const vec arm=hinge.states[sz(t.hinge*hinge.states.size())].xyz[i]-center;
                max_shift=std::max({max_shift,(t.rotate(probe.arm)+t.translation-probe.arm).norm(),
                    (arm-probe.arm).norm(), (t.rotate(arm)+t.translation-probe.arm).norm()});
            }
            const fl radius2 = sqr(sf.get_cutoff()+max_shift+1e-8);
            for (const atom &b : receptor.grid_atoms()) {
                if (b.xs >= XS_TYPE_SIZE || vec_distance_sqr(xyz,b.coords) > radius2) continue;
                const sz idx=get_type_pair_index(atom_type::XS,a,b);
                probe.neighbors.push_back({b.coords,&tables.element(idx)});
            }
            atoms.push_back(std::move(probe));
        }
        // Atom-major accumulation keeps its exact, displacement-bounded neighbor list hot.
        fl e0=0;
        std::vector<std::array<fl,3>> energies(nodes.size(), {0,0,0});
        for (const probe_atom &a : atoms) {
            const vec starting=center+a.arm;
            e0+=a.energy(starting,tables);
            flv internal_only;
            for (const auto &state : hinge.states)
                internal_only.push_back(a.energy(state.xyz[a.index],tables));
            for (sz n=0; n<nodes.size(); ++n) {
                const transform &t=nodes[n];
                const vec rotated=center+t.rotate(a.arm);
                const sz bin=sz(t.hinge*hinge.states.size());
                const vec bent=hinge.states[bin].xyz[a.index];
                const fl rigid_energy=a.energy(rotated+t.translation,tables);
                energies[n][0]+=hinge.states.size()==1 ? rigid_energy : a.energy(center+t.rotate(bent-center)+t.translation,tables);
                energies[n][1]+=rigid_energy;
                energies[n][2]+=internal_only[bin];
            }
        }
        accumulator accum[4];
        std::cout << std::setprecision(11) << "nodes,n_heavy,inter,normalization,joint,rigid,hinge,no_free,min_delta,free,hinge_moment,first_index,second_index,wall_s\n";
        for (sz n=0; n<nodes.size(); ++n) {
            const fl di=hinge.states[sz(nodes[n].hinge*hinge.states.size())].intra_delta;
            accum[0].add(scale*(energies[n][0]+di-e0),rt);
            accum[1].add(scale*(energies[n][1]-e0),rt);
            accum[2].add(scale*(energies[n][2]+di-e0),rt);
            accum[3].add(scale*di,rt);
            const sz count=n+1;
            if (!(count & (count-1)) && count>=256) {
                const fl wall=std::chrono::duration<fl>(std::chrono::steady_clock::now()-began).count();
                std::cout << count << ',' << atoms.size() << ',' << e0 << ',' << scale;
                const fl free=accum[3].work(count,rt);
                std::cout << ',' << accum[0].work(count,rt)-free << ',' << accum[1].work(count,rt)
                    << ',' << accum[2].work(count,rt)-free << ',' << accum[0].work(count,rt);
                std::cout << ',' << accum[0].minimum << ',' << free << ',' << hinge.reduced_moment
                    << ',' << hinge.start << ',' << hinge.end << ',' << wall << '\n';
            }
        }
    } catch(const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
