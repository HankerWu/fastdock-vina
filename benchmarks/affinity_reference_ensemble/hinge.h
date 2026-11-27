#ifndef FASTDOCK_BENCH_HINGE_H
#define FASTDOCK_BENCH_HINGE_H
#include "precalculate.h"
#include <map>

// Extract declared torsional cuts from the prepared PDBQT tree. Atom indices
// are the engine's internal topology, never serial numbers inferred from XYZ.
struct torsion_cut {
    sz start, end;
    std::vector<bool> moving;
};
static std::vector<torsion_cut> cuts(const model &m) {
    std::map<unsigned,sz> index;
    std::vector<torsion_cut> result;
    szv open;
    const ligand lig=m.get_ligand(0);
    for (const auto &entry : lig.cont) {
        const std::string &line=entry.first;
        if (entry.second) {
            const auto i=*entry.second;
            if (!index.emplace(std::stoul(line.substr(6,5)),i).second)
                throw std::runtime_error("Duplicate atom serial");
            for (sz j : open) result[j].moving[i]=true;
        } else if (line.compare(0,7,"BRANCH ")==0) {
            std::istringstream fields(line.substr(7));
            sz first,second;
            if (!(fields>>first>>second)) throw std::runtime_error("Invalid BRANCH");
            open.push_back(result.size());
            result.push_back({first,second,std::vector<bool>(m.num_movable_atoms())});
        } else if (line.compare(0,9,"ENDBRANCH")==0) {
            if (open.empty()) throw std::runtime_error("Invalid ENDBRANCH");
            open.pop_back();
        }
    }
    if (!open.empty()) throw std::runtime_error("Unclosed BRANCH");
    for (auto &cut : result) {
        cut.start=index.at(cut.start);
        cut.end=index.at(cut.end);
    }
    return result;
}
struct hinge_state {
    vecv xyz;
    fl intra_delta;
};
struct local_hinge {
    std::vector<hinge_state> states;
    fl reduced_moment=0;
    int start=-1,end=-1;
    explicit local_hinge(const model &m,const ScoringFunction &sf,unsigned bins) {
        vecv origin;
        bool macro=false;
        for (sz i=0;i<m.num_movable_atoms();++i) {
            origin.push_back(m.movable_coords(i));
            macro |= is_glue_type(m.movable_atom(i).xs);
        }
        if (!macro) {
            const auto options=cuts(m);
            const torsion_cut *best=nullptr;
            fl best_left=0,best_right=0;
            for (const auto &cut:options) {
                const vec axis=unit(origin[cut.end]-origin[cut.start]);
                fl moment[2]={0,0};
                for (sz i=0;i<origin.size();++i) {
                    if (m.movable_atom(i).is_hydrogen()) continue;
                    moment[cut.moving[i]]+=sqr(cross_product(axis,origin[i]-origin[cut.start]));
                }
                const fl total=moment[0]+moment[1];
                const fl reduced=total>0 ? moment[0]*moment[1]/total : 0;
                if (reduced>1e-8 && reduced>reduced_moment) {
                    reduced_moment=reduced;
                    best=&cut;best_left=moment[0];best_right=moment[1];
                }
            }
            if (best) {
                start=best->start;end=best->end;
                const vec axis=unit(origin[end]-origin[start]);
                const vec center=m.ligand_heavy_centroid(0);
                const auto pairs=m.get_ligand(0).pairs;
                const precalculate_byatom tables(sf,m);
                const fl i0=eval_interacting_pairs(tables,1000,pairs,origin);
                for (unsigned k=0;k<bins;++k) {
                    const fl theta=pi/6*(2*(k+0.5)/bins-1);
                    vecv next;
                    vec centroid(0,0,0);sz heavy=0;
                    for (sz i=0;i<origin.size();++i) {
                        const fl fraction=(best->moving[i]? best_left : -best_right)/(best_left+best_right);
                        const fl angle=fraction*theta;
                        const vec v=origin[i]-origin[start];
                        const vec xyz=origin[start]+std::cos(angle)*v+std::sin(angle)*cross_product(axis,v)
                            +((1-std::cos(angle))*(axis*v))*axis;
                        next.push_back(xyz);
                        if (!m.movable_atom(i).is_hydrogen()) {centroid+=xyz;++heavy;}
                    }
                    const vec shift=center-(1.0/heavy)*centroid;
                    for (vec &x:next) x+=shift;
                    const fl intra=eval_interacting_pairs(tables,1000,pairs,next);
                    states.push_back({std::move(next),intra-i0});
                }
                return;
            }
        }
        states.push_back({std::move(origin),0});
    }
private:
    static vec unit(const vec &v) {
        if (v.norm()<1e-8) throw std::runtime_error("Degenerate rotation axis");
        return (1/v.norm())*v;
    }
};
#endif
