#ifndef FASTDOCK_RESEARCH_CH_PI_SEARCH_GRID_H
#define FASTDOCK_RESEARCH_CH_PI_SEARCH_GRID_H

#include "grid.h"
#include "igrid.h"
#include "model.h"
#include "parse_pdbqt.h"
#include <set>
#include <sstream>

// A search-only grid for ligand CH -> receptor aromatic carbon. The Python
// graph reader supplies the typed serials and centers and the published kernel.
class ch_pi_search_grid : public igrid {
    const igrid& native;
    grid field;
    szv donors;

public:
    ch_pi_search_grid(const igrid& base, const grid_dims& dims, const model& m,
                      std::istream& input) : native(base), field(dims) {
        fl amplitude, center, width, cutoff;
        if (!(input >> amplitude >> center >> width >> cutoff) || !(amplitude >= 0) ||
            !(center > 0 && width > 0 && cutoff > 0) ||
            !std::isfinite(amplitude + center + width + cutoff))
            throw std::runtime_error("Invalid CH-pi kernel");
        const auto indices = ligand_heavy_indices_by_serial(m);
        sz serial;
        std::set<sz> seen;
        std::string line;
        std::getline(input, line);
        std::getline(input, line);
        std::istringstream row(line);
        while (row >> serial) {
            auto item = indices.find(serial);
            if (item == indices.end() || m.movable_atom(item->second).el != EL_TYPE_C ||
                !seen.insert(serial).second)
                throw std::runtime_error("Invalid CH-pi ligand carbon serial");
            donors.push_back(item->second);
        }
        if (!row.eof() || donors.empty()) throw std::runtime_error("Invalid CH-pi donor list");
        while (std::getline(input, line)) {
            vec atom;
            std::istringstream xyz(line);
            std::string extra;
            if (!(xyz >> atom[0] >> atom[1] >> atom[2]) || xyz >> extra ||
                !std::isfinite(atom[0]+atom[1]+atom[2]))
                throw std::runtime_error("Invalid aromatic carbon position");
            int lo[3], hi[3];
            for (sz d = 0; d < 3; ++d) {
                const fl spacing = field.m_factor_inv[d];
                lo[d] = std::max(0, int(std::ceil((atom[d] - cutoff - dims[d].begin) / spacing)));
                hi[d] = std::min(int(dims[d].n_voxels),
                                 int(std::floor((atom[d] + cutoff - dims[d].begin) / spacing)));
            }
            for (int k = lo[2]; k <= hi[2]; ++k)
                for (int j = lo[1]; j <= hi[1]; ++j)
                    for (int i = lo[0]; i <= hi[0]; ++i) {
                        const vec delta = field.index_to_argument(i,j,k) - atom;
                        const fl r2 = delta * delta;
                        if (r2 < cutoff * cutoff) {
                            const fl z = (std::sqrt(r2) - center) / width;
                            field.m_data(i,j,k) -= amplitude * std::exp(-0.5*z*z);
                        }
                    }
        }
    }

    fl at(const vec& x, vec* gradient = nullptr) const {
        // Native search supplies the wall penalty once; this additional
        // attractive field uses a clamped extension and contributes no wall.
        return gradient ? field.evaluate(x, 0, max_fl, *gradient)
                        : field.evaluate(x, 0, max_fl);
    }
    fl eval(const model& m, fl cap) const override {
        fl total = native.eval(m, cap);
        for (sz i : donors) total += at(m.coords[i]);
        return total;
    }
    fl eval_deriv(model& m, fl cap) const override {
        fl total = native.eval_deriv(m, cap);
        for (sz i : donors) {
            vec g;
            total += at(m.coords[i], &g);
            m.minus_forces[i] += g;
        }
        return total;
    }
    fl eval_intra(model& m, fl cap) const override { return native.eval_intra(m, cap); }
};
#endif
