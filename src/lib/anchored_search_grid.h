#ifndef VINA_ANCHORED_SEARCH_GRID_H
#define VINA_ANCHORED_SEARCH_GRID_H

#include <set>

#include "igrid.h"
#include "parse_pdbqt.h"
#include "model.h"

struct search_anchor_site {
	sz index;
	vec target;
	fl radius;
};

inline std::vector<search_anchor_site> read_search_anchors(std::istream& file, const model& ligand,
                                                       const vec& lower, const vec& upper) {
	std::vector<search_anchor_site> sites;
	std::set<sz> assigned;
	const auto heavy = ligand_heavy_indices_by_serial(ligand);
	std::string line;
	while (std::getline(file, line)) {
		if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
		std::istringstream row(line);
		std::string serial_text;
		fl x, y, z, radius = 0;
		std::string extra;
		if (!(row >> serial_text >> x >> y >> z) || serial_text.size() > 5 ||
		    serial_text.find_first_not_of("0123456789") != std::string::npos)
			throw std::invalid_argument("search_anchor rows require a positive atom serial and three coordinates");
		if (row >> extra) {
			if (extra != "shell" || !(row >> radius) || !std::isfinite(radius) || radius <= 0 || (row >> extra))
				throw std::invalid_argument("search_anchor optional suffix must be 'shell' and one positive finite radius");
		}
		const vec point(x, y, z);
		VINA_FOR(i, 3)
			if (!std::isfinite(point[i]) || point[i] < lower[i] + radius ||
			    point[i] > upper[i] - radius)
				throw std::invalid_argument("search_anchor point or shell must be finite and inside the search box");
		const auto it = heavy.find(std::stoul(serial_text));
		if (it == heavy.end())
			throw std::invalid_argument("search_anchor must name a real movable ligand heavy atom");
		const sz index = it->second;
		if (!assigned.insert(index).second)
			throw std::invalid_argument("search_anchor contains a repeated atom serial");
		sites.push_back({index, point, radius});
	}
	if (file.bad() || sites.empty()) throw std::invalid_argument("Cannot read any sites from search_anchor");
	return sites;
}

// Guide the grid search only. Workers share immutable sites and apply
// derivatives to their own model. Normalize so ring size cannot amplify the prior.
class anchored_search_grid : public igrid {
public:
	anchored_search_grid(const igrid& base, std::vector<search_anchor_site> sites)
		: base(base), sites(std::move(sites)) {
		if (this->sites.empty()) throw std::invalid_argument("search_anchor needs at least one site");
	}

	// Pseudo-Huber potential with a 1 Å length scale and 4 raw energy units.
	static fl penalty(const vec& offset, vec* derivative = nullptr, fl radius = 0) {
		const fl d2 = sqr(offset);
		if (radius > 0) {
			// Smooth the norm at the origin. The minimum remains at exactly radius.
			constexpr fl epsilon2 = 1e-8;
			const fl length = std::sqrt(d2 + epsilon2);
			const fl error = length - std::sqrt(radius * radius + epsilon2);
			const fl scale = std::sqrt(1 + error * error);
			if (derivative) *derivative = (4 * error / (scale * length)) * offset;
			return 4 * error * error / (scale + 1);
		}
		const fl scale = std::sqrt(1 + d2);
		if (derivative) *derivative = (4 / scale) * offset;
		return 4 * d2 / (scale + 1);
	}

	fl eval(const model& m, fl cap) const override {
		return base.eval(m, cap) + bias(m, nullptr);
	}
	fl eval_deriv(model& m, fl cap) const override {
		const fl native = base.eval_deriv(m, cap);
		return native + bias(m, &m.minus_forces);
	}
	fl eval_intra(model& m, fl cap) const override { return base.eval_intra(m, cap); }

private:
	fl bias(const model& m, vecv* gradients) const {
		fl energy = 0;
		const fl weight = 1.0 / sites.size();
		for (const auto& site : sites) {
			vec derivative;
			energy += penalty(m.coords[site.index] - site.target,
			                  gradients ? &derivative : nullptr, site.radius);
			if (gradients) (*gradients)[site.index] += weight * derivative;
		}
		return weight * energy;
	}
	const igrid& base;
	std::vector<search_anchor_site> sites;
};

#endif
