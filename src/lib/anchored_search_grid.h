#ifndef VINA_ANCHORED_SEARCH_GRID_H
#define VINA_ANCHORED_SEARCH_GRID_H

#include <map>
#include <set>
#include <tuple>

#include "igrid.h"
#include "parse_pdbqt.h"
#include "model.h"

struct search_anchor_site {
	sz index;
	vec target;
	fl radius;
};

using anchor_assignment = std::vector<search_anchor_site>;

inline std::vector<anchor_assignment> read_search_anchors(std::istream& file, const model& ligand,
                                                       const vec& lower, const vec& upper) {
	std::vector<anchor_assignment> alternatives(1);
	std::set<sz> assigned;
	const auto heavy = ligand_heavy_indices_by_serial(ligand);
	std::string line;
	while (std::getline(file, line)) {
		if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
		std::istringstream row(line);
		std::string serial_text;
		fl x, y, z, radius = 0;
		std::string extra;
		row >> serial_text;
		if (serial_text == "OR") {
			if ((row >> extra) || alternatives.back().empty())
				throw std::invalid_argument("search_anchor requires a nonempty assignment before OR");
			alternatives.emplace_back();
			assigned.clear();
			continue;
		}
		if (!(row >> x >> y >> z) || serial_text.size() > 5 ||
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
		alternatives.back().push_back({index, point, radius});
	}
	if (file.bad() || alternatives.back().empty()) throw std::invalid_argument("Cannot read any sites from search_anchor");
	return alternatives;
}

// Guide the grid search only. Workers share immutable sites and apply
// derivatives to their own model. Normalize so ring size cannot amplify the prior.
class anchored_search_grid : public igrid {
public:
	anchored_search_grid(const igrid& base, anchor_assignment sites)
		: anchored_search_grid(base, std::vector<anchor_assignment>{std::move(sites)}) {}

	anchored_search_grid(const igrid& base, std::vector<anchor_assignment> alternatives)
		: base(base) {
		if (alternatives.empty()) throw std::invalid_argument("search_anchor needs at least one assignment");
		// Shared terms are evaluated once even if a ring has many valid orientations.
		std::map<std::tuple<sz, fl, fl, fl, fl>, sz> indices;
		std::set<szv> seen;
		for (const auto& alternative : alternatives) {
			if (alternative.empty()) throw std::invalid_argument("search_anchor needs nonempty assignments");
			szv ids;
			for (const auto& site : alternative) {
				const auto key = std::make_tuple(site.index, site.target[0], site.target[1], site.target[2], site.radius);
				const auto entry = indices.emplace(key, sites.size());
				if (entry.second) sites.push_back(site);
				ids.push_back(entry.first->second);
			}
			std::sort(ids.begin(), ids.end());
			if (seen.insert(ids).second) assignments.push_back(std::move(ids));
		}
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
	fl single_bias(const model& m, vecv* gradients) const {
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
	fl bias(const model& m, vecv* gradients) const {
		if (assignments.size() == 1) return single_bias(m, gradients);
		// Immutable grids are shared by all workers. Scratch is private to each calling thread.
		struct workspace { flv energies, weights, coefficients; vecv derivatives; };
		thread_local workspace w;
		w.energies.resize(sites.size());
		if (gradients) w.derivatives.resize(sites.size());
		VINA_FOR_IN(i, sites) {
			const auto& site = sites[i];
			w.energies[i] = penalty(m.coords[site.index] - site.target,
			                        gradients ? &w.derivatives[i] : nullptr, site.radius);
		}
		w.weights.resize(assignments.size());
		fl lowest = max_fl;
		VINA_FOR_IN(j, assignments) {
			fl total = 0;
			for (sz i : assignments[j]) total += w.energies[i];
			w.weights[j] = total / assignments[j].size();
			lowest = std::min(lowest, w.weights[j]);
		}
		// A fixed smoothing scale in raw search units, not an assay temperature.
		constexpr fl tau = 0.6;
		fl sum = 0;
		for (auto& value : w.weights) {
			value = std::exp((lowest - value) / tau);
			sum += value;
		}
		if (gradients) {
			w.coefficients.assign(sites.size(), 0);
			VINA_FOR_IN(j, assignments) {
				const fl fraction = w.weights[j] / (sum * assignments[j].size());
				for (sz i : assignments[j]) w.coefficients[i] += fraction;
			}
			VINA_FOR_IN(i, sites)
				(*gradients)[sites[i].index] += w.coefficients[i] * w.derivatives[i];
		}
		return lowest - tau * std::log(sum / assignments.size());
	}

	const igrid& base;
	anchor_assignment sites;
	std::vector<szv> assignments;
};

#endif
