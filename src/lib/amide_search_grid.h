#ifndef VINA_AMIDE_SEARCH_GRID_H
#define VINA_AMIDE_SEARCH_GRID_H

#include <array>
#include <set>
#include <sstream>

#include "igrid.h"
#include "model.h"
#include "parse_pdbqt.h"

struct amide_search_site {
	sz oxygen, carbon, nitrogen;
	std::array<sz, 2> arms;
	sz arm_count;
	fl barrier;

	// Project both terminal bonds perpendicular to the C--N unit axis.
	// For a single arm, E = B(1 - (u.v)^2) with u and v unit projections.
	fl bias(const vecv& coords, vecv* gradients = nullptr) const {
		const vec axis = coords[nitrogen] - coords[carbon];
		const fl length = axis.norm();
		if (!(length > 1e-8)) throw std::invalid_argument("search_amide has a degenerate C-N axis");
		const vec w = (1 / length) * axis;
		const vec oc = coords[oxygen] - coords[carbon];
		const fl oc_axial = oc * w;
		const vec op = oc - oc_axial * w;
		const fl op_len = op.norm();
		if (!(op_len > 1e-8)) throw std::invalid_argument("search_amide has a collinear O-C-N angle");
		const vec u = (1 / op_len) * op;
		const fl weight = barrier / arm_count;
		fl energy = 0;
		for (sz i = 0; i < arm_count; ++i) {
			const sz x = arms[i];
			const vec xn = coords[x] - coords[nitrogen];
			const fl xn_axial = xn * w;
			const vec xp = xn - xn_axial * w;
			const fl xp_len = xp.norm();
			if (!(xp_len > 1e-8)) throw std::invalid_argument("search_amide has a collinear C-N-X angle");
			const vec v = (1 / xp_len) * xp;
			const fl sine = w * cross_product(u, v);
			energy += weight * sqr(sine);
			if (gradients) {
				const fl cosine = u * v;
				const fl scale = -2 * weight * cosine;
				const vec d_oc = (scale / op_len) * (v - cosine * u);
				const vec d_xn = (scale / xp_len) * (u - cosine * v);
				const vec d_axis = (-1 / length) * (oc_axial * d_oc + xn_axial * d_xn);
				(*gradients)[oxygen] += d_oc;
				(*gradients)[carbon] -= d_oc + d_axis;
				(*gradients)[nitrogen] += d_axis - d_xn;
				(*gradients)[x] += d_xn;
			}
		}
		return energy;
	}
};

// PDBQT supplies elements, bonds and active rotors. The declared bond orders
// and formal amide chemistry come from the prepared molecule's chemical graph.
inline std::vector<amide_search_site> read_search_amides(std::istream& file, const model& m) {
	const auto indices = ligand_heavy_indices_by_serial(m);
	vecv initial;
	initial.reserve(m.num_movable_atoms());
	for (sz i = 0; i < m.num_movable_atoms(); ++i) initial.push_back(m.movable_coords(i));
	std::vector<amide_search_site> sites;
	std::set<std::pair<sz, sz>> bonds;
	std::string line;
	while (std::getline(file, line)) {
		std::istringstream row(line);
		std::vector<std::string> fields;
		for (std::string token; row >> token;) fields.push_back(token);
		if (fields.empty()) continue;
		if (fields.size() != 5 && fields.size() != 6)
			throw std::invalid_argument("search_amide rows require O C N X [X2] B");
		fl barrier;
		std::istringstream value(fields.back());
		char extra;
		if (!(value >> barrier) || (value >> extra) || !std::isfinite(barrier) || barrier < 0)
			throw std::invalid_argument("search_amide requires a finite nonnegative barrier");
		fields.pop_back();
		szv ids;
		std::set<sz> distinct;
		for (const auto& field : fields) {
			if (field.size() > 5 || field.find_first_not_of("0123456789") != std::string::npos)
				throw std::invalid_argument("search_amide requires ligand atom serials");
			const auto it = indices.find(std::stoul(field));
			if (it == indices.end() || !distinct.insert(it->second).second)
				throw std::invalid_argument("search_amide requires distinct movable ligand heavy atoms");
			ids.push_back(it->second);
		}
		const sz o = ids[0], c = ids[1], n = ids[2];
		if (m.movable_atom(o).el != EL_TYPE_O || m.movable_atom(c).el != EL_TYPE_C ||
		    m.movable_atom(n).el != EL_TYPE_N || !bonds.emplace(c, n).second)
			throw std::invalid_argument("search_amide requires a unique O-C-N site");
		const auto oxygen = m.get_atom(o);
		if (oxygen.bonds.size() != 1 || !(oxygen.bonds[0].connected_atom_index == atom_index(c, false)))
			throw std::invalid_argument("search_amide oxygen must be terminal and bonded to the carbonyl carbon");
		bool rotor = false;
		for (const auto& bond : m.get_atom(c).bonds) {
			const auto other = bond.connected_atom_index;
			if (other.in_grid) throw std::invalid_argument("search_amide requires a ligand site");
			if (other.i == n) rotor = bond.rotatable;
			else if (other.i != o && m.movable_atom(other.i).el != EL_TYPE_C && m.movable_atom(other.i).el != EL_TYPE_H)
				throw std::invalid_argument("search_amide requires a carboxamide carbonyl");
		}
		if (!rotor) throw std::invalid_argument("search_amide C-N must be a bonded active rotor");
		const std::set<sz> arms(ids.begin() + 3, ids.end());
		std::set<sz> neighbors;
		for (const auto& bond : m.get_atom(n).bonds) {
			const auto other = bond.connected_atom_index;
			if (other.in_grid) throw std::invalid_argument("search_amide requires a ligand site");
			if (other.i == c || m.movable_atom(other.i).el == EL_TYPE_H) continue;
			if (m.movable_atom(other.i).el != EL_TYPE_C)
				throw std::invalid_argument("search_amide nitrogen substituents must be carbon atoms");
			neighbors.insert(other.i);
		}
		if (arms != neighbors)
			throw std::invalid_argument("search_amide requires all carbon arms bonded to the amide nitrogen");
		const sz second = ids.size() == 5 ? ids[4] : ids[3];
		amide_search_site site{o, c, n, {{ids[3], second}}, ids.size() - 3, barrier};
		site.bias(initial);  // Catch undefined dihedrals before launching worker threads.
		sites.push_back(site);
	}
	if (file.bad() || sites.empty()) throw std::invalid_argument("Cannot read any sites from search_amide");
	return sites;
}

class amide_search_grid : public igrid {
public:
	amide_search_grid(const igrid& base, std::vector<amide_search_site> sites)
		: base(base), sites(std::move(sites)) {}
	static fl penalty(const model& m, const std::vector<amide_search_site>& sites, vecv* gradients = nullptr) {
		fl total = 0;
		for (const auto& site : sites) total += site.bias(m.coords, gradients);
		return total;
	}
	fl eval(const model& m, fl cap) const override {
		return base.eval(m, cap) + penalty(m, sites);
	}
	fl eval_deriv(model& m, fl cap) const override {
		const fl native = base.eval_deriv(m, cap);
		return native + penalty(m, sites, &m.minus_forces);
	}
	fl eval_intra(model& m, fl cap) const override { return base.eval_intra(m, cap); }
private:
	const igrid& base;
	std::vector<amide_search_site> sites;
};

#endif
