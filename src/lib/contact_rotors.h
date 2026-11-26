#ifndef VINA_CONTACT_ROTORS_H
#define VINA_CONTACT_ROTORS_H

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "model.h"

// An experimental terminal score. It never participates in the search potential.
class contact_rotors {
public:
	explicit contact_rotors(const model& m) {
		if (m.num_ligands() != 1 || m.num_flex() || m.grid_atoms().empty())
			throw std::invalid_argument("contact_rotor_score requires one ligand and an explicit rigid receptor");
		const sz n = m.num_movable_atoms();
		for (sz i = 0; i < n; ++i) {
			const auto& a = m.get_atom(atom_index(i, false));
			if (a.el == EL_TYPE_Dummy || (a.ad >= AD_TYPE_CG0 && a.ad <= AD_TYPE_CG3))
				throw std::invalid_argument("contact_rotor_score requires a ligand without hydration or macrocycle pseudo-atoms");
			if (!a.is_hydrogen())
				heavy.push_back(i);
		}
		for (sz i : heavy) {
			for (const auto& bond : m.get_atom(atom_index(i, false)).bonds) {
				const auto j = bond.connected_atom_index;
				if (!bond.rotatable || j.in_grid || j.i <= i || j.i >= n || m.get_atom(j).is_hydrogen())
					continue;
				// Cut the bond in the chemical graph. A different covalent path
				// between its endpoints means this is not an independent rotor.
				std::vector<bool> side(n, false);
				szv visit(1, j.i);
				side[j.i] = true;
				for (sz k = 0; k < visit.size(); ++k) {
					const sz current = visit[k];
					for (const auto& edge : m.get_atom(atom_index(current, false)).bonds) {
						const auto dst = edge.connected_atom_index;
						if (dst.in_grid || dst.i >= n || side[dst.i] || (current == j.i && dst.i == i))
							continue;
						side[dst.i] = true;
						visit.push_back(dst.i);
					}
				}
				if (side[i])
					continue;
				const sz side_size = std::count_if(heavy.begin(), heavy.end(), [&](sz a) { return side[a]; });
				if (side_size > 1 && heavy.size() - side_size > 1)
					cuts.push_back({i, j.i, std::move(side)});
			}
		}
	}

	std::string remark(const model& m, fl raw_score) const {
		szv contacts;
		if (!cuts.empty()) {
			for (sz i : heavy) {
				const vec& xyz = m.movable_coords(i);
				for (const atom& rec : m.grid_atoms()) {
					if (!rec.is_hydrogen() && rec.el != EL_TYPE_Dummy && vec_distance_sqr(xyz, rec.coords) <= 16) {
						contacts.push_back(i);
						break;
					}
				}
			}
		}
		sz count = 0;
		for (const auto& cut : cuts) {
			const vec origin = m.movable_coords(cut.a);
			const vec axis = m.movable_coords(cut.b) - origin;
			const fl norm2 = sqr(axis);
			if (norm2 < epsilon_fl)
				throw std::invalid_argument("contact_rotor_score found a zero-length bond");
			bool sides[2] = {false, false};
			for (sz i : contacts) {
				if (sqr(cross_product(m.movable_coords(i) - origin, axis)) >= norm2)
					sides[cut.side[i]] = true;
				if (sides[0] && sides[1]) {
					++count;
					break;
				}
			}
		}
		static const fl cost = 0.00198720425864 * 298.15 * std::log(3.0);
		std::ostringstream out;
		out << "REMARK VINA CONTACT ROTOR SCORE: " << std::fixed << std::setprecision(6)
		    << raw_score + cost * count << ' ' << count << ' ' << cuts.size() << '\n';
		return out.str();
	}

private:
	struct cut { sz a, b; std::vector<bool> side; };
	szv heavy;
	std::vector<cut> cuts;
};

#endif
