#ifndef VINA_POLAR_GEOMETRY_H
#define VINA_POLAR_GEOMETRY_H

#include "model.h"

// Extra final-ranking evidence from complementary N/O sites. Geometry and
// saturation are independent of the scoring function's reporting transform.
class polar_geometry {
public:
	struct terms {
		fl radial = 0;
		fl angular = 0;
	};

	explicit polar_geometry(const model& m)
		: ligand(collect(m.atoms, false)), receptor(collect(m.grid_atoms(), true)) {}

	struct burial_terms {
		fl burial = 0;
		fl geometry = 0;
		fl uncompensated = 0;
	};

	terms evaluate(const model& m) const {
		terms sum;
		for (const auto& lig : ligand) {
			const terms best = best_partner(m, lig);
			sum.radial += best.radial;
			sum.angular += best.angular;
		}
		return sum;
	}

	burial_terms evaluate_burial(const model& m) const {
		burial_terms sum;
		for (const auto& lig : ligand) {
			const fl buried = burial_fraction(m, m.coords[lig.index]);
			const terms best = best_partner(m, lig);
			sum.burial += buried;
			sum.geometry += best.radial - best.angular;
			sum.uncompensated += buried * (1 - best.angular);
		}
		return sum;
	}

private:
	friend class net_polar;
	struct site {
		sz index;
		bool donor;
		szv neighbors;
	};

	static fl burial_fraction(const model& m, const vec& xyz) {
		fl density = 0;
		for (const atom& rec : m.grid_atoms()) {
			if (rec.el == EL_TYPE_H || rec.el == EL_TYPE_Dummy || rec.el == EL_TYPE_SIZE)
				continue;
			const fl r2 = (rec.coords - xyz).norm_sqr();
			if (r2 >= 5.5 * 5.5)
				continue;
			density += r2 <= 3.5 * 3.5 ? 1 : (5.5 - std::sqrt(r2)) / 2;
			if (density >= 8)
				return 1;
		}
		return density / 8;
	}

	terms best_partner(const model& m, const site& lig) const {
		const vec& lig_xyz = m.coords[lig.index];
		terms best;
		for (const auto& rec : receptor) {
			const bool ligand_accepts = !lig.neighbors.empty() && rec.donor;
			const bool receptor_accepts = !rec.neighbors.empty() && lig.donor;
			if (!ligand_accepts && !receptor_accepts)
				continue;
			const vec& rec_xyz = m.grid_atoms()[rec.index].coords;
			const vec lig_to_rec = rec_xyz - lig_xyz;
			const fl distance_sqr = lig_to_rec.norm_sqr();
			if (distance_sqr <= 0 || distance_sqr > 16.0)
				continue;
			const fl distance = std::sqrt(distance_sqr);
			const fl deviation = (distance - 2.9) / 0.45;
			const fl weight = std::exp(-0.5 * deviation * deviation);
			const vec unit = (1 / distance) * lig_to_rec;
			if (ligand_accepts)
				update(best, weight, lig, lig_xyz, unit, [&](sz n) -> const vec& { return m.coords[n]; });
			if (receptor_accepts)
				update(best, weight, rec, rec_xyz, -1 * unit, [&](sz n) -> const vec& { return m.grid_atoms()[n].coords; });
		}
		return best;
	}

	static std::vector<site> collect(const atomv& atoms, bool in_grid) {
		std::vector<site> sites;
		VINA_FOR_IN(i, atoms) {
			const atom& a = atoms[i];
			if (a.el != EL_TYPE_N && a.el != EL_TYPE_O)
				continue;
			site polar{i, xs_is_donor(a.xs), {}};
			if (xs_is_acceptor(a.xs)) {
				for (const bond& edge : a.bonds) {
					// Receptor/ligand topology is built separately. Cross-body
					// indexes are never covalent neighbors of a polar acceptor.
					const atom_index& neighbor = edge.connected_atom_index;
					if (neighbor.in_grid != in_grid)
						continue;
					const sz element = atoms[neighbor.i].el;
					if (element != EL_TYPE_H && element != EL_TYPE_Dummy &&
					    element != EL_TYPE_Met && element != EL_TYPE_SIZE)
						polar.neighbors.push_back(neighbor.i);
				}
			}
			if (polar.donor || !polar.neighbors.empty() || (!in_grid && xs_is_acceptor(a.xs)))
				sites.push_back(std::move(polar));
		}
		return sites;
	}

	template<typename Coordinates>
	static fl angle_quality(const site& acceptor, const vec& xyz,
	                   const vec& towards_donor, const Coordinates& coordinate) {
		static const fl cos80 = std::cos(80 * pi / 180);
		static const fl cos110 = std::cos(110 * pi / 180);
		fl largest_cosine = -1;
		bool has_direction = false;
		for (sz index : acceptor.neighbors) {
			const vec towards_neighbor = coordinate(index) - xyz;
			const fl length = towards_neighbor.norm();
			if (length > 0) {
				has_direction = true;
				largest_cosine = std::max(largest_cosine, (towards_neighbor * towards_donor) / length);
			}
		}
		return has_direction ? std::max(fl(0), std::min(fl(1), (cos80 - largest_cosine) / (cos80 - cos110))) : -1;
	}

	template<typename Coordinates>
	static void update(terms& best, fl weight, const site& acceptor, const vec& xyz,
	                   const vec& towards_donor, const Coordinates& coordinate) {
		const fl gate = angle_quality(acceptor, xyz, towards_donor, coordinate);
		if (gate >= 0) {
			best.radial = std::max(best.radial, weight);
			best.angular = std::max(best.angular, weight * gate);
		}
	}

	std::vector<site> ligand;
	std::vector<site> receptor;
};

#endif
