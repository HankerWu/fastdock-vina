#ifndef VINA_NET_POLAR_H
#define VINA_NET_POLAR_H

#include "polar_geometry.h"

// Uses only ligand coordinates, fixed receptor coordinates, their Vina atom
// typing and covalent topology, and optional external water oxygen positions.
class net_polar {
public:
	static constexpr fl area_cost() { return 0.02; }

	net_polar(const model& m, const vecv& waters)
		: polar(m), ligand_site(index_sites(polar.ligand, m.atoms.size())),
		  receptor_site(index_sites(polar.receptor, m.grid_atoms().size())),
		  ligand_radii(radii(m.atoms)), receptor_radii(radii(m.grid_atoms())),
		  bridges(prepare_waters(m, waters)) {}

	fl direction(const model& m, sz lig, sz rec, const vec& lig_xyz, const vec& rec_xyz) const {
		const vec displacement = rec_xyz - lig_xyz;
		const fl distance = displacement.norm();
		if (distance == 0)
			return 1;
		const vec unit = (1 / distance) * displacement;
		const auto& a = m.atoms[lig];
		const auto& b = m.grid_atoms()[rec];
		fl gate = 0;
		if (xs_donor_acceptor(b.xs, a.xs))
			gate = direct_gate(ligand_site[lig], lig_xyz, unit,
			                   [&](sz i) -> const vec& { return m.coords[i]; });
		if (xs_donor_acceptor(a.xs, b.xs))
			gate = std::max(gate, direct_gate(receptor_site[rec], rec_xyz, -1 * unit,
			                                 [&](sz i) -> const vec& { return m.grid_atoms()[i].coords; }));
		return gate;
	}

	void add_surface(const model& m, net_polar_terms& terms, bool satisfaction) const {
		std::vector<sphere> ligand_neighbors, receptor_neighbors;
		for (const auto& site : polar.ligand) {
			const sz index = site.index;
			const vec& center = m.coords[index];
			const fl radius = ligand_radii[index];
			ligand_neighbors.clear();
			receptor_neighbors.clear();
			VINA_FOR_IN(i, ligand_radii) {
				if (i != index)
					add_neighbor(ligand_neighbors, center, radius, m.coords[i], ligand_radii[i]);
			}
			VINA_FOR_IN(i, receptor_radii)
				add_neighbor(receptor_neighbors, center, radius, m.grid_atoms()[i].coords, receptor_radii[i]);

			if (receptor_neighbors.empty())
				continue;
			sz lost = 0;
			for (const auto& unit : directions()) {
				const vec point = radius * unit;
				if (!occluded(point, ligand_neighbors) && occluded(point, receptor_neighbors))
					++lost;
			}
			const fl area = 4 * pi * sqr(radius) * lost / directions().size();
			terms.buried_area += area;
			if (area == 0)
				continue;
			const fl water = bridges.empty() ? 0 : best_bridge(m, site);
			terms.water_area += area * water;
			if (satisfaction) {
				const fl direct = polar.best_partner(m, site).angular;
				terms.direct_credit_area += area * direct;
				// Charge each site for the uncredited fraction of its lost area.
				// A bridge can only improve upon that site's direct credit.
				terms.additional_water_credit_area += area * std::max(fl(0), water - direct);
			}
		}
	}

private:
	using site = polar_geometry::site;
	struct sphere {
		vec relative_center;
		fl radius_squared;
	};
	struct water_leg {
		vec direction;
		fl quality;
	};
	struct water_bridge {
		vec oxygen;
		std::vector<water_leg> legs;
	};

	static std::vector<const site*> index_sites(const std::vector<site>& sites, sz atom_count) {
		std::vector<const site*> indices(atom_count, nullptr);
		for (const auto& entry : sites)
			indices[entry.index] = &entry;
		return indices;
	}

	template<class Coordinates>
	static fl direct_gate(const site* acceptor, const vec& center, const vec& unit,
	                      const Coordinates& coordinates) {
		if (!acceptor)
			return 1;
		const fl gate = polar_geometry::angle_quality(*acceptor, center, unit, coordinates);
		// A missing covalent direction carries no new geometric evidence.
		return gate < 0 ? 1 : gate;
	}

	static flv radii(const atomv& atoms) {
		// Heavy-atom van der Waals radii by EL_TYPE; H and dummy sites
		// carry no separate sphere in the united-atom surface model.
		constexpr fl vdw[] = {0, 1.70, 1.55, 1.52, 1.80, 1.80, 1.47,
		                      1.75, 1.85, 1.98, 2.10, 2.02, 1.80, 0};
		flv result;
		result.reserve(atoms.size());
		for (const auto& a : atoms) {
			const fl radius = a.el < EL_TYPE_SIZE ? vdw[a.el] : 0;
			result.push_back(radius > 0 ? radius + 1.4 : 0);
		}
		return result;
	}

	static const std::array<vec, 256>& directions() {
		static const auto points = [] {
			std::array<vec, 256> samples;
			const fl golden_angle = pi * (3 - std::sqrt(5.0));
			VINA_FOR_IN(i, samples) {
				const fl z = 1 - 2 * (i + 0.5) / samples.size();
				const fl rho = std::sqrt(1 - sqr(z));
				const fl phi = i * golden_angle;
				samples[i] = vec(rho * std::cos(phi), rho * std::sin(phi), z);
			}
			return samples;
		}();
		return points;
	}

	static void add_neighbor(std::vector<sphere>& neighbors, const vec& center, fl radius,
	                         const vec& other, fl other_radius) {
		if (other_radius == 0)
			return;
		const vec relative = other - center;
		if (relative.norm_sqr() < sqr(radius + other_radius))
			neighbors.push_back({relative, sqr(other_radius)});
	}

	static bool occluded(const vec& point, const std::vector<sphere>& spheres) {
		return std::any_of(spheres.begin(), spheres.end(), [&](const sphere& s) {
			return (point - s.relative_center).norm_sqr() < s.radius_squared;
		});
	}

	template<class Coordinates>
	static fl contact(const site& polar_site, const vec& center, const vec& to_water,
	                  fl distance, const Coordinates& coordinates) {
		if (distance < 2.4 || distance > 3.5)
			return 0;
		const fl radial = std::exp(-0.5 * sqr((distance - 2.9) / 0.35));
		const fl gate = polar_site.donor ? 1 : std::max(fl(0),
			polar_geometry::angle_quality(polar_site, center, (1 / distance) * to_water, coordinates));
		return radial * gate;
	}

	std::vector<water_bridge> prepare_waters(const model& m, const vecv& waters) const {
		std::vector<water_bridge> prepared;
		for (const auto& oxygen : waters) {
			bool clashes = false;
			VINA_FOR_IN(i, receptor_radii) {
				if (receptor_radii[i] > 0 && (oxygen - m.grid_atoms()[i].coords).norm_sqr() < sqr(2.2)) {
					clashes = true;
					break;
				}
			}
			if (clashes)
				continue;
			water_bridge bridge{oxygen, {}};
			for (const auto& protein : polar.receptor) {
				const vec& center = m.grid_atoms()[protein.index].coords;
				const vec to_water = oxygen - center;
				const fl distance = to_water.norm();
				const fl quality = contact(protein, center, to_water, distance,
					[&](sz i) -> const vec& { return m.grid_atoms()[i].coords; });
				if (quality > 0)
					bridge.legs.push_back({(-1 / distance) * to_water, quality});
			}
			if (!bridge.legs.empty())
				prepared.push_back(std::move(bridge));
		}
		return prepared;
	}

	fl best_bridge(const model& m, const site& ligand) const {
		static const fl cos75 = std::cos(75 * pi / 180);
		static const fl cos135 = std::cos(135 * pi / 180);
		fl best = 0;
		const vec& center = m.coords[ligand.index];
		for (const auto& water : bridges) {
			const vec to_water = water.oxygen - center;
			const fl distance = to_water.norm();
			const fl quality = contact(ligand, center, to_water, distance,
			                           [&](sz i) -> const vec& { return m.coords[i]; });
			if (quality <= best)
				continue;
			bool clashes = false;
			VINA_FOR_IN(i, ligand_radii) {
				if (ligand_radii[i] > 0 && (water.oxygen - m.coords[i]).norm_sqr() < sqr(2.2)) {
					clashes = true;
					break;
				}
			}
			if (clashes)
				continue;
			for (const auto& leg : water.legs) {
				const fl cosine = ((-1 / distance) * to_water) * leg.direction;
				if (cosine >= cos135 && cosine <= cos75)
					best = std::max(best, quality * leg.quality);
			}
		}
		return best;
	}

	const polar_geometry polar;
	const std::vector<const site*> ligand_site;
	const std::vector<const site*> receptor_site;
	const flv ligand_radii;
	const flv receptor_radii;
	const std::vector<water_bridge> bridges;
};

#endif
