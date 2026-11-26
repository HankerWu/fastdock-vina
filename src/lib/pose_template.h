#ifndef VINA_POSE_TEMPLATE_H
#define VINA_POSE_TEMPLATE_H

#include <array>
#include <limits>
#include <map>
#include <boost/optional.hpp>
#include <stdexcept>

#include "atom.h"

// Two receptor-indexed contact channels. Each receptor atom contributes the
// strongest distance contact with any heavy ligand atom, or any N/O for the
// polar channel. A receptor atom can contribute to either molecule alone.
class receptor_contact_template {
public:
	struct point { vec position; bool polar; };
	using points = std::vector<point>;
	using ligand_atoms = std::vector<std::pair<sz, std::string>>;

	receptor_contact_template(const atomv& receptor, const points& reference, const ligand_atoms& ligand) {
		for (const atom& a : receptor) {
			if (a.el != EL_TYPE_H && a.el != EL_TYPE_Dummy)
				sites.push_back({a.coords, a.el == EL_TYPE_N || a.el == EL_TYPE_O});
		}
		for (const auto& a : ligand)
			ligand_indices.emplace_back(a.first, is_polar(a.second));
		visit_contacts(reference, [&](sz i, const channels& w) {
			sites[i].reference = w;
			accumulate_product(reference_norm, w, w);
		});
	}

	static bool is_polar(const std::string& element) { return element == "N" || element == "O"; }

	fl overlap(const vecv& coords) const {
		points ligand;
		ligand.reserve(ligand_indices.size());
		for (const auto& a : ligand_indices)
			ligand.push_back({coords[a.first], a.second});
		channels dot{}, ligand_norm{};
		visit_contacts(ligand, [&](sz i, const channels& w) {
			accumulate_product(dot, w, sites[i].reference);
			accumulate_product(ligand_norm, w, w);
		});
		fl sum = 0;
		sz active = 0;
		for (sz k = 0; k < dot.size(); ++k) {
			const fl denominator = ligand_norm[k] + reference_norm[k] - dot[k];
			if (denominator > 0) {
				sum += dot[k] / denominator;
				++active;
			}
		}
		return active ? sum / active : 0;
	}

private:
	using channels = std::array<fl, 2>;
	struct site { vec position; bool polar; channels reference{}; };
	std::vector<site> sites;
	std::vector<std::pair<sz, bool>> ligand_indices;
	channels reference_norm{};

	static void accumulate_product(channels& sum, const channels& a, const channels& b) {
		for (sz k = 0; k < sum.size(); ++k)
			sum[k] += a[k] * b[k];
	}

	static fl contact_weight(fl nearest_sqr, fl inner, fl outer) {
		if (nearest_sqr <= sqr(inner)) return 1;
		if (nearest_sqr >= sqr(outer)) return 0;
		return (outer - std::sqrt(nearest_sqr)) / (outer - inner);
	}

	template<typename Visitor>
	void visit_contacts(const points& ligand, Visitor visit) const {
		// The constructor already requires a heavy atom on each molecule.
		vec lower = ligand.front().position, upper = lower;
		for (const point& a : ligand) {
			for (sz k = 0; k < 3; ++k) {
				lower[k] = std::min(lower[k], a.position[k]);
				upper[k] = std::max(upper[k], a.position[k]);
			}
		}
		lower -= 6;
		upper += 6;
		// Crop separately around each pose. The cached reference norm still
		// includes all reference-only sites, including those outside this box.
		VINA_FOR_IN(i, sites) {
			const site& s = sites[i];
			const vec& r = s.position;
			if (r[0] < lower[0] || r[0] > upper[0] || r[1] < lower[1] || r[1] > upper[1] ||
			    r[2] < lower[2] || r[2] > upper[2]) continue;
			fl general_min = sqr(6.0), polar_min = sqr(4.5);
			for (const point& a : ligand) {
				const fl d = vec_distance_sqr(r, a.position);
				general_min = std::min(general_min, d);
				if (s.polar && a.polar) polar_min = std::min(polar_min, d);
			}
			visit(i, channels{{contact_weight(general_min, 4, 6), contact_weight(polar_min, 3, 4.5)}});
		}
	}
};

// The raw term is added before Vina's ligand-specific torsion transformation.
class pose_template {
public:
	using ligand_atoms = std::vector<std::pair<sz, std::string>>;
	static constexpr fl strength = 6.0;

	pose_template(std::istream& input, const ligand_atoms& ligand, const atomv* contact_receptor = nullptr) {
		std::map<std::string, group> by_element;
		sz reference_count = 0;
		sz number = 0;
		std::string line;
		while (std::getline(input, line)) {
			++number;
			if (line.find_first_not_of(" \t\r") == std::string::npos)
				continue;
			std::istringstream row(line);
			std::string element, extra;
			fl x, y, z;
			if (!(row >> element >> x >> y >> z) || row >> extra ||
			    !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
				throw std::invalid_argument("Expected element and three finite coordinates at pose template line " + std::to_string(number));
			if (element == "H" || element == "D" || element == "T" ||
			    element == "G0" || element == "G1" || element == "G2" || element == "G3")
				continue;
			if (!element_symbol(element))
				throw std::invalid_argument("Invalid element at pose template line " + std::to_string(number) + ": " + element);
			by_element[element].reference.emplace_back(x, y, z);
			++reference_count;
		}
		if (input.bad())
			throw std::invalid_argument("Could not read pose template");
		if (reference_count == 0 || ligand.empty())
			throw std::invalid_argument("Pose template and ligand must each contain real heavy atoms");
		if (contact_receptor) {
			receptor_contact_template::points reference;
			for (const auto& entry : by_element) {
				for (const vec& position : entry.second.reference)
					reference.push_back({position, receptor_contact_template::is_polar(entry.first)});
			}
			contacts.emplace(*contact_receptor, reference, ligand);
			return;
		}
		ligand_weight = 0.5 / ligand.size();
		reference_weight = 0.5 / reference_count;
		for (const auto& atom : ligand) {
			auto element = by_element.find(atom.second);
			if (element != by_element.end())
				element->second.ligand.push_back(atom.first);
		}
		for (auto& entry : by_element) {
			group& atoms = entry.second;
			if (!atoms.ligand.empty()) {
				atoms.reverse_distance.resize(atoms.reference.size());
				atoms.reverse_nearest.resize(atoms.reference.size());
				groups.push_back(std::move(atoms));
			}
		}
	}

	// Returns overlap. If requested, add the derivative of -strength * overlap
	// to the Cartesian gradients before the model projects them onto pose DOFs.
	fl overlap(const vecv& coords, vecv* gradient = nullptr) {
		if (contacts) {
			VINA_CHECK(!gradient); // contact_rank only scores the native retained poses
			return contacts->overlap(coords);
		}
		constexpr fl sigma_sqr = 1.5 * 1.5;
		const fl infinity = std::numeric_limits<fl>::infinity();
		fl ligand_sum = 0;
		fl reference_sum = 0;
		for (group& atoms : groups) {
			std::fill(atoms.reverse_distance.begin(), atoms.reverse_distance.end(), infinity);
			for (sz ligand_index : atoms.ligand) {
				const vec& position = coords[ligand_index];
				fl forward_distance = infinity;
				sz forward_nearest = 0;
				VINA_FOR_IN(j, atoms.reference) {
					const fl distance = vec_distance_sqr(position, atoms.reference[j]);
					if (distance < forward_distance) {
						forward_distance = distance;
						forward_nearest = j;
					}
					if (distance < atoms.reverse_distance[j]) {
						atoms.reverse_distance[j] = distance;
						atoms.reverse_nearest[j] = ligand_index;
					}
				}
				const fl similarity = std::exp(-forward_distance / (2 * sigma_sqr));
				ligand_sum += similarity;
				if (gradient && similarity > 0)
					(*gradient)[ligand_index] += (strength * ligand_weight * similarity / sigma_sqr)
					                             * (position - atoms.reference[forward_nearest]);
			}
			VINA_FOR_IN(j, atoms.reference) {
				const fl similarity = std::exp(-atoms.reverse_distance[j] / (2 * sigma_sqr));
				reference_sum += similarity;
				if (gradient && similarity > 0) {
					const sz ligand_index = atoms.reverse_nearest[j];
					(*gradient)[ligand_index] += (strength * reference_weight * similarity / sigma_sqr)
					                             * (coords[ligand_index] - atoms.reference[j]);
				}
			}
		}
		return ligand_sum * ligand_weight + reference_sum * reference_weight;
	}

private:
	struct group {
		szv ligand;
		vecv reference;
		flv reverse_distance;
		szv reverse_nearest;
	};

	static bool element_symbol(const std::string& name) {
		static const char* const symbols[] = {
			"He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
			"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
			"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba",
			"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
			"W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
			"U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
			"Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
		};
		return std::find(std::begin(symbols), std::end(symbols), name) != std::end(symbols);
	}

	boost::optional<receptor_contact_template> contacts;
	fl ligand_weight = 0;
	fl reference_weight = 0;
	std::vector<group> groups;
};

#endif
