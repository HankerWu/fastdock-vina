#ifndef VINA_POSE_TEMPLATE_H
#define VINA_POSE_TEMPLATE_H

#include <limits>
#include <map>
#include <stdexcept>

#include "common.h"

// The raw term is added before Vina's ligand-specific torsion transformation.
class pose_template {
public:
	using ligand_atoms = std::vector<std::pair<sz, std::string>>;
	static constexpr fl strength = 6.0;

	pose_template(std::istream& input, const ligand_atoms& ligand) {
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

	fl ligand_weight;
	fl reference_weight;
	std::vector<group> groups;
};

#endif
