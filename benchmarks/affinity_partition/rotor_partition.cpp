// Fixed-pose thermodynamic diagnostic. All coordinate-independent input topology
// and native inter/intramolecular potential tables are shared across angular samples.
#include "non_cache.h"
#include "parse_pdbqt.h"
#include "precalculate.h"

#include <map>

struct rotation {
	sz first, second;
	std::vector<szv> sides;
};

class rotor_partition {
  public:
	explicit rotor_partition(model &system) : m(system) {
		// BRANCH records describe the actual kinematic tree. Distance-inferred
		// covalent neighbours may also contain spurious cross-fragment contacts.
		std::map<sz, sz> serial_to_index;
		std::vector<rotation> branches;
		szv open;
		for (const auto &entry : m.ligands[0].cont) {
			const auto &line = entry.first;
			if (entry.second) {
				const sz index = *entry.second;
				const sz serial = std::stoul(line.substr(6, 5));
				if (!serial_to_index.emplace(serial, index).second)
					throw std::runtime_error("Duplicate ligand serial");
				for (sz branch : open)
					branches[branch].sides.front().push_back(index);
			} else if (line.compare(0, 7, "BRANCH ") == 0) {
				std::istringstream fields(line.substr(7));
				sz first, second;
				if (!(fields >> first >> second))
					throw std::runtime_error("Invalid BRANCH");
				open.push_back(branches.size());
				branches.push_back(rotation{first, second, {szv{}}});
			} else if (line.compare(0, 9, "ENDBRANCH") == 0) {
				if (open.empty())
					throw std::runtime_error("Unbalanced ENDBRANCH");
				open.pop_back();
			}
		}
		if (!open.empty() || serial_to_index.size() != m.num_movable_atoms())
			throw std::runtime_error("Incomplete ligand torsion tree");

		for (const auto &branch : branches) {
			szv right = branch.sides.front(), left;
			std::sort(right.begin(), right.end());
			for (sz i = 0; i < m.num_movable_atoms(); ++i)
				if (!std::binary_search(right.begin(), right.end(), i))
					left.push_back(i);
			const auto heavy = [&](const szv &side) {
				return std::count_if(side.begin(), side.end(),
									 [&](sz i) { return !m.atoms[i].is_hydrogen(); });
			};
			const auto l = heavy(left), r = heavy(right);
			if (std::min(l, r) < 2)
				continue;
			rotation rotor{serial_to_index.at(branch.first), serial_to_index.at(branch.second), {}};
			if (l <= r)
				rotor.sides.push_back(left);
			if (r <= l)
				rotor.sides.push_back(right);
			rotors.push_back(rotor);
		}
	}

	static fl free_energy(const flv &energies, fl rt) {
		const fl best = *std::min_element(energies.begin(), energies.end());
		fl weight = 0;
		for (fl e : energies)
			weight += std::exp(-(e - best) / rt);
		return best - rt * std::log(weight / energies.size());
	}

	void evaluate(const non_cache &grid, const precalculate_byatom &pairs, sz steps, fl rt, sz mode) {
		const vecv origin = m.coords;
		const vec cap(1000, 1000, 1000);
		const fl e0 = grid.eval(m, 1000);
		const fl i0 = m.evali(pairs, cap);
		fl correction = 0, coarse = 0;
		for (const auto &rotor : rotors) {
			const vec pivot = origin[rotor.first];
			const vec arm = origin[rotor.second] - pivot;
			if (arm.norm() < 1e-6)
				throw std::runtime_error("Degenerate rotor axis");
			const vec axis = (1.0 / arm.norm()) * arm;
			flv full, free, full_half, free_half;
			for (const auto &side : rotor.sides) {
				for (sz step = 0; step < steps; ++step) {
					const fl angle = 2 * pi * step / steps;
					for (sz atom : side) {
						const vec v = origin[atom] - pivot;
						m.coords[atom] = pivot + std::cos(angle) * v +
										 std::sin(angle) * cross_product(axis, v) +
										 ((1 - std::cos(angle)) * (axis * v)) * axis;
					}
					const fl intra = m.evali(pairs, cap);
					const fl total = intra + grid.eval(m, 1000);
					free.push_back(intra);
					full.push_back(total);
					if (!(step % 2)) {
						free_half.push_back(intra);
						full_half.push_back(total);
					}
				}
				m.coords = origin;
			}
			correction += free_energy(full, rt) - free_energy(free, rt) - e0;
			coarse += free_energy(full_half, rt) - free_energy(free_half, rt) - e0;
		}
		std::cout << mode << ',' << rotors.size() << ',' << e0 << ',' << i0 << ',' << correction << ','
				  << coarse << '\n';
	}

	static grid_dims scan_bounds(const model &m, const std::vector<std::string> &poses) {
		grid_dims bounds;
		for (auto &dim : bounds) {
			dim.begin = max_fl;
			dim.end = -max_fl;
			dim.n_voxels = 1;
		}
		model probe = m;
		for (const auto &pose : poses) {
			replace_fixed_pose_coordinates(probe, pose);
			fl diameter = 0;
			for (const vec &a : probe.coords)
				for (const vec &b : probe.coords)
					diameter = std::max(diameter, vec_distance_sqr(a, b));
			const fl margin = 2 * std::sqrt(diameter) + 1;
			for (const vec &xyz : probe.coords)
				for (sz i = 0; i < 3; ++i) {
					bounds[i].begin = std::min(bounds[i].begin, xyz[i] - margin);
					bounds[i].end = std::max(bounds[i].end, xyz[i] + margin);
				}
		}
		return bounds;
	}

  private:
	model &m;
	std::vector<rotation> rotors;
};

int main(int argc, char **argv) {
	try {
		if (argc != 6)
			throw std::runtime_error("Usage: rotor_partition receptor.pdbqt ligand.pdbqt native_pool.pdbqt "
									 "even_angle_count temperature_K");
		const int steps = std::stoi(argv[4]);
		const fl temperature = std::stod(argv[5]);
		if (steps < 4 || steps % 2 || !std::isfinite(temperature) || temperature <= 0)
			throw std::runtime_error(
				"Use a positive temperature and an even angular grid of at least four points");
		std::ifstream file(argv[3]);
		std::vector<std::string> models;
		std::string line, block;
		while (std::getline(file, line)) {
			if (line.compare(0, 6, "MODEL ") == 0) {
				if (!block.empty())
					throw std::runtime_error("Nested MODEL");
				block = line + "\n";
			} else if (!block.empty()) {
				block += line + "\n";
				if (line == "ENDMDL") {
					models.push_back(block);
					block.clear();
				}
			}
		}
		if (models.empty() || !block.empty())
			throw std::runtime_error("No complete MODEL records");
		model rec = parse_receptor_pdbqt(argv[1]), m = rec;
		m.append(parse_ligand_pdbqt_from_file(argv[2], atom_type::XS));
		for (sz i = 0; i < m.num_movable_atoms(); ++i)
			if (m.movable_atom(i).el == EL_TYPE_Dummy ||
				(m.movable_atom(i).ad >= AD_TYPE_CG0 && m.movable_atom(i).ad <= AD_TYPE_CG3))
				throw std::runtime_error("Unsupported macrocycle or dummy atom");
		const auto scoring = std::make_shared<ScoringFunction>(
			SF_VINA, flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5 * 0.05846 / 0.1 - 1});
		const precalculate_byatom pairs(*scoring, m);
		non_cache grid(rec, rotor_partition::scan_bounds(m, models), scoring, 0);
		grid.prepare_for_model(m);
		rotor_partition sampler(m);
		std::cout << std::setprecision(12) << "model,rotors,inter,intra,work,half_grid_work\n";
		for (sz i = 0; i < models.size(); ++i) {
			replace_fixed_pose_coordinates(m, models[i]);
			sampler.evaluate(grid, pairs, steps, 0.00198720425864 * temperature, i + 1);
		}
	} catch (const std::exception &e) {
		std::cerr << e.what() << '\n';
		return 1;
	}
}
