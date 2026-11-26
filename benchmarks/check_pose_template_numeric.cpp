#include "non_cache.h"
#include "parse_pdbqt.h"
#include "pose_template.h"
#include "precalculate.h"
#include "quasi_newton.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

void require(bool condition, const char* message) {
	if (!condition)
		throw std::runtime_error(message);
}

std::string atom_record(sz serial, const char* type, fl x) {
	std::ostringstream atom;
	atom << "ATOM  " << std::setw(5) << serial << "  X   UNL     1    "
	     << std::fixed << std::setprecision(3) << std::setw(8) << x
	     << "   0.000   0.000  1.00  0.00     0.000 " << type << '\n';
	return atom.str();
}

model parse_root(const std::string& atoms) {
	return parse_ligand_pdbqt_from_string("ROOT\n" + atoms + "ENDROOT\nTORSDOF 0\n", atom_type::XS);
}

non_cache empty_receptor_grid(const model& ligand, const std::shared_ptr<const ScoringFunction>& scoring) {
	grid_dims bounds;
	for (auto& axis : bounds) {
		axis.begin = -80;
		axis.end = 80;
		axis.n_voxels = 400;
	}
	non_cache grid(ligand, bounds, scoring, 1e6);
	grid.prepare_for_model(ligand);
	return grid;
}

void check_cartesian() {
	const pose_template::ligand_atoms atoms{{0, "C"}, {1, "C"}, {2, "O"}, {3, "N"}};
	std::istringstream source("C 1 0 0\nC 3 1 0\nO 0 0 0\nCl 0 0 0\nH 0 0 0\nG0 4 4 4\n");
	pose_template reference(source, atoms);
	vecv xyz{{0, 0, 0}, {2.5, 0.8, 0.1}, {0.5, 0.4, 0.3}, {0, 0, 0}};
	vecv gradient(xyz.size(), vec(0.7, -0.5, 0.3));
	const fl sigma = 1.5;
	const auto gaussian = [sigma](fl squared) { return std::exp(-squared / (2 * sigma * sigma)); };
	const fl expected = (gaussian(1) + gaussian(0.30) + gaussian(0.50)) / 4;
	require(std::abs(reference.overlap(xyz, &gradient) - expected) < 1e-14,
	        "Bidirectional overlap or unmatched-element normalization differs from the formula");
	const fl step = 1e-5;
	const vec prior(0.7, -0.5, 0.3);
	for (sz i = 0; i < xyz.size(); ++i) {
		for (sz axis = 0; axis < 3; ++axis) {
			const fl origin = xyz[i][axis];
			xyz[i][axis] = origin + step;
			const fl plus = -6 * reference.overlap(xyz);
			xyz[i][axis] = origin - step;
			const fl minus = -6 * reference.overlap(xyz);
			xyz[i][axis] = origin;
			require(std::abs((plus - minus) / (2 * step) - (gradient[i][axis] - prior[axis])) < 1e-9,
			        "Template Cartesian derivative disagrees with centered finite differences");
		}
	}

	std::istringstream coincident("C 0 0 0\nC 0 0 0\n");
	pose_template tied(coincident, {{0, "C"}, {1, "C"}});
	const vecv duplicate{{0, 0, 0}, {0, 0, 0}};
	vecv tie_gradient(duplicate.size(), zero_vec);
	require(tied.overlap(duplicate, &tie_gradient) == 1 && tie_gradient[0].norm_sqr() == 0
	        && tie_gradient[1].norm_sqr() == 0, "Coincident atoms must give a stationary maximum");
	std::istringstream unmatched("Xe 0 0 0\n");
	pose_template no_match(unmatched, atoms);
	require(no_match.overlap(xyz) == 0, "A template with no matching elements must give zero overlap");
}

void check_contact_footprints() {
	auto site = [](fl x, sz element) {
		atom a;
		a.coords = vec(x, 0, 0);
		a.el = element;
		return a;
	};
	const atomv receptor{site(0, EL_TYPE_C), site(12, EL_TYPE_C), site(24, EL_TYPE_C),
	                     site(12, EL_TYPE_H), site(12, EL_TYPE_Dummy)};
	std::istringstream nearby("C 0 0 0\nC 12 0 0\n");
	pose_template comparison(nearby, {{0, "C"}, {1, "C"}}, &receptor);
	require(std::abs(comparison.overlap({{12,0,0},{24,0,0}}) - 1./3) < 1e-14,
	        "Reference-only and query-only contacts must both enter the denominator");

	const atomv polar{site(0, EL_TYPE_O)};
	std::istringstream donor("N 0 0 0\n");
	pose_template two_channels(donor, {{0, "O"}}, &polar);
	require(std::abs(two_channels.overlap({{3.75,0,0}}) - 5./6) < 1e-14,
	        "Polar contact ramp or continuous Tanimoto normalization is wrong");

	std::istringstream repeated("N 0 0 0\n");
	pose_template one_sided(repeated, {{0, "C"}}, &polar);
	require(one_sided.overlap({{0,0,0}}) == 0.5, "A one-sided polar channel is active and scores zero");
	std::istringstream distant("N 100 0 0\n");
	pose_template empty(distant, {{0, "O"}}, &polar);
	require(empty.overlap({{100,0,0}}) == 0, "No contacts in either channel must score zero");
}

void check_parsed_elements() {
	const auto molecule = parse_root(atom_record(1, "C", 0) + atom_record(2, "CG0", 1.5)
	      + atom_record(3, "G0", 2) + atom_record(4, "W", 3) + atom_record(5, "HD", 4)
	      + atom_record(6, "Se", 5) + atom_record(7, "S", 6) + atom_record(8, "Zn", 7)
	      + atom_record(9, "Fe", 8) + atom_record(10, "NA", 9));
	const pose_template::ligand_atoms expected{{0, "C"}, {1, "C"}, {5, "Se"}, {6, "S"},
	                                          {7, "Zn"}, {8, "Fe"}, {9, "N"}};
	require(molecule.chemical_heavy_atom_elements() == expected,
	        "Chemical types must omit dummies/hydrogen and distinguish metals, Se, and CG0");
}

void check_optimizer(const std::string& ligand_path) {
	const auto vina = std::make_shared<ScoringFunction>(SF_VINA, flv{-0.035579, -0.005156, 0.840245, -0.035069,
	                                    -0.587439, 50, 5 * 0.05846 / 0.1 - 1});
	const vec cap(1000, 1000, 1000);
	model ligand = parse_ligand_pdbqt_from_file(ligand_path, atom_type::XS);
	const precalculate_byatom pairs(*vina, ligand);
	const non_cache grid = empty_receptor_grid(ligand, vina);
	const conf initial = ligand.get_initial_conf();
	ligand.set(initial);
	std::ostringstream template_rows;
	for (const auto& atom : ligand.chemical_heavy_atom_elements()) {
		const vec& pos = ligand.movable_coords(atom.first);
		template_rows << atom.second << ' ' << pos[0] + 0.25 << ' ' << pos[1] + 0.34
		              << ' ' << pos[2] - 0.45 << '\n';
	}
	std::istringstream source(template_rows.str());
	pose_template reference(source, ligand.chemical_heavy_atom_elements());
	flv native_gradient(ligand.get_size().num_degrees_of_freedom());
	flv full_gradient(native_gradient.size());
	const fl native = ligand.eval_deriv(pairs, grid, cap, native_gradient);
	const fl guided = ligand.eval_deriv(pairs, grid, cap, full_gradient, &reference);
	require(std::abs(guided - native + 6 * ligand.template_overlap(reference)) < 1e-12,
	        "Template gradient path must report the same raw term as rank/select");
	const fl step = 1e-5;
	for (sz dof = 0; dof < full_gradient.size(); ++dof) {
		flv basis(full_gradient.size(), 0);
		basis[dof] = 1;
		conf perturbed = initial;
		perturbed.set_incremented(initial, basis, step);
		ligand.set(perturbed);
		const fl plus = -6 * ligand.template_overlap(reference);
		perturbed.set_incremented(initial, basis, -step);
		ligand.set(perturbed);
		const fl minus = -6 * ligand.template_overlap(reference);
		require(std::abs((plus - minus) / (2 * step) - (full_gradient[dof] - native_gradient[dof])) < 1e-5,
		        "Template contribution to pose translation/orientation/torsions disagrees with finite differences");
	}

	model point = parse_root(atom_record(1, "C", 0));
	const precalculate_byatom point_pairs(*vina, point);
	const non_cache point_grid = empty_receptor_grid(point, vina);
	std::istringstream target("C 2 0 0\n");
	pose_template guide(target, point.chemical_heavy_atom_elements());
	output_type result(point.get_initial_conf(), 0);
	quasi_newton refine(point.get_size(), 15);
	int calls = 0;
	refine(point, point_pairs, point_grid, result, cap, calls, &guide);
	require(std::abs(point.movable_coords(0)[0] - 2) < 1e-5 && std::abs(result.e + 6) < 1e-8 && calls > 1,
	        "Native BFGS with the template term must reach the attractive minimum");
}

int main(int argc, char** argv) {
	require(argc == 2, "Pass the 1IEP ligand PDBQT path");
	check_cartesian();
	check_contact_footprints();
	check_parsed_elements();
	check_optimizer(argv[1]);
	std::cout << "Template overlap, exact elements, Cartesian and pose derivatives, and moving BFGS minimum passed.\n";
}
