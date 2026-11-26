#include "anchored_search_grid.h"
#include "non_cache.h"
#include "parse_pdbqt.h"
#include "precalculate.h"

#include <iostream>
#include <stdexcept>

void check(bool value, const char* message) {
	if (!value) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
	check(argc == 2, "Pass the 1IEP ligand path");
	const fl step = 1e-5;
	for (const fl radius : {0., 2.}) {
		for (const vec offset : {vec(0, 0, 0), vec(2, 0, 0), vec(0.3, -0.8, 0.2), vec(10, 2, -4)}) {
			vec analytic;
			anchored_search_grid::penalty(offset, &analytic, radius);
			for (sz axis = 0; axis < 3; ++axis) {
				vec plus = offset, minus = offset;
				plus[axis] += step;
				minus[axis] -= step;
				const fl numeric = (anchored_search_grid::penalty(plus, nullptr, radius) -
				                    anchored_search_grid::penalty(minus, nullptr, radius)) / (2 * step);
				check(std::abs(analytic[axis] - numeric) < 1e-7, "Cartesian gradient");
			}
		}
		check(anchored_search_grid::penalty(vec(radius, 0, 0), nullptr, radius) == 0, "minimum");
	}

	model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
	const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
		flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5 * 0.05846 / 0.1 - 1});
	const precalculate_byatom pairs(*sf, m);
	grid_dims bounds;
	for (auto& axis : bounds) { axis.begin = -80; axis.end = 80; axis.n_voxels = 400; }
	non_cache base(m, bounds, sf, 1e6);
	base.prepare_for_model(m);
	const conf origin = m.get_initial_conf();
	m.set(origin);
	const auto indices = ligand_heavy_indices_by_serial(m);
	const sz first = indices.at(3), second = indices.at(39);
	const vec a(0.25, -0.52, 0.81), b(-1.9, 0.4, 0.23);
	const anchor_assignment first_map = {{first, m.movable_coords(first) - a, 0},
	                                     {second, m.movable_coords(second) - b, 2}};
	const vec c(1.3, 0.6, -0.3);
	const anchor_assignment second_map = {first_map[0], {second, m.movable_coords(second) - c, 0}};
	const anchor_assignment permuted = {first_map[1], first_map[0]};
	const fl mean1 = (anchored_search_grid::penalty(a) + anchored_search_grid::penalty(b, nullptr, 2)) / 2;
	const fl mean2 = (anchored_search_grid::penalty(a) + anchored_search_grid::penalty(c)) / 2;
	const fl ensemble = -0.6 * std::log((std::exp(-mean1/0.6) + std::exp(-mean2/0.6)) / 2);
	const std::vector<std::vector<anchor_assignment>> cases = {
		{first_map}, {first_map, permuted}, {first_map, second_map},
		{second_map, permuted, first_map}};
	const vec cap(1000, 1000, 1000);
	flv g0(m.get_size().num_degrees_of_freedom()), g1(g0.size());
	const fl e0 = m.eval_deriv(pairs, base, cap, g0);
	for (sz k = 0; k < cases.size(); ++k) {
		const anchored_search_grid guide(base, cases[k]);
		m.set(origin);
		const fl expected = k < 2 ? mean1 : ensemble;
		const fl e1 = m.eval_deriv(pairs, guide, cap, g1);
		check(std::abs(e1 - e0 - expected) < 1e-12, "ensemble derivative energy");
		check(std::abs(guide.eval(m, 1000) - base.eval(m, 1000) - expected) < 1e-12, "ensemble scalar energy");
		check(guide.eval_intra(m, 1000) == base.eval_intra(m, 1000), "ligand-only reference");
		for (sz i = 0; i < g0.size(); ++i) {
			flv basis(g0.size(), 0);
			basis[i] = 1;
			conf shifted = origin;
			shifted.set_incremented(origin, basis, step);
			m.set(shifted);
			const fl plus = guide.eval(m, 1000) - base.eval(m, 1000);
			shifted.set_incremented(origin, basis, -step);
			m.set(shifted);
			const fl minus = guide.eval(m, 1000) - base.eval(m, 1000);
			check(std::abs((plus - minus)/(2*step) - (g1[i]-g0[i])) < 1e-6, "torsion/rigid gradient");
		}
	}
	m.set(origin);
	const anchor_assignment distant = {{first, m.movable_coords(first) + vec(1e5, 0, 0), 0}};
	const anchored_search_grid separated(base, std::vector<anchor_assignment>{first_map, distant});
	check(std::abs(separated.eval(m, 1000) - base.eval(m, 1000) - mean1 - 0.6*std::log(2.)) < 1e-12,
	      "stable large-energy separation");
	std::cout << "Single and alternative point/shell energies, order, multiplicity, stability and pose gradients passed\n";
}
