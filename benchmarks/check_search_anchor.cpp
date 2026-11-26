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
	const anchored_search_grid guide(base, {{first, m.movable_coords(first) - a, 0},
	                                        {second, m.movable_coords(second) - b, 2}});
	const vec cap(1000, 1000, 1000);
	flv g0(m.get_size().num_degrees_of_freedom()), g1(g0.size());
	const fl e0 = m.eval_deriv(pairs, base, cap, g0);
	const fl e1 = m.eval_deriv(pairs, guide, cap, g1);
	const fl mean = (anchored_search_grid::penalty(a) + anchored_search_grid::penalty(b, nullptr, 2)) / 2;
	check(std::abs(e1 - e0 - mean) < 1e-12, "mixed-site normalization/energy");
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
	std::cout << "Point, shell and mixed-site energies, Cartesian and pose gradients passed\n";
}
