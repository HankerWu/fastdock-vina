#include "amide_search_grid.h"
#include "non_cache.h"
#include "precalculate.h"

#include <iostream>

static void require(bool condition, const char* message) {
	if (!condition) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
	require(argc == 2, "Pass the prepared amide fixture fdv-second/prepared/6WQ8_ligand.pdbqt");
	constexpr fl step = 1e-5;
	const amide_search_site two{0, 1, 2, {{3, 4}}, 2, 5};
	const amide_search_site swapped{0, 1, 2, {{4, 3}}, 2, 5};
	const amide_search_site one{0, 1, 2, {{3, 3}}, 1, 5};
	vecv xyz{vec(1, 0, -0.5), vec(0, 0, 0), vec(0, 0, 1.35),
	         vec(1, 0, 1.9), vec(-1, 0, 1.92)};
	require(two.bias(xyz) == 0, "both planar wells should have zero energy");
	xyz[3] = vec(0, 1, 1.9);
	require(std::abs(one.bias(xyz) - 5) < 1e-12 && std::abs(two.bias(xyz) - 2.5) < 1e-12,
	        "one-arm barrier and two-arm normalization");
	xyz[0] += vec(-0.2, 0.3, 0.1);
	xyz[2] += vec(0.4, -0.2, 0.07);
	xyz[4] += vec(0.1, -0.5, 0.3);
	for (const auto& site : {one, two, swapped}) {
		vecv gradient(xyz.size(), zero_vec);
		const fl energy = site.bias(xyz, &gradient);
		require(std::abs(energy - site.bias(xyz)) < 1e-12, "scalar vs gradient energy");
		for (sz atom = 0; atom < xyz.size(); ++atom)
			for (sz k = 0; k < 3; ++k) {
				auto plus = xyz, minus = xyz;
				plus[atom][k] += step;
				minus[atom][k] -= step;
				const fl numeric = (site.bias(plus) - site.bias(minus)) / (2 * step);
				require(std::abs(numeric - gradient[atom][k]) < 1e-7, "Cartesian gradient");
			}
	}
	require(two.bias(xyz) == swapped.bias(xyz), "N-arm permutation changed the energy");
	auto m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
	const conf origin = m.get_initial_conf();
	m.set(origin);
	const std::string row = "2 1 3 4 17 5";
	for (const auto& invalid : std::vector<std::string>{"", row + "\n" + row, "2 1 3 4 5",
	         "1 2 3 4 17 5", "2 1 4 3 17 5", "2 1 3 4 4 5", "2 1 3 4 17 -1", "2 1 3 4 17 5x",
	         "2 1 3 4 17 inf", "2 1 3 4 99999 5", "2 1 3 4 17 nan"}) {
		std::istringstream input(invalid);
		bool rejected = false;
		try { read_search_amides(input, m); }
		catch (const std::invalid_argument&) { rejected = true; }
		require(rejected, "accepted an invalid amide record");
	}
	std::istringstream good(row);
	auto sites = read_search_amides(good, m);
	require(sites.size() == 1 && sites.front().arm_count == 2, "6WQ8 site recognition");

	const auto sf = std::make_shared<ScoringFunction>(SF_VINA,
		flv{-0.035579, -0.005156, 0.840245, -0.035069, -0.587439, 50, 5 * 0.05846 / 0.1 - 1});
	const precalculate_byatom pairs(*sf, m);
	grid_dims bounds;
	for (auto& axis : bounds) { axis.begin = -80; axis.end = 80; axis.n_voxels = 400; }
	non_cache base(m, bounds, sf, 1e6);
	base.prepare_for_model(m);
	const amide_search_grid guide(base, sites);
	const vec cap(1000, 1000, 1000);
	flv g0(m.get_size().num_degrees_of_freedom()), g1(g0.size());
	const fl e0 = m.eval_deriv(pairs, base, cap, g0);
	const fl e1 = m.eval_deriv(pairs, guide, cap, g1);
	const auto bias = [&] { return guide.eval(m, 1000) - base.eval(m, 1000); };
	require(std::abs(bias() - (e1 - e0)) < 1e-10, "integrated scalar vs gradient energy");
	require(base.eval_intra(m, 1000) == guide.eval_intra(m, 1000), "ligand reference should be native");
	for (sz i = 0; i < g0.size(); ++i) {
		flv basis(g0.size(), 0);
		basis[i] = 1;
		conf shifted = origin;
		shifted.set_incremented(origin, basis, step);
		m.set(shifted);
		const fl plus = bias();
		shifted.set_incremented(origin, basis, -step);
		m.set(shifted);
		const fl minus = bias();
		require(std::abs((plus - minus) / (2 * step) - (g1[i] - g0[i])) < 1e-6,
		        "rigid body or torsion gradient");
	}
	std::cout << "Amide Cartesian and conformation gradients, planar wells, atom validation and native reference passed\n";
}
