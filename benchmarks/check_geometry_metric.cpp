// Verify every ligand coordinate against Cartesian finite differences,
// including branches, macrocycle closure dummies, and transformed conformers.
#include "parse_pdbqt.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>

namespace {
flv finite_difference(model& m, const conf& base) {
	const auto size = m.get_size();
	flv result(size.num_degrees_of_freedom(), 1), step(result.size(), 0);
	conf trial(size);
	const fl h = 1e-4;
	// Translation entries have analytic value 1 in both algorithms.
	for (sz k = 3; k < result.size(); ++k) {
		step[k] = 1;
		trial.set_incremented(base, step, h);
		m.set(trial);
		const auto plus = m.get_chemical_heavy_atom_movable_coords();
		trial.set_incremented(base, step, -h);
		m.set(trial);
		const auto minus = m.get_chemical_heavy_atom_movable_coords();
		fl sum = 0;
		VINA_FOR_IN(i, plus) sum += vec_distance_sqr(plus[i], minus[i]);
		const fl g = sum / (4 * h*h * plus.size());
		result[k] = 1 / std::min(fl(64), std::max(fl(0.25), g));
		step[k] = 0;
	}
	m.set(base);
	return result;
}

template<typename F>
double microseconds_per_metric(F compute, model& m, const conf& c) {
	const unsigned repeats = 2000;
	volatile fl checksum = 0;
	const auto start = std::chrono::steady_clock::now();
	for (unsigned i = 0; i < repeats; ++i)
		checksum += compute(m, c)[i % (6 + c.ligands.front().torsions.size())];
	if (!std::isfinite(checksum)) throw std::runtime_error("Non-finite metric");
	const auto elapsed = std::chrono::steady_clock::now() - start;
	return std::chrono::duration<double, std::micro>(elapsed).count() / repeats;
}
}

int main(int argc, char** argv) {
	if (argc < 2) throw std::runtime_error("Supply at least one ligand PDBQT fixture");
	for (int f = 1; f < argc; ++f) {
		model m = parse_ligand_pdbqt_from_file(argv[f], atom_type::XS);
		const conf initial = m.get_initial_conf();
		const auto analytic = [](model& model, const conf& c){ return model.geometry_inverse_metric(c); };
		fl max_abs = 0;
		for (unsigned i = 0; i < 24; ++i) {
			conf c = initial;
			auto& l = c.ligands.front();
			l.rigid.position += vec(10 * i, -4.5 * i, 1.5 * i);
			quaternion_increment(l.rigid.orientation, vec(0.13*i, -0.19*i, 0.045*i));
			VINA_FOR_IN(k, l.torsions) l.torsions[k] = (i ? 2.7*std::sin((i+1)*(k+1.)) : 0.);
			const flv expected = finite_difference(m,c);
			const auto expected_coords = m.get_chemical_heavy_atom_movable_coords();
			m.set(initial); // simulate a stale model from a rejected pose
			const flv computed = analytic(m,c);
			const auto actual_coords = m.get_chemical_heavy_atom_movable_coords();
			if (actual_coords.size() != expected_coords.size() || computed.size() != expected.size())
				throw std::runtime_error("Wrong coordinate or heavy-atom count");
			VINA_FOR_IN(j, actual_coords)
				if (vec_distance_sqr(actual_coords[j],expected_coords[j]) != 0)
					throw std::runtime_error("Metric did not leave model at requested conformation");
			VINA_FOR_IN(k, computed) {
				const fl diff = std::abs(computed[k] - expected[k]);
				if (!std::isfinite(diff) || diff > 2e-7)
					throw std::runtime_error("Analytic metric disagrees with coordinate probes");
				max_abs = std::max(max_abs, diff);
			}
		}
		const double fd_us = microseconds_per_metric(finite_difference, m, initial);
		const double exact_us = microseconds_per_metric(analytic, m, initial);
		std::cout << argv[f] << " dof=" << m.get_size().num_degrees_of_freedom()
		          << " max_abs=" << max_abs << " finite_us=" << fd_us
		          << " analytic_us=" << exact_us << '\n';
	}
}
