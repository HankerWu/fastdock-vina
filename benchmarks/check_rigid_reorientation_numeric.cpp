#include "mutate.h"
#include "parse_pdbqt.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
	if (!condition)
		throw std::runtime_error(message);
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 2)
		throw std::runtime_error("Pass a ligand PDBQT fixture with chemical heavy atoms and torsions");
	model m = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
	conf accepted = m.get_initial_conf();
	require(!accepted.ligands.front().torsions.empty(), "Fixture must have ligand torsions");
	accepted.ligands.front().torsions.front() = 1.1;
	accepted.ligands.front().rigid.position += vec(8, 2, -4);
	quaternion_increment(accepted.ligands.front().rigid.orientation, vec(0, 0.7, 0.9));
	m.set(accepted);
	const auto before = m.get_chemical_heavy_atom_movable_coords();
	const auto centroid = [](const vecv& coordinates) {
		vec center(0, 0, 0);
		for (const auto& atom : coordinates)
			center += atom;
		return (1.0 / fl(coordinates.size())) * center;
	};
	require(before.size() >= 3, "Fixture needs an extended ligand");
	const vec center = centroid(before);
	conf stale = accepted;
	stale.ligands.front().rigid.position += vec(-11, -15, 17);
	stale.ligands.front().torsions.front() = -1.2;
	quaternion_increment(stale.ligands.front().rigid.orientation, vec(0.6, -0.5, 0));

	unsigned global_rotations = 0;
	fl max_displacement = 0;
	fl max_center_error = 0;
	fl max_pair_error = 0;
	rng generator(42);
	for (unsigned step = 0; step < 2048; ++step) {
		m.set(stale); // Emulate a rejected trial's model storage.
		conf proposal = accepted;
		if (!mutate_conf(proposal, m, 2, generator, true))
			continue;
		++global_rotations;
		require(proposal.ligands.front().torsions == accepted.ligands.front().torsions,
		        "Rigid reorientation changed the accepted torsions");
		m.set(proposal);
		const auto after = m.get_chemical_heavy_atom_movable_coords();
		max_center_error = std::max(max_center_error, std::sqrt(vec_distance_sqr(center, centroid(after))));
		VINA_FOR_IN(i, before) {
			max_displacement = std::max(max_displacement, std::sqrt(vec_distance_sqr(before[i], after[i])));
			VINA_FOR(j, i)
				max_pair_error = std::max(max_pair_error,
				    std::abs(vec_distance_sqr(before[i], before[j]) - vec_distance_sqr(after[i], after[j])));
		}
	}
	require(global_rotations > 0 && max_displacement > 2, "Fixture failed to exercise nonlocal rotation");
	require(max_center_error < 1e-10, "Proposal moved the accepted chemical heavy centroid");
	require(max_pair_error < 1e-9, "Proposal deformed the accepted internal geometry");
	std::cout << "Global rotations=" << global_rotations << " max centroid drift=" << max_center_error
	          << " A max squared pair-distance drift=" << max_pair_error << " A^2\n";

	unsigned local_rotations = 0;
	fl native_centroid_drift = 0;
	fl centered_centroid_drift = 0;
	rng native_generator(43), centered_generator(43);
	for (unsigned step = 0; step < 2048; ++step) {
		conf native = accepted;
		conf centered = accepted;
		m.set(stale);
		require(!mutate_conf(native, m, 2, native_generator, false), "Native local proposal was global");
		m.set(stale);
		require(!mutate_conf(centered, m, 2, centered_generator, false, false, true),
		        "Centered local proposal was global");
		const auto& source = accepted.ligands.front();
		const auto& old_pivot = native.ligands.front();
		const auto& new_pivot = centered.ligands.front();
		require(native_generator == centered_generator, "Centered proposal changed the RNG stream");
		require(old_pivot.rigid.orientation == new_pivot.rigid.orientation &&
		        old_pivot.torsions == new_pivot.torsions, "Centered proposal changed rotation or torsions");
		if (old_pivot.rigid.orientation == source.rigid.orientation) {
			require(vec_distance_sqr(old_pivot.rigid.position, new_pivot.rigid.position) == 0,
			        "Centered mode changed a non-rotational proposal");
			continue;
		}
		++local_rotations;
		m.set(native);
		native_centroid_drift = std::max(native_centroid_drift,
		    std::sqrt(vec_distance_sqr(center, m.ligand_heavy_centroid(0))));
		m.set(centered);
		centered_centroid_drift = std::max(centered_centroid_drift,
		    std::sqrt(vec_distance_sqr(center, m.ligand_heavy_centroid(0))));
	}
	require(local_rotations > 0 && native_centroid_drift > 0.1, "Fixture did not expose root-rotation drift");
	require(centered_centroid_drift < 1e-10, "Centered local proposal drifted from the accepted centroid");
	std::cout << "Local rotations=" << local_rotations << " max native centroid drift="
	          << native_centroid_drift << " A max centered centroid drift="
	          << centered_centroid_drift << " A\n";
}
