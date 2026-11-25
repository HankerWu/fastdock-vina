#include "cache.h"
#include "monte_carlo.h"
#include "parse_pdbqt.h"
#include "precalculate.h"
#include "scoring_function.h"
#include "trajectory_archive.h"

#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
	if (!condition)
		throw std::runtime_error(message);
}

void check_queue(model m) {
	trajectory_archive archive;
	const conf initial = m.get_initial_conf();
	const auto capture = [&](fl shift, fl energy) {
		output_type pose(initial, energy);
		pose.c.ligands[0].rigid.position[0] += shift;
		m.set(pose.c);
		archive.capture(m, pose, 1.0, 2);
	};

	capture(0, -1);
	capture(1.5, -2);
	capture(4, 0);
	require(archive.accepted.size() == 2 && archive.accepted.back().e == -1,
	        "A full archive accepted a worse distinct minimum");

	// This better pose replaces both conflicting minima, keeping the queue separated.
	capture(0.75, -3);
	require(archive.accepted.size() == 1 && archive.accepted.front().e == -3,
	        "A better overlapping pose did not displace both worse minima");
	capture(0.75, -3);
	capture(4, -2);
	capture(8, -4);
	require(archive.accepted.size() == 2 && archive.accepted[0].e == -4 && archive.accepted[1].e == -3,
	        "The archive did not retain the two best separated minima");

	output_container historical;
	historical.push_back(new output_type(archive.accepted[0].c, -8));
	std::vector<trajectory_archive> archives{archive, archive};
	require(append_accepted_minima(historical, archives) == 4 && historical.size() == 2,
	        "Raw accounting or exact duplicate merge is incorrect");
	require(!historical[0].archive_from_accepted && historical[0].e == -8 &&
	        historical[1].archive_from_accepted && historical[1].e == -3,
	        "The merge lost historical priority or accepted-seed provenance");

	archives.resize(3);
	archives[0].samples = {{vec(2.0, 0, 0)}, {vec(1, 0, 0)}};
	archives[1].samples = {{vec(2.0, 0, 0)}};
	archives[2].samples = {{vec(2.000001, 0, 0)}};
	require(archived_trajectory_support({vec(0, 0, 0)}, archives) == 2,
	        "Votes were repeated, excluded at 2.0 Å, or included beyond 2.0 Å");
}

void check_search(const model& original) {
	const ScoringFunction scoring(SF_VINA, {-0.035579, -0.005156, 0.840245, -0.035069,
	                                       -0.587439, 50, 0.05846});
	const precalculate scores(scoring);
	const precalculate_byatom ligand_scores(scoring, original);
	grid_dims dims;
	const vec center(15.190, 53.903, 16.917);
	VINA_FOR(i, 3) {
		dims[i].begin = center[i] - 10;
		dims[i].end = center[i] + 10;
		dims[i].n_voxels = 40;
	}
	cache grid(dims, 1e6);
	szv_grid index(original, szv_grid_dims(dims), scores.cutoff_sqr());
	grid.populate(original, scores, index, original.get_movable_atom_types(atom_type::XS));

	monte_carlo search;
	search.max_evals = 1800;
	search.local_steps = unsigned((25 + original.num_movable_atoms()) / 3);
	search.hunt_cap = vec(10, 10, 10);
	search.num_saved_mins = 3;
	search.min_rmsd = 1;

	model native(original);
	model observed(original);
	output_container native_poses;
	output_container observed_poses;
	trajectory_archive archive;
	rng native_rng(42);
	rng observed_rng(42);
	search(native, native_poses, ligand_scores, grid, grid.corner1(), grid.corner2(), nullptr, native_rng);
	search(observed, observed_poses, ligand_scores, grid, grid.corner1(), grid.corner2(), nullptr, observed_rng, &archive);
	require(native_poses.size() == observed_poses.size(), "Archiving changed the historical pool size");
	VINA_FOR_IN(i, native_poses) {
		const output_type& first = native_poses[i];
		const output_type& second = observed_poses[i];
		require(first.e == second.e && identical_archive_conformation(first.c, second.c) &&
		        rmsd_upper_bound(first.coords, second.coords) == 0,
		        "Archiving changed a historical Monte Carlo minimum");
	}
	require(rmsd_upper_bound(native.get_chemical_heavy_atom_movable_coords(),
	                         observed.get_chemical_heavy_atom_movable_coords()) == 0,
	        "Archive snapshots changed the terminal search model");
	require(native_rng() == observed_rng(), "Archiving advanced the trajectory RNG");
	require(archive.native_evaluations > 0 && archive.native_refinements > 0 &&
	        !archive.accepted.empty() && archive.accepted.size() <= search.num_saved_mins,
	        "Archiving did not record the bounded observed search");

	for (const auto& pose : archive.accepted) {
		observed.set(pose.c);
		require(rmsd_upper_bound(pose.coords, observed.get_chemical_heavy_atom_movable_coords()) == 0,
		        "Archived coordinates do not match the accepted conformation");
	}
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 3)
		throw std::runtime_error("Pass the receptor and ligand PDBQT fixture paths");
	model m = parse_receptor_pdbqt(argv[1]);
	m.append(parse_ligand_pdbqt_from_file(argv[2], atom_type::XS));
	check_queue(m);
	check_search(m);
	std::cout << "Archive capacity, duplicate priority, inclusive unique votes, search and RNG parity passed.\n";
}
