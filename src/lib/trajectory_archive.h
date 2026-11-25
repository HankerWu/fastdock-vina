#ifndef VINA_TRAJECTORY_ARCHIVE_H
#define VINA_TRAJECTORY_ARCHIVE_H

#include "coords.h"
#include "model.h"

// This object owns only recorded evidence. It never mutates the search model,
// writes to the historical output container, evaluates energy or advances RNG.
struct trajectory_archive {
	output_container accepted;
	std::vector<vecv> samples;
	sz native_evaluations = 0;
	sz native_refinements = 0;

	void capture(const model& m, const output_type& pose, fl min_rmsd, sz limit) {
		if (accepted.size() == limit && pose.e >= accepted.back().e)
			return;

		vecv heavy = m.get_chemical_heavy_atom_movable_coords();
		szv displaced;
		VINA_FOR_IN(i, accepted) {
			const fl rmsd = rmsd_upper_bound(heavy, accepted[i].coords);
			if (rmsd != 0 && rmsd >= min_rmsd)
				continue;
			if (accepted[i].e <= pose.e)
				return;
			displaced.push_back(i);
		}
		for (auto it = displaced.rbegin(); it != displaced.rend(); ++it)
			accepted.erase(accepted.begin() + *it);
		if (accepted.size() == limit)
			accepted.pop_back();

		const auto position = std::upper_bound(accepted.begin(), accepted.end(), pose.e,
			[](fl energy, const output_type& existing) { return energy < existing.e; });
		auto item = accepted.insert(position, new output_type(pose.c, pose.e));
		item->coords = std::move(heavy);
		item->archive_from_accepted = true;
	}
};

inline bool identical_archive_conformation(const conf& a, const conf& b) {
	// Archive docking accepts one ligand and no flexible receptor. Equal
	// conformation parameters also preserve any explicit hydrogen coordinates.
	const ligand_conf& first = a.ligands.front();
	const ligand_conf& second = b.ligands.front();
	const vec& x = first.rigid.position;
	const vec& y = second.rigid.position;
	return x[0] == y[0] && x[1] == y[1] && x[2] == y[2] &&
	       first.rigid.orientation == second.rigid.orientation && first.torsions == second.torsions;
}

inline sz append_accepted_minima(output_container& poses, std::vector<trajectory_archive>& archives) {
	sz total = 0;
	for (auto& archive : archives) {
		auto& accepted = archive.accepted;
		total += accepted.size();
		for (sz i = 0; i < accepted.size();) {
			const conf& candidate = accepted[i].c;
			if (std::any_of(poses.begin(), poses.end(), [&](const output_type& pose) {
			    return identical_archive_conformation(pose.c, candidate);
			})) {
				++i;
			} else {
				poses.transfer(poses.end(), accepted.begin() + i, accepted);
			}
		}
	}
	return total;
}

inline sz archived_trajectory_support(const vecv& coords, const std::vector<trajectory_archive>& archives) {
	return std::count_if(archives.begin(), archives.end(), [&](const trajectory_archive& history) {
		return std::any_of(history.samples.begin(), history.samples.end(), [&](const vecv& sample) {
			return rmsd_upper_bound(coords, sample) <= 2.0;
		});
	});
}

#endif
