#include "parse_pdbqt.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <tuple>

namespace {

struct signature {
	std::set<std::pair<unsigned, unsigned>> bonds;
	std::set<std::tuple<unsigned, unsigned, sz>> internal;
	std::vector<sz> xs;
};

signature topology(const model& ligand_model) {
	signature result;
	const auto ligand = ligand_model.get_ligand(0);
	std::vector<unsigned> serial(ligand_model.num_atoms());
	for (const auto& line : ligand.cont)
		if (line.second)
			serial[*line.second] = std::stoul(line.first.substr(6, 5));
	VINA_FOR(i, ligand_model.num_atoms()) {
		const auto atom = ligand_model.get_atom(i);
		result.xs.push_back(atom.xs);
		for (const auto& bond : atom.bonds) {
			const unsigned other = serial[bond.connected_atom_index.i];
			if (serial[i] < other)
				result.bonds.emplace(serial[i], other);
		}
	}
	for (const auto& pair : ligand.pairs) {
		const auto a = serial[pair.a], b = serial[pair.b];
		result.internal.emplace(std::min(a, b), std::max(a, b), pair.type_pair_index);
	}
	return result;
}

bool same(const signature& first, const signature& second) {
	return first.bonds == second.bonds && first.internal == second.internal && first.xs == second.xs;
}

void require(bool condition, const char* message) {
	if (!condition)
		throw std::runtime_error(message);
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 4)
		throw std::runtime_error("Pass prepared 9FA3, saved bare pose 9FA3/43/model2 and the complete serial bond table");
	model prepared = parse_ligand_pdbqt_from_file(argv[1], atom_type::XS);
	const auto old = topology(prepared);
	const model legacy_reloaded = parse_ligand_pdbqt_from_file(argv[2], atom_type::XS);
	const auto reloaded = topology(legacy_reloaded);
	require(old.internal.size() == 41 && reloaded.internal.size() == 41,
	        "9FA3 regression inputs no longer exercise equal pair counts");
	require(old.bonds.count({13, 14}) == 1 && old.bonds.count({13, 15}) == 0 &&
	        reloaded.bonds.count({13, 14}) == 0 && reloaded.bonds.count({13, 15}) == 1,
	        "9FA3 regression inputs no longer exercise the cyclopropyl edge swap");
	std::ifstream pose_file(argv[2]);
	const std::string pose((std::istreambuf_iterator<char>(pose_file)), std::istreambuf_iterator<char>());
	replace_fixed_pose_coordinates(prepared, pose);
	require(same(old, topology(prepared)), "Coordinate replacement changed inferred scoring topology");
	VINA_FOR(i, prepared.num_atoms())
		require(vec_distance_sqr(prepared.get_coords(i), legacy_reloaded.get_coords(i)) == 0,
		        "A saved atom coordinate was not transferred accurately");

	model corrected = parse_ligand_pdbqt_with_bonds(argv[1], atom_type::XS, argv[3]);
	const auto explicit_graph = topology(corrected);
	require(explicit_graph.bonds.count({13, 14}) == 1 && explicit_graph.bonds.count({13, 15}) == 1 &&
	        explicit_graph.bonds.count({14, 15}) == 1,
	        "Explicit graph did not preserve all cyclopropyl chemical edges");
	const auto graph_after_reload = topology(parse_ligand_pdbqt_with_bonds(argv[2], atom_type::XS, argv[3]));
	require(same(explicit_graph, graph_after_reload), "Explicit topology changed after coordinate rounding");
	replace_fixed_pose_coordinates(corrected, pose);
	require(same(explicit_graph, topology(corrected)) && corrected.num_internal_pairs() == 40,
	        "Transferred coordinates changed the explicit graph or its chemically complete internal pairs");
	std::cout << "Coordinate transfer preserved 41 inferred pairs; complete chemical graph retained all "
	          << "cyclopropyl edges and its 40 internal pairs across both parsing paths.\n";
}
