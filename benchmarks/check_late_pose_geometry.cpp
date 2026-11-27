// Exercise both native and chemical-only metrics on the parsed atom types.
#include "coords.h"
#include "parse_pdbqt.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

void require(bool condition, const char* message) {
	if (!condition) throw std::runtime_error(message);
}

std::string atom_record(sz index, const char* type, fl x) {
	std::ostringstream output;
	output << "ATOM  " << std::right << std::setw(5) << index << "  X   UNL     1    "
	       << std::fixed << std::setprecision(3) << std::setw(8) << x
	       << "   0.000   0.000  1.00  0.00     0.000 " << type << "\n";
	return output.str();
}

model parse_root(const std::string& atoms) {
	return parse_ligand_pdbqt_from_string("ROOT\n" + atoms + "ENDROOT\nTORSDOF 0\n", atom_type::XS);
}

model closure_ligand(fl closure_carbon, fl dummy_offset) {
	std::string atoms = atom_record(1, "C", 0) + atom_record(2, "CG0", closure_carbon);
	const char* types[] = {"G0", "G1", "G2", "G3", "W", "HD"};
	for (sz i = 0; i < 6; ++i)
		atoms += atom_record(i + 3, types[i], 2 + i + dummy_offset);
	return parse_root(atoms);
}

output_type* pose_for(const model& input, fl energy) {
	auto* pose = new output_type(input.get_initial_conf(), energy);
	pose->coords = input.get_chemical_heavy_atom_movable_coords();
	return pose;
}

int main() {
	const model a = closure_ligand(1.5, 0);
	const model shifted_dummies = closure_ligand(1.5, 5);
	const model shifted_carbon = closure_ligand(3.1, 5);

	grid_dims gd;
	gd[0].begin = -2;
	gd[0].end = 2;
	gd[0].n_voxels = 8;
	for (const auto typing : {atom_type::AD, atom_type::XS}) {
		for (const char* type : {"G0", "G1", "G2", "G3", "H", "HD", "CG0", "CG1", "CG2", "CG3", "W", "C"}) {
			const model m = parse_ligand_pdbqt_from_string(
				"ROOT\n" + atom_record(1, "C", 0) + atom_record(2, type, 3) + "ENDROOT\nTORSDOF 0\n", typing);
			const std::string name(type);
			const bool bounded = name[0] != 'G' && (typing == atom_type::AD ||
			                     (name != "H" && name != "HD" && name != "W"));
			require(m.within_grid(gd, 0) == !bounded, "Incorrect force-field grid-domain atom selection");
			require(m.get_grid_atom_movable_coords().size() == 1 + sz(bounded),
			        "Printable box and validity checks must select the same atoms");
		}
	}
	const auto native = a.get_heavy_atom_movable_coords();
	const auto chemical = a.get_chemical_heavy_atom_movable_coords();
	require(native.size() == 3, "Native inter-pose geometry must retain C, CG0 and W only");
	require(chemical.size() == 2 && chemical[0][0] == 0 && chemical[1][0] == 1.5,
	        "The fixed chemical index list must retain C and CG0, omit G0-G3/W/H, and preserve order");
	require(rmsd_upper_bound(native, shifted_dummies.get_heavy_atom_movable_coords()) > 1,
	        "Native hydrated ensemble RMSD should see the displaced W");

	for (const auto typing : {atom_type::XS, atom_type::AD}) {
		const auto placed = [&](fl shift, bool only_markers) {
			std::string records;
			if (!only_markers) records += atom_record(1, "CG2", 0) + atom_record(2, "W", 1);
			for (sz k = 0; k < 4; ++k)
				records += atom_record(k + 3, ("G" + std::to_string(k)).c_str(), 2 + k + shift);
			return parse_ligand_pdbqt_from_string("ROOT\n" + records + "ENDROOT\nTORSDOF 0\n", typing);
		};
		for (bool only_markers : {false, true}) {
			const auto anchor = placed(0, only_markers), moved = placed(7, only_markers);
			require(anchor.get_heavy_atom_movable_coords().size() == (only_markers ? 0u : 2u),
			        "Native coordinates must exclude all four G labels but retain CG and W");
			require(rmsd_upper_bound(anchor.get_heavy_atom_movable_coords(), moved.get_heavy_atom_movable_coords()) == 0 &&
			        anchor.rmsd_upper_bound(moved) == 0 && anchor.rmsd_lower_bound(moved) == 0 &&
			        anchor.rmsd_ligands_upper_bound(moved) == 0, "Native RMSD counted a closure marker");
		}
	}

	output_container candidates;
	candidates.push_back(pose_for(a, -1));
	candidates.push_back(pose_for(shifted_dummies, -2));
	candidates.push_back(pose_for(shifted_carbon, -3));
	const auto selected = select_distinct_poses(candidates, 1, 3);
	require(selected.size() == 2 && selected[0].e == -3 && selected[1].e == -2,
	        "A changed CG0 must be distinct while dummy-only changes are merged");

	const model ethane = parse_root(atom_record(1, "C", -0.77) + atom_record(2, "C", 0.77));
	const model exchanged = parse_root(atom_record(1, "C", 0.77) + atom_record(2, "C", -0.77));
	output_container equivalents;
	equivalents.push_back(pose_for(ethane, -1));
	equivalents.push_back(pose_for(exchanged, -1));
	const auto selected_equivalents = select_distinct_poses(equivalents, 1, 2);
	require(selected_equivalents.size() == 2 && selected_equivalents[0].coords[0][0] == -0.77,
	        "Equivalent atoms must not be remapped and exact score ties must preserve input order");

	std::cout << "Parsed C/CG0, G0-G3/W, hydrogen and index-paired symmetry regressions passed.\n";
}
