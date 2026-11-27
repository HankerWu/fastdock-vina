/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
		   The Olson Lab, 
		   The Scripps Research Institute

*/

#include "ad4cache.h"


namespace fs = boost::filesystem;

std::string get_adtype_str(sz& t) {
	switch(t) {
		case AD_TYPE_C : return "C";
		case AD_TYPE_A : return "A";
		case AD_TYPE_N : return "N";
		case AD_TYPE_O : return "O";
		case AD_TYPE_P : return "P";
		case AD_TYPE_S : return "S";
		case AD_TYPE_H : return "H";
		case AD_TYPE_F : return "F";
		case AD_TYPE_I : return "I";
		case AD_TYPE_NA: return "NA";
		case AD_TYPE_OA: return "OA";
		case AD_TYPE_SA: return "SA";
		case AD_TYPE_HD: return "HD";
		case AD_TYPE_Mg: return "Mg";
		case AD_TYPE_Mn: return "Mn";
		case AD_TYPE_Zn: return "Zn";
		case AD_TYPE_Ca: return "Ca";
		case AD_TYPE_Fe: return "Fe";
		case AD_TYPE_Cl: return "Cl";
		case AD_TYPE_Br: return "Br";
		case AD_TYPE_Si: return "Si";
		case AD_TYPE_At: return "At";
		case AD_TYPE_W : return "W";
		default: VINA_CHECK(false);
	}
}

fl ad4cache::eval_atom(const atom& a, const vec& location, fl cap, vec* gradient) const {
	sz type = a.get(atom_type::AD);
	switch (type) {
		case AD_TYPE_G0: case AD_TYPE_G1: case AD_TYPE_G2: case AD_TYPE_G3:
			if (gradient) gradient->assign(0);
			return 0;
		case AD_TYPE_CG0: case AD_TYPE_CG1: case AD_TYPE_CG2: case AD_TYPE_CG3:
			type = AD_TYPE_C;
			break;
	}
	return m_grids[type].evaluate_ad4(location, m_slope, cap, a.charge,
	                                  m_grids[AD_TYPE_SIZE], m_grids[AD_TYPE_SIZE + 1], gradient);
}

fl ad4cache::eval(const model& m, fl cap) const {
	fl energy = 0;
	VINA_FOR(i, m.num_movable_atoms())
		if (m.is_atom_in_ligand(i))
			energy += eval_atom(m.atoms[i], m.coords[i], cap, nullptr);
	return energy;
}

fl ad4cache::eval_intra(model& m, fl cap) const {
	fl energy = 0;
	VINA_FOR(i, m.num_movable_atoms())
		if (!m.is_atom_in_ligand(i))
			energy += eval_atom(m.atoms[i], m.coords[i], cap, nullptr);
	return energy;
}

fl ad4cache::eval_deriv(model& m, fl cap) const {
	fl energy = 0;
	VINA_FOR(i, m.num_movable_atoms())
		energy += eval_atom(m.atoms[i], m.coords[i], cap, &m.minus_forces[i]);
	return energy;
}

bool ad4cache::is_in_grid(const model& m, fl margin) const {
	return m.within_grid(m_gd, margin);
}

bool ad4cache::are_atom_types_grid_initialized(szv atom_types) const {
	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (!is_atom_type_grid_initialized(t)) {
			std::cerr << "ERROR: Affinity map for atom type " << get_adtype_str(t) << " is not present.\n";
			return false;
		}
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE)) {
		std::cerr << "ERROR: Electrostatic map is not present.\n";
		return false;
	}

	if (!is_atom_type_grid_initialized(AD_TYPE_SIZE + 1)) {
		std::cerr << "ERROR: Desolvation map is not present.\n";
		return false;
	}

	return true;
}

void ad4cache::read(const std::string& map_prefix) {
	// AD4 combines atom, electrostatic and desolvation fields. They must
	// describe the same lattice, and failed loads must not alter the cache.
	std::vector<grid> loaded(m_grids.size());
	const path electrostatic(map_prefix + ".e.map");
	const path desolvation(map_prefix + ".d.map");
	const auto reference = loaded[AD_TYPE_SIZE].read_map(electrostatic);
	require_matching_map_dimensions(reference,
	                                loaded[AD_TYPE_SIZE + 1].read_map(desolvation),
	                                desolvation.string());
	VINA_FOR(t, AD_TYPE_SIZE) {
		switch (t) {
			case AD_TYPE_G0: case AD_TYPE_G1: case AD_TYPE_G2: case AD_TYPE_G3:
			case AD_TYPE_CG0: case AD_TYPE_CG1: case AD_TYPE_CG2: case AD_TYPE_CG3:
				continue; // Closure carbons use the ordinary carbon map.
		}
		const path filename(map_prefix + "." + get_adtype_str(t) + ".map");
		if (!fs::exists(filename)) continue;
		const auto dimensions = loaded[t].read_map(filename);
		require_matching_map_dimensions(reference, dimensions, filename.string());
	}
	m_grids.swap(loaded);
	m_gd = reference;
}

void ad4cache::write(const std::string& out_prefix, const szv& atom_types, const std::string& gpf_filename,
					      const std::string& fld_filename, const std::string& receptor_filename) {
	std::string atom_type;
	std::string filename;
	bool got_C_already = false;

	VINA_FOR_IN(i, atom_types) {
		sz t = atom_types[i];

		switch (t)
		{
			case AD_TYPE_G0:
			case AD_TYPE_G1:
			case AD_TYPE_G2:
			case AD_TYPE_G3:
				continue;
			case AD_TYPE_CG0:
			case AD_TYPE_CG1:
			case AD_TYPE_CG2:
			case AD_TYPE_CG3:
				t = AD_TYPE_C;
				break;
		}

		if (t == AD_TYPE_C && got_C_already)
			continue;
		if (t == AD_TYPE_C)
			got_C_already = true;

		if (m_grids[t].initialized()) {
			if (t < AD_TYPE_SIZE)
				atom_type = get_adtype_str(t);
			else if (t == AD_TYPE_SIZE)
				atom_type = "e";
			else if (t == AD_TYPE_SIZE + 1)
				atom_type = "d";

			filename = out_prefix + "." + atom_type + ".map";

			path p(filename);
			ofile out(p);

			// write header
			out << "GRID_PARAMETER_FILE " << gpf_filename << "\n";
			out << "GRID_DATA_FILE " << fld_filename << "\n";
			out << "MACROMOLECULE " << receptor_filename << "\n";

			// m_factor_inv is spacing
			// check that it's the same in every dimension (it must be)
			// check that == operator is OK
			if ((m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[1]) & (m_grids[t].m_factor_inv[0] != m_grids[t].m_factor_inv[2])) {
				printf("m_factor_inv x=%f, y=%f, z=%f\n", m_grids[t].m_factor_inv[0], m_grids[t].m_factor_inv[1], m_grids[t].m_factor_inv[2]);
				return;
			}

			out << "SPACING " << m_grids[t].m_factor_inv[0] << "\n";

			// The number of elements in the grid is an odd number. But NELEMENTS has to be an even number.
			int size_x = (m_grids[t].m_data.dim0() % 2 == 0) ? m_grids[t].m_data.dim0() : m_grids[t].m_data.dim0() - 1;
			int size_y = (m_grids[t].m_data.dim1() % 2 == 0) ? m_grids[t].m_data.dim1() : m_grids[t].m_data.dim1() - 1;
			int size_z = (m_grids[t].m_data.dim2() % 2 == 0) ? m_grids[t].m_data.dim2() : m_grids[t].m_data.dim2() - 1;
			out << "NELEMENTS " << size_x << " " << size_y  << " " << size_z << "\n";

			// center
			fl cx = m_grids[t].m_init[0] + m_grids[t].m_range[0] * 0.5;
			fl cy = m_grids[t].m_init[1] + m_grids[t].m_range[1] * 0.5;
			fl cz = m_grids[t].m_init[2] + m_grids[t].m_range[2] * 0.5;
			out << "CENTER " << cx << " " << cy << " " << cz << "\n";

			// write data
			VINA_FOR(z, m_grids[t].m_data.dim2()) {
				VINA_FOR(y, m_grids[t].m_data.dim1()) {
					VINA_FOR(x, m_grids[t].m_data.dim0()) {
						out << std::setprecision(4) << m_grids[t].m_data(x, y, z) << "\n"; // slow?
					} // x
				} // y
			} // z
		} // map initialized
	} // map atom type
} // cache::write
