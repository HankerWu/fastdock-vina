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

/* Modified for FastDock Vina: reuse type-pair table pointers during grid
   population, share the atom-type mapping and reuse receptor candidates across grid rows.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#include <bitset>

#include "cache.h"
#include "model.h"
#include "precalculate.h"


namespace fs = boost::filesystem;

std::string convert_XS_to_string(sz t) {
	switch(t) {
		case XS_TYPE_C_H     : return "C_H";
		case XS_TYPE_C_P     : return "C_P";
		case XS_TYPE_N_P     : return "N_P";
		case XS_TYPE_N_D     : return "N_D";
		case XS_TYPE_N_A     : return "N_A";
		case XS_TYPE_N_DA    : return "N_DA";
		case XS_TYPE_O_P     : return "O_P";
		case XS_TYPE_O_D     : return "O_D";
		case XS_TYPE_O_A     : return "O_A";
		case XS_TYPE_O_DA    : return "O_DA";
		case XS_TYPE_S_P     : return "S_P";
		case XS_TYPE_P_P     : return "P_P";
		case XS_TYPE_F_H     : return "F_H";
		case XS_TYPE_Cl_H    : return "Cl_H";
		case XS_TYPE_Br_H    : return "Br_H";
		case XS_TYPE_I_H     : return "I_H";
		case XS_TYPE_Si      : return "Si";
		case XS_TYPE_At      : return "At";
		case XS_TYPE_Met_D   : return "Met_D";
		case XS_TYPE_W       : return "W";
		default: VINA_CHECK(false);
	}
}

namespace {

// Ring closure pseudoatoms have no receptor map. The two carbon aliases
// use the same maps as ordinary hydrophobic and polar carbons.
sz affinity_map_type(sz type) {
	switch (type) {
		case XS_TYPE_G0:
		case XS_TYPE_G1:
		case XS_TYPE_G2:
		case XS_TYPE_G3:
			return XS_TYPE_SIZE;
		case XS_TYPE_C_H_CG0:
		case XS_TYPE_C_H_CG1:
		case XS_TYPE_C_H_CG2:
		case XS_TYPE_C_H_CG3:
			return XS_TYPE_C_H;
		case XS_TYPE_C_P_CG0:
		case XS_TYPE_C_P_CG1:
		case XS_TYPE_C_P_CG2:
		case XS_TYPE_C_P_CG3:
			return XS_TYPE_C_P;
		default:
			return type;
	}
}

} // namespace

std::bitset<XS_TYPE_SIZE> cache::map_types(const szv& atom_types) {
	std::bitset<XS_TYPE_SIZE> types;
	for (sz type : atom_types) {
		const sz t = affinity_map_type(type);
		if (t < XS_TYPE_SIZE) types.set(t);
	}
	return types;
}

std::bitset<XS_TYPE_SIZE> cache::missing_map_types(const szv& atom_types) const {
	auto types = map_types(atom_types);
	VINA_FOR(t, XS_TYPE_SIZE)
		if (m_grids[t].initialized()) types.reset(t);
	return types;
}

template<bool Cubic>
fl cache::eval_atoms(const model& m, fl cap, vecv* gradient, bool flex_only) const {
	if (flex_only && m.num_flex() == 0) return 0;
	fl energy = 0;
	VINA_FOR(i, m.num_movable_atoms()) {
		if (flex_only && m.is_atom_in_ligand(i)) continue;
		const sz t = affinity_map_type(m.atoms[i].xs);
		if (t >= XS_TYPE_SIZE) {
			if (gradient) (*gradient)[i].assign(0);
			continue;
		}
		const grid& field = i < m_flex_grids.size() && m_flex_grids[i].initialized()
		                  ? m_flex_grids[i] : m_grids[t];
		if (Cubic)
			energy += field.evaluate_cubic(m.coords[i], m_slope, cap, gradient ? &(*gradient)[i] : nullptr);
		else if (gradient)
			energy += field.evaluate(m.coords[i], m_slope, cap, (*gradient)[i]);
		else
			energy += field.evaluate(m.coords[i], m_slope, cap);
	}
	return energy;
}

template fl cache::eval_atoms<true>(const model&, fl, vecv*, bool) const;

fl cache::eval(const model& m, fl v) const { return eval_atoms<false>(m, v, nullptr, false); }
fl cache::eval_intra(model& m, fl v) const { return eval_atoms<false>(m, v, nullptr, true); }
fl cache::eval_deriv(model& m, fl v) const { return eval_atoms<false>(m, v, &m.minus_forces, false); }

bool cache::is_in_grid(const model& m, fl margin) const {
	return m.within_grid(m_gd, margin);
}

bool cache::are_atom_types_grid_initialized(const szv& atom_types) const {
	for (sz type : atom_types) {
		const sz t = affinity_map_type(type);
		if (t >= XS_TYPE_SIZE) continue;
		if (!is_atom_type_grid_initialized(t)) {
			std::cerr << "ERROR: Affinity map for atom type " <<  convert_XS_to_string(t) << " is not present.\n";
			return false;
		}
	}

	return true;
}

std::string cache::flex_definition(const model& receptor) {
	// Use the supplied domain and atom order, independent of the current pose
	// or any ligand appended later. Comments do not define movable atoms.
	std::string definition;
	for (const auto& line : receptor.flex_context) {
		if (line.first.empty() || line.first.compare(0, 6, "REMARK") == 0) continue;
		definition += line.first + '\n';
	}
	return definition;
}

void cache::read(const std::string& map_prefix, const model* receptor) {
	// Load a complete replacement rather than mixing old and new channels.
	std::vector<grid> loaded(m_grids.size());
	grid_dims reference;
	bool have_maps = false;
	VINA_FOR(t, XS_TYPE_SIZE) {
		if (affinity_map_type(t) != t) continue;
		const path filename(map_prefix + "." + convert_XS_to_string(t) + ".map");
		if (!fs::exists(filename)) continue;
		const auto dimensions = loaded[t].read_map(filename);
		if (!have_maps)
			reference = dimensions;
		else
			require_matching_map_dimensions(reference, dimensions, filename.string());
		have_maps = true;
	}
	if (!have_maps) {
		std::cerr << "\nERROR: No *.map files with prefix \"" << map_prefix << "\"\n";
		exit(EXIT_FAILURE);
	}
	std::vector<grid> flex_grids;
	std::string definition;
	const path manifest(map_prefix + ".flex-fields");
	if (receptor && receptor->num_flex() && fs::exists(manifest)) {
		ifile input(manifest);
		std::string version;
		std::getline(input, version);
		sz count;
		const sz movable = receptor->num_movable_atoms();
		if (version != "FastDock Vina flex fields 1" || !(input >> count) || count > movable)
			throw std::runtime_error("Invalid flex field manifest: " + manifest.string());
		szv indices;
		for (sz n = 0; n < count; ++n) {
			sz i;
			if (!(input >> i) || i >= movable || (!indices.empty() && i <= indices.back()))
				throw std::runtime_error("Invalid flex atom index: " + manifest.string());
			indices.push_back(i);
		}
		input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		definition.assign(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
		if (definition != flex_definition(*receptor))
			throw std::runtime_error("Flexible receptor does not match " + manifest.string());
		for (sz i : indices) {
			const path field(map_prefix + ".flex-" + std::to_string(i) + ".map");
			flex_grids.resize(i + 1);
			require_matching_map_dimensions(reference, flex_grids[i].read_map(field), field.string());
		}
	} else if (receptor && receptor->num_flex()) {
		std::cerr << "Using legacy shared affinity maps for flexible residues: the maps contain no covalent exclusions.\n";
	}
	m_grids.swap(loaded);
	m_flex_grids.swap(flex_grids);
	m_flex_definition.swap(definition);
	m_gd = reference;
}

static void write_xs_field(const grid& field, const path& file, const std::string& gpf,
                           const std::string& fld, const std::string& receptor) {
	const sz nx = field.m_data.dim0() - 1;
	const sz ny = field.m_data.dim1() - 1;
	const sz nz = field.m_data.dim2() - 1;
	if (nx % 2 || ny % 2 || nz % 2)
		throw std::runtime_error("Can't write maps. Number of voxels (NELEMENTS) is odd. Use --force_even_voxels.");
	ofile out(file);
	out << "GRID_PARAMETER_FILE " << gpf << "\nGRID_DATA_FILE " << fld
	    << "\nMACROMOLECULE " << receptor << "\nSPACING " << field.m_factor_inv[0]
	    << "\nNELEMENTS " << nx << " " << ny << " " << nz << "\nCENTER";
	VINA_FOR(axis, 3) out << ' ' << field.m_init[axis] + field.m_range[axis] * 0.5;
	out << "\n" << std::setprecision(4);
	VINA_FOR(z, field.m_data.dim2())
	VINA_FOR(y, field.m_data.dim1())
	VINA_FOR(x, field.m_data.dim0())
		out << field.m_data(x,y,z) << '\n';
}

void cache::write(const std::string& prefix, const szv& atom_types, const std::string& gpf,
                  const std::string& fld, const std::string& receptor) {
	const auto types = map_types(atom_types);
	VINA_FOR(t, XS_TYPE_SIZE) {
		if (types[t] && m_grids[t].initialized())
			write_xs_field(m_grids[t], prefix + "." + convert_XS_to_string(t) + ".map", gpf, fld, receptor);
	}
	szv indices;
	VINA_FOR_IN(i, m_flex_grids) {
		if (!m_flex_grids[i].initialized()) continue;
		write_xs_field(m_flex_grids[i], prefix + ".flex-" + std::to_string(i) + ".map", gpf, fld, receptor);
		indices.push_back(i);
	}
	const path manifest(prefix + ".flex-fields");
	if (indices.empty()) {
		// A reused prefix must not apply fields from a previous flexible receptor.
		if (fs::exists(manifest)) fs::remove(manifest);
	} else {
		ofile out(manifest);
		out << "FastDock Vina flex fields 1\n" << indices.size();
		for (sz i : indices) out << ' ' << i;
		out << '\n' << m_flex_definition;
	}
}

void cache::populate(const model& m, const precalculate& p, const szv_grid& receptor_index,
                     const std::bitset<XS_TYPE_SIZE>& map_types_needed, const std::vector<szv>& rigid_exclusions) {
	szv needed;
	VINA_FOR(t, XS_TYPE_SIZE) {
		if (map_types_needed[t] && !m_grids[t].initialized())
			needed.push_back(t);
	}
	if(needed.empty())
		return;
	// Publish these maps only once every new grid has been filled.
	std::vector<grid> new_grids(needed.size());
	for (auto& g : new_grids)
		g.init(m_gd);
	const sz num_maps = needed.size();

	sz nat = num_atom_types(atom_type::XS);

	const grid& g = new_grids.front();

	const fl cutoff_sqr = p.cutoff_sqr();

	// A grid point uses one distance for all requested ligand types.
	// Keep these energies adjacent for a contiguous read and accumulation.
	const sz sample_count = p.fast_index(cutoff_sqr) + 1;
	std::vector<flv> tables(nat);
	VINA_FOR(t1, nat) {
		if (!receptor_index.atom_types()[t1]) continue;
		flv& samples = tables[t1];
		samples.resize(sample_count * num_maps);
		VINA_FOR_IN(j, needed) {
			const auto& pair = p.element(triangular_matrix_index_permissive(nat, t1, needed[j]));
			VINA_FOR(i, sample_count)
				samples[i * num_maps + j] = pair.fast_sample(i);
		}
	}

	struct x_span {
		sz cell;
		sz end;
	};
	flv x_coords(g.m_data.dim0());
	std::vector<x_span> spans;
	VINA_FOR_IN(x, x_coords) {
		const fl coord = g.index_to_argument(x, 0, 0)[0];
		const sz cell = receptor_index.cell_index(0, coord);
		if (spans.empty() || spans.back().cell != cell)
			spans.push_back({cell, x + 1});
		else
			spans.back().end = x + 1;
		x_coords[x] = coord;
	}

	const auto& receptor = m.grid_atoms();
	// Project the same full-row split into each spatial cell's x span.
	szv atom_x_split(receptor.size());
	VINA_FOR_IN(i, receptor) {
		const auto split = std::lower_bound(x_coords.begin(), x_coords.end(), receptor[i].coords[0]);
		atom_x_split[i] = split - x_coords.begin();
	}

	flv affinities;
	VINA_FOR(z, g.m_data.dim2()) {
		const fl z_coord = g.index_to_argument(0, 0, z)[2];
		const sz z_cell = receptor_index.cell_index(2, z_coord);
		VINA_FOR(y, g.m_data.dim1()) {
			const fl y_coord = g.index_to_argument(0, y, 0)[1];
			const sz y_cell = receptor_index.cell_index(1, y_coord);
			sz start = 0;
			for (const x_span& span : spans) {
				affinities.resize((span.end - start) * num_maps);
				std::fill(affinities.begin(), affinities.end(), 0);
				// Receptor order determines each point's sum. Visit only the
				// in-range run of points on either side of this atom.
				for (sz i : receptor_index.possibilities(span.cell, y_cell, z_cell)) {
					const atom& a = receptor[i];
					const fl dy2 = sqr(a.coords[1] - y_coord);
					const fl dz2 = sqr(a.coords[2] - z_coord);
					if (dy2 + dz2 > cutoff_sqr) continue;
					const sz middle = (std::min)((std::max)(atom_x_split[i], start), span.end);
					const fl* samples = tables[a.xs].data();
					const auto accumulate = [&](sz x) {
						const fl r2 = sqr(a.coords[0] - x_coords[x]) + dy2 + dz2;
						if (!(r2 <= cutoff_sqr)) return false;
						const fl* energies = samples + p.fast_index(r2) * num_maps;
						fl* point = affinities.data() + (x - start) * num_maps;
						VINA_FOR(j, num_maps) point[j] += energies[j];
						return true;
					};
					for (sz x = middle; x < span.end && accumulate(x); ++x) {}
					for (sz x = middle; x > start && accumulate(x - 1); --x) {}
				}
				VINA_FOR_IN(j, needed)
					VINA_RANGE(x, start, span.end)
						new_grids[j].m_data(x, y, z) = affinities[(x - start) * num_maps + j];
				start = span.end;
			}
		}
	}
	// Form each atom's raw field before interpolation, the positive-energy
	// cap and the box penalty. Ligand atoms continue to share the type map.
	m_flex_grids.resize(rigid_exclusions.size());
	if (!rigid_exclusions.empty()) m_flex_definition = flex_definition(m);
	VINA_FOR_IN(i, rigid_exclusions) {
		if (rigid_exclusions[i].empty()) continue;
		const sz t = affinity_map_type(m.atoms[i].xs);
		const auto channel = std::find(needed.begin(), needed.end(), t);
		if (channel == needed.end()) continue;
		grid& field = m_flex_grids[i];
		field = new_grids[channel - needed.begin()];
		for (const sz j : rigid_exclusions[i]) {
			const atom& fixed = receptor[j];
			if (!receptor_index.atom_types()[fixed.xs]) continue;
			const auto& pair = p.element(triangular_matrix_index_permissive(nat, fixed.xs, t));
			const fl cutoff = std::sqrt(cutoff_sqr);
			sz first[3], end[3];
			VINA_FOR(axis, 3) {
				const fl offset = fixed.coords[axis] - field.m_init[axis];
				const fl step = field.m_factor_inv[axis];
				const fl count = field.m_data.dim(axis);
				first[axis] = static_cast<sz>(std::max(0.0, std::min(count, std::floor((offset - cutoff) / step))));
				end[axis] = static_cast<sz>(std::max(0.0, std::min(count, std::ceil((offset + cutoff) / step) + 1)));
			}
			VINA_RANGE(z, first[2], end[2])
			VINA_RANGE(y, first[1], end[1])
			VINA_RANGE(x, first[0], end[0]) {
				const vec point = field.index_to_argument(x, y, z);
				const fl r2 = vec_distance_sqr(point, fixed.coords);
				if (r2 <= cutoff_sqr)
					field.m_data(x, y, z) -= pair.fast_sample(p.fast_index(r2));
			}
		}
	}
	VINA_FOR_IN(j, needed)
		m_grids[needed[j]] = std::move(new_grids[j]);
}
