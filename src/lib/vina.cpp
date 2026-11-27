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

/* Modified for FastDock Vina: move precomputed tables and grid objects
   to avoid deep copies and build deferred affinity maps at use.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#include "rigid_window.h"
#include "printable_box.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <limits>

#include "vina.h"
#include "contact_rotors.h"
#include "scoring_function.h"
#include "precalculate.h"
#include "pose_template.h"
#include "anchored_search_grid.h"
#include "amide_search_grid.h"
#include "polar_geometry.h"
#include "net_polar.h"
#include "trajectory_archive.h"


void Vina::cite() {
	const std::string cite_message = "\
#################################################################\n\
# If you used AutoDock Vina in your work, please cite:          #\n\
#                                                               #\n\
# J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli  #\n\
# AutoDock Vina 1.2.0: New Docking Methods, Expanded Force      #\n\
# Field, and Python Bindings, J. Chem. Inf. Model. (2021)       #\n\
# DOI 10.1021/acs.jcim.1c00203                                  #\n\
#                                                               #\n\
# O. Trott, A. J. Olson,                                        #\n\
# AutoDock Vina: improving the speed and accuracy of docking    #\n\
# with a new scoring function, efficient optimization and       #\n\
# multithreading, J. Comp. Chem. (2010)                         #\n\
# DOI 10.1002/jcc.21334                                         #\n\
#                                                               #\n\
# Please see https://github.com/ccsb-scripps/AutoDock-Vina for  #\n\
# more information.                                             #\n\
#################################################################\n";

	std::cout << cite_message << '\n';
}

int Vina::generate_seed(const int seed) {
	// Seed generator, if the global seed (m_seed) was defined to 0
	// it seems that we want to generate random seed, otherwise it means
	// that we want to use a particular seed.
	if (seed == 0) {
		return auto_seed();
	} else {
		return seed;
	}
}

void Vina::set_receptor(const std::string& rigid_name, const std::string& flex_name) {
	// Read the receptor PDBQT file
	/* CONDITIONS:
		- 1. AD4/Vina rigid  NO, flex  NO: FAIL
		- 2. AD4      rigid YES, flex YES: FAIL
		- 3. AD4      rigid YES, flex  NO: FAIL
		- 4. AD4      rigid  NO, flex YES: SUCCESS (need to read maps later)
		- 5. Vina     rigid YES, flex YES: SUCCESS
		- 6. Vina     rigid YES, flex  NO: SUCCESS
		- 7. Vina     rigid  NO, flex YES: SUCCESS (need to read maps later)
	*/
	if (rigid_name.empty() && flex_name.empty() && m_sf_choice == SF_VINA) {
		// CONDITION 1
		std::cerr << "ERROR: No (rigid) receptor or flexible residues were specified. (vina.cpp)\n";
		exit(EXIT_FAILURE);
	} else if (m_sf_choice == SF_AD42 && !rigid_name.empty()) {
		// CONDITIONS 2, 3
		std::cerr << "ERROR: Only flexible residues allowed with the AD4 scoring function. No (rigid) receptor.\n";
		exit(EXIT_FAILURE);
	}

	// CONDITIONS 4, 5, 6, 7 (rigid_name and flex_name are empty strings per default)
	m_receptor = parse_receptor_pdbqt(rigid_name, flex_name, m_scoring_function->get_atom_typing());

	m_model = m_receptor;
	invalidate_poses();
	m_receptor_initialized = true;
	// If we are reading another receptor we should not consider the ligand and the map as initialized anymore
	m_ligand_initialized = false;
	m_map_initialized = false;
	m_use_explicit_receptor = false;
	if (m_deferred_map_population && m_no_refine)
		m_non_cache = non_cache();
	m_deferred_map_population = false;
}

void Vina::reset_ligand_model(atom_type::t atom_typing) {
	if (m_receptor_initialized)
		m_model = m_receptor;
	else
		m_model = model(atom_typing);
}

void Vina::finish_ligand_setup(atom_type::t atom_typing) {
	precalculate_byatom precalculated_byatom(*m_scoring_function, m_model);
	if (m_map_initialized && !m_deferred_map_population) {
		const szv atom_types = m_model.get_movable_atom_types(atom_typing);
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if (!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if (!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}
	if (m_map_initialized && m_use_explicit_receptor)
		m_non_cache.prepare_for_model(m_model);

	invalidate_poses();
	m_precalculated_byatom = std::move(precalculated_byatom);
	m_ligand_initialized = true;
}

void Vina::ensure_vina_maps() {
	if (!m_deferred_map_population) return;
	const auto missing = m_grid.missing_map_types(m_model.get_movable_atom_types(atom_type::XS));
	if (missing.any())
		m_non_cache.populate_maps(m_grid, m_receptor, missing, !m_no_refine);
}

void Vina::set_ligand_from_string(const std::string& ligand_string) {
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand file. Ligand string is empty.\n";
		exit(EXIT_FAILURE);
	}

	set_single_ligand_model(parse_ligand_pdbqt_from_string(ligand_string, m_scoring_function->get_atom_typing()));
}

void Vina::set_single_ligand_model(model ligand) {
	const auto atom_typing = m_scoring_function->get_atom_typing();
	reset_ligand_model(atom_typing);
	m_model.append(ligand);
	finish_ligand_setup(atom_typing);
}

void Vina::set_ligand_from_string(const std::vector<std::string>& ligand_string) {
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand list. Ligands list is empty.\n";
		exit(EXIT_FAILURE);
	}
	const atom_type::t atom_typing = m_scoring_function->get_atom_typing();
	// A later ligand can fail type or topology validation. Keep the previous
	// prepared model and its pair tables together until every input parses.
	model combined = m_receptor_initialized ? m_receptor : model(atom_typing);
	for (const auto& ligand : ligand_string)
		combined.append(parse_ligand_pdbqt_from_string(ligand, atom_typing));
	m_model = std::move(combined);
	finish_ligand_setup(atom_typing);
}

void Vina::set_ligand_from_file(const std::string& ligand_name) {
	set_ligand_from_string(get_file_contents(ligand_name));
}

void Vina::set_ligand_from_file_with_bonds(const std::string& ligand_name, const std::string& bonds_name) {
	if (!m_receptor_initialized || m_receptor.num_flex() != 0 || m_sf_choice == SF_AD42 ||
	    (m_map_initialized && !m_use_explicit_receptor))
		throw vina_runtime_error("ligand_bonds requires one Vina or Vinardo ligand and a rigid receptor");
	set_single_ligand_model(parse_ligand_pdbqt_with_bonds(ligand_name, m_scoring_function->get_atom_typing(), bonds_name));
}

void Vina::set_pose_coordinates_for_scoring(const std::string& pose_name) {
	if (!m_ligand_initialized || !m_receptor_initialized || m_no_refine ||
	    (m_map_initialized && !m_use_explicit_receptor) || m_sf_choice == SF_AD42 ||
	    m_model.num_ligands() != 1 || m_model.num_flex() != 0)
		throw vina_runtime_error("pose_coordinates requires a prepared single ligand and an explicit rigid receptor");
	replace_fixed_pose_coordinates(m_model, get_file_contents(pose_name));
}

void Vina::show_rigid_window_pool(const std::string &path, unsigned samples) const {
	if (!m_ligand_initialized || !m_receptor_initialized || m_no_refine ||
		(m_map_initialized && !m_use_explicit_receptor) || m_sf_choice != SF_VINA || m_model.num_ligands() != 1 ||
		m_model.num_flex() != 0)
		throw vina_runtime_error(
			"rigid_window_pool requires one prepared Vina ligand and an explicit rigid receptor");
	std::ifstream input(path);
	if (!input)
		throw vina_runtime_error("Cannot open rigid_window_pool: " + path);
	model working = m_model;
	const rigid_window::evaluator window(working, *m_scoring_function, samples);
	const fl slope =
		m_scoring_function->conf_independent(working, 1) - m_scoring_function->conf_independent(working, 0);
	const vec cap(1000, 1000, 1000);
	std::string line, block;
	std::set<sz> identifiers;
	sz number = 0;
	std::ostringstream rows;
	rows << std::setprecision(12);
	while (std::getline(input, line)) {
		if (!line.empty() && line.back() == '\r')
			line.pop_back();
		std::istringstream fields(line);
		std::string tag;
		fields >> tag;
		if (tag == "MODEL") {
			if (!block.empty())
				throw vina_runtime_error("Nested MODEL in rigid_window_pool");
			long long id;
			std::string extra;
			if (!(fields >> id) || id <= 0 || (fields >> extra) || !identifiers.insert(static_cast<sz>(id)).second)
				throw vina_runtime_error("rigid_window_pool requires unique positive MODEL identifiers");
			number = static_cast<sz>(id);
			block = "MODEL " + std::to_string(number) + "\n";
		} else if (!block.empty()) {
			if (tag != "ENDMDL") {
				block += line + '\n';
				continue;
			}
			std::string extra;
			if (fields >> extra)
				throw vina_runtime_error("Unexpected fields after ENDMDL in rigid_window_pool");
			block += "ENDMDL\n";
			replace_fixed_pose_coordinates(working, block);
			const auto r = window.evaluate(working, slope);
			rows << number << ',' << r.inter << ',' << working.evali(m_precalculated_byatom, cap) << ',' << slope
				 << ',' << r.correction << ',' << r.half << ',' << r.uniform << ',' << r.effective << '\n';
			block.clear();
		} else if (line.find_first_not_of(" \t") != std::string::npos)
			throw vina_runtime_error("Expected MODEL in rigid_window_pool");
	}
	if (input.bad() || !block.empty() || !number)
		throw vina_runtime_error("Incomplete or empty rigid_window_pool");
	std::cout << "model,inter,intra,slope,window_delta,half_delta,uniform_delta,effective_samples\n" << rows.str();
}

void Vina::set_ligand_from_file(const std::vector<std::string>& ligand_name) {
	std::vector<std::string> ligand_string;

	VINA_RANGE(i, 0, ligand_name.size())
		ligand_string.push_back(get_file_contents(ligand_name[i]));

	set_ligand_from_string(ligand_string);
}

/*
void Vina::set_ligand(OpenBabel::OBMol* mol) {
	// Add OBMol to the model
	OpenBabel::OBConversion conv;
	conv.SetOutFormat("PDBQT");
	set_ligand_from_string(conv.WriteString(mol));
}

void Vina::set_ligand(std::vector<OpenBabel::OBMol*> mol) {
	// Add OBMols to the model
	std::vector<std::string> ligand_string;

	OpenBabel::OBConversion conv;
	conv.SetOutFormat("PDBQT");

	VINA_RANGE(i, 0, ligand_name.size())
		ligand_string.push_back(conv.WriteString(mol[i]));

	set_ligand_from_string(ligand_string);
}
*/

void Vina::set_vina_weights(double weight_gauss1, double weight_gauss2, double weight_repulsion,
							double weight_hydrophobic, double weight_hydrogen, double weight_glue,
							double weight_rot) {
	flv weights;

	if (m_sf_choice == SF_VINA) {
		weights.push_back(weight_gauss1);
		weights.push_back(weight_gauss2);
		weights.push_back(weight_repulsion);
		weights.push_back(weight_hydrophobic);
		weights.push_back(weight_hydrogen);
		weights.push_back(weight_glue);
		weights.push_back(5 * weight_rot / 0.1 - 1);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_vinardo_weights(double weight_gauss1, double weight_repulsion,
							   double weight_hydrophobic, double weight_hydrogen, double weight_glue,
							   double weight_rot) {
	flv weights;

	if (m_sf_choice == SF_VINARDO) {
		weights.push_back(weight_gauss1);
		weights.push_back(weight_repulsion);
		weights.push_back(weight_hydrophobic);
		weights.push_back(weight_hydrogen);
		weights.push_back(weight_glue);
		weights.push_back(5 * weight_rot / 0.1 - 1);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_ad4_weights(double weight_ad4_vdw , double weight_ad4_hb,
						   double weight_ad4_elec, double weight_ad4_dsolv,
						   double weight_glue, double weight_ad4_rot) {
	flv weights;

	if (m_sf_choice == SF_AD42) {
		weights.push_back(weight_ad4_vdw);
		weights.push_back(weight_ad4_hb);
		weights.push_back(weight_ad4_elec);
		weights.push_back(weight_ad4_dsolv);
		weights.push_back(weight_glue);
		weights.push_back(weight_ad4_rot);

		// Store in Vina object
		m_weights = weights;

		// Since we set (different) weights, we automatically initialize the forcefield
		set_forcefield();
	}
}

void Vina::set_forcefield() {
    const fl cutoff = m_scoring_function && m_sf_choice != SF_AD42
        ? m_scoring_function->get_cutoff() : default_xs_interaction_cutoff;
    const bool smooth = m_scoring_function && m_scoring_function->smooth_cutoff();
    m_scoring_function = std::make_shared<ScoringFunction>(m_sf_choice, m_weights, cutoff, smooth);
}

// FastDock: bind the cutoff before any ligand or map lookup tables are created.
void Vina::set_interaction_cutoff(double cutoff, bool smooth_cutoff) {
	if (m_sf_choice == SF_AD42)
		throw vina_runtime_error("interaction_cutoff is supported by Vina and Vinardo");
	if (m_ligand_initialized || m_map_initialized)
		throw vina_runtime_error("set interaction_cutoff before loading ligands or maps");
	m_scoring_function = std::make_shared<ScoringFunction>(m_sf_choice, m_weights, cutoff, smooth_cutoff);
}

void Vina::set_rerank_gauss2_scale(double scale) {
	if (!(scale >= 0 && scale <= 1))
		throw vina_runtime_error("rerank_gauss2_scale must be finite and between 0 and 1");
	if (scale != 1.0 && (m_sf_choice != SF_VINA || m_no_refine))
		throw vina_runtime_error("rerank_gauss2_scale below 1 requires Vina scoring with explicit receptor refinement");
	if (scale != m_rerank_gauss2_scale && m_map_initialized)
		throw vina_runtime_error("set rerank_gauss2_scale before building search maps");
	m_rerank_gauss2_scale = scale;
}

std::vector<double> Vina::grid_dimensions_from_ligand(double buffer_size) {
	std::vector<double> box_dimensions(6, 0);
	std::vector<double> box_center(3, 0);
	std::vector<double> max_distance(3, 0);

	// The center of the ligand will be the center of the box
	box_center = m_model.center();

	// Get the furthest atom coordinates from the center in each dimensions
	VINA_FOR(i, m_model.num_movable_atoms()) {
		const vec& atom_coords = m_model.get_coords(i);

		VINA_FOR_IN(j, atom_coords) {
			double distance = std::fabs(box_center[j] - atom_coords[j]);

			if (max_distance[j] < distance)
				max_distance[j] = distance;
		}
	}

	// Get the final dimensions of the box
	box_dimensions[0] = box_center[0];
	box_dimensions[1] = box_center[1];
	box_dimensions[2] = box_center[2];
	box_dimensions[3] = std::ceil((max_distance[0] + buffer_size) * 2);
	box_dimensions[4] = std::ceil((max_distance[1] + buffer_size) * 2);
	box_dimensions[5] = std::ceil((max_distance[2] + buffer_size) * 2);

	return box_dimensions;
}

void Vina::compute_vina_maps(double center_x, double center_y, double center_z,
                             double size_x, double size_y, double size_z,
                             double granularity, bool force_even_voxels, bool defer_population) {
	build_vina_maps(vec(center_x, center_y, center_z), vec(size_x, size_y, size_z),
	                granularity, force_even_voxels, defer_population);
}

void Vina::initialize_vina_batch_maps(double center_x, double center_y, double center_z,
                                     double size_x, double size_y, double size_z, double granularity) {
	compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z,
	                  granularity, false, true);
}

void Vina::build_vina_maps(const vec& center, const vec& span, fl granularity,
                           bool force_even_voxels, bool defer_population) {
	if (m_sf_choice == SF_AD42) {
		throw vina_runtime_error("Cannot compute Vina affinity maps using the AD4 scoring function.");
	} else if (!m_receptor_initialized) {
		throw vina_runtime_error("Cannot compute Vina or Vinardo affinity maps. The (rigid) receptor was not initialized.");
	} else if (span[0] <= 0 || span[1] <= 0 || span[2] <= 0) {
		throw vina_runtime_error("Grid box dimensions must be greater than 0 Angstrom");
	} else if (span[0] * span[1] * span[2] > 27e3) {
		std::cerr << "WARNING: Search space volume is greater than 27000 Angstrom^3 (See FAQ)\n";
	}

	grid_dims gd;
	const fl slope = 1e6;

	// Grid dimensions
	VINA_FOR_IN(i, gd) {
		gd[i].n_voxels = sz(std::ceil(span[i] / granularity));

		// If odd n_voxels increment by 1
		if (force_even_voxels && (gd[i].n_voxels % 2 == 1))
			// because sample points (npts) == n_voxels + 1
			gd[i].n_voxels += 1;

		fl real_span = granularity * gd[i].n_voxels;
		gd[i].begin = center[i] - real_span / 2;
		gd[i].end = gd[i].begin + real_span;
	}

	const std::string scoring_name = m_sf_choice == SF_VINA ? "Vina" : "Vinardo";
	doing((defer_population ? "Preparing " : "Computing ") + scoring_name + " grid", m_verbosity, 0);

	cache grid(gd, slope);
	non_cache receptor(m_receptor, gd, m_scoring_function, slope);
	if (!m_no_refine && m_ligand_initialized)
		receptor.prepare_for_model(m_model);
	if (!defer_population) {
		const auto map_types = cache::map_types(m_ligand_initialized
			? m_model.get_movable_atom_types(atom_type::XS)
			: m_scoring_function->get_atom_types());
		if (map_types.any())
			receptor.populate_maps(grid, m_receptor, map_types, !m_no_refine);
	}
	done(m_verbosity, 0);

	if (!m_no_refine || defer_population) {
		m_non_cache = std::move(receptor);
	} else {
		m_non_cache = non_cache();
	}
	m_grid = std::move(grid);
	m_mixed_search_grid = boost::none;
	m_deferred_map_population = defer_population;
	if (m_rerank_gauss2_scale != 1.0) {
		flv weights = m_weights;
		weights[1] *= m_rerank_gauss2_scale; // Vina's second Gaussian
		m_rerank_scoring_function = std::make_shared<ScoringFunction>(SF_VINA, weights, m_scoring_function->get_cutoff(),
		                                                             m_scoring_function->smooth_cutoff());
	}

	m_use_explicit_receptor = !m_no_refine;
	invalidate_poses();
	m_map_initialized = true;
}

void Vina::load_maps(std::string maps) {
	if (m_armijo_refinement)
		throw vina_runtime_error("Armijo refinement requires search maps built from a rigid receptor");
	if (m_late_pose_selection)
		throw vina_runtime_error("late pose selection requires search maps built from a rigid receptor");
	if (m_rerank_gauss2_scale != 1.0)
		throw vina_runtime_error("rerank_gauss2_scale below 1 requires search maps built from a rigid receptor");
	if (m_scoring_function->smooth_cutoff())
		throw vina_runtime_error("map files cannot encode smooth_cutoff");
	if (m_sf_choice != SF_AD42 && m_scoring_function->get_cutoff() != default_xs_interaction_cutoff)
		throw vina_runtime_error("map files require the default interaction_cutoff");
	const fl slope = 1e6; // FIXME: too large? used to be 100

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		doing("Reading Vina maps", m_verbosity, 0);
		cache grid(slope);
		grid.read(maps, &m_receptor);
		done(m_verbosity, 0);
		m_grid = std::move(grid);
	} else {
		doing("Reading AD4.2 maps", m_verbosity, 0);
		ad4cache grid(slope);
		grid.read(maps);
		done(m_verbosity, 0);
		m_ad4grid = std::move(grid);
	}

	// Check that all the affinity map are present for ligands/flex residues (if initialized already)
	if (m_ligand_initialized) {
		atom_type::t atom_typing = m_scoring_function->get_atom_typing();
		szv atom_types = m_model.get_movable_atom_types(atom_typing);

		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if(!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if(!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}

	// Imported maps cannot be extended by a previously prepared receptor.
	if (m_deferred_map_population && m_no_refine)
		m_non_cache = non_cache();
	m_deferred_map_population = false;
	m_use_explicit_receptor = false;
	invalidate_poses();
	m_map_initialized = true;
}

void Vina::write_maps(const std::string& map_prefix, const std::string& gpf_filename,
					  const std::string& fld_filename, const std::string& receptor_filename) {
	if (m_scoring_function->smooth_cutoff())
		throw vina_runtime_error("map files cannot encode smooth_cutoff");
	if (m_sf_choice != SF_AD42 && m_scoring_function->get_cutoff() != default_xs_interaction_cutoff)
		throw vina_runtime_error("map files require the default interaction_cutoff");
	if (!m_map_initialized) {
		throw vina_runtime_error("Cannot write affinity maps. Affinity maps were not initialized.");
	}

	szv atom_types;
	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	if (m_ligand_initialized)
		atom_types = m_model.get_movable_atom_types(atom_typing);
	else if (m_deferred_map_population)
		throw vina_runtime_error("Load a ligand before writing deferred affinity maps.");
	else
		atom_types = m_scoring_function->get_atom_types();

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		ensure_vina_maps();
		doing("Writing Vina maps", m_verbosity, 0);
		m_grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
		done(m_verbosity, 0);
	} else {
		// Add electrostatics and desolvation maps
		atom_types.push_back(AD_TYPE_SIZE);
		atom_types.push_back(AD_TYPE_SIZE + 1);
		doing("Writing AD4.2 maps", m_verbosity, 0);
		m_ad4grid.write(map_prefix, atom_types, gpf_filename, fld_filename, receptor_filename);
		done(m_verbosity, 0);
	}
}

static bool uses_native_energy_window(const output_type& pose) {
	return pose.amide_penalty >= 0 || pose.template_overlap >= 0 || pose.polar_radial >= 0 || pose.polar_burial >= 0 || pose.net_polar_score || pose.archive_votes >= 0;
}

fl Vina::output_energy_reference() const {
	if (!uses_native_energy_window(m_poses[0]))
		return m_poses[0].e;
	return std::min_element(m_poses.begin(), m_poses.end())->e;
}

std::vector< std::vector<double> > Vina::get_poses_coordinates(int how_many, double energy_range) {
	std::vector< std::vector<double> > coordinates;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		const auto selected = ranked_output_indices(m_poses, how_many, m_output_min_rmsd, output_energy_reference() + energy_range);
		for (sz i : selected) {
			// Push the current pose to model
			m_model.set(m_poses[i].c);
			coordinates.push_back(m_model.get_ligand_coords());
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);
	} else {
		std::cerr << "WARNING: Could not find any pose coordinaates.\n";
	}

	return coordinates;
}

std::vector< std::vector<double> > Vina::get_poses_energies(int how_many, double energy_range) {
	std::vector< std::vector<double> > energies;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		const auto selected = ranked_output_indices(m_poses, how_many, m_output_min_rmsd, output_energy_reference() + energy_range);
		for (sz i : selected) {
			// Push the current pose to model
			energies.push_back({m_poses[i].e,
								m_poses[i].inter, m_poses[i].intra,
								m_poses[i].conf_independent, m_poses[i].unbound});
		}
	} else {
		std::cerr << "WARNING: Could not find any pose energies.\n";
	}

	return energies;
}

static void write_net_polar_satisfaction(std::ostream& output, const char* label,
                                         const net_polar_terms& terms,
                                         const net_polar_satisfaction_score& score) {
	output << label << ' ' << terms.direct_credit_area << ' ' << terms.additional_water_credit_area
	       << ' ' << score.direct_credit << ' ' << score.additional_water_credit
	       << ' ' << score.dry_delta << ' ' << score.water_delta << '\n';
}

std::string Vina::vina_remarks(output_type &pose, fl lb, fl ub) {
	std::ostringstream remark;

	remark.setf(std::ios::fixed, std::ios::floatfield);
	remark.setf(std::ios::showpoint);

	if (pose.amide_penalty >= 0 || pose.net_polar_score || pose.polar_burial >= 0) {
		remark << "REMARK VINA RESULT: " << std::setw(9) << std::setprecision(3) << pose.e
		       << "  " << std::setw(9) << lb << "  " << std::setw(9) << ub << '\n';
		remark << std::defaultfloat << std::noshowpoint
		       << std::setprecision(std::numeric_limits<fl>::max_digits10);
		remark << "REMARK VINA NATIVE: " << pose.e << ' ' << pose.native_rank << '\n';
	}
	if (pose.amide_penalty >= 0) {
		remark << "REMARK VINA AMIDE RESULT: " << pose.amide_rank_score() << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA AMIDE PENALTY: " << pose.amide_penalty << '\n';
		remark << std::fixed << std::showpoint;
	} else if (pose.net_polar_score) {
		const auto& score = *pose.net_polar_score;
		remark << "REMARK VINA NET POLAR RESULT: " << score.rank_score << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA NET POLAR TERMS: " << score.terms.native_hbond << ' '
		       << score.terms.directed_hbond << ' ' << score.terms.direction_delta << ' '
		       << score.terms.buried_area << ' ' << score.terms.water_area << '\n';
		remark << "REMARK VINA NET POLAR DELTAS:";
		for (const fl delta : score.deltas)
			remark << ' ' << delta;
		remark << '\n' << "REMARK VINA NET POLAR MODE: " << net_polar_name() << '\n';
		remark << "REMARK VINA NET POLAR WATER POINTS: " << m_net_polar_waters.size() << '\n';
		if (score.satisfaction)
			write_net_polar_satisfaction(remark, "REMARK VINA NET POLAR SATISFACTION:",
			                             score.terms, *score.satisfaction);
		remark << std::fixed << std::showpoint;
	} else if (pose.polar_burial >= 0) {
		remark << "REMARK VINA POLAR BURIAL NATIVE RANK: " << pose.native_rank << '\n';
		remark << "REMARK VINA POLAR BURIAL RESULT: " << pose.polar_burial_rank_score
		       << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA POLAR BURIAL TERMS: " << pose.polar_burial << ' '
		       << pose.polar_geometry_defect << ' ' << pose.polar_uncompensated << '\n';
		remark << "REMARK VINA POLAR BURIAL DELTAS:";
		for (const fl delta : pose.polar_burial_deltas)
			remark << ' ' << delta;
		remark << '\n' << "REMARK VINA POLAR BURIAL MODE: " << polar_burial_name() << '\n';
		remark << std::fixed << std::showpoint;
	} else if (pose.archive_votes >= 0) {
		remark << std::defaultfloat << std::noshowpoint
		       << std::setprecision(std::numeric_limits<fl>::max_digits10);
		remark << "REMARK VINA ARCHIVE RESULT: " << pose.e << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA NATIVE: " << pose.e << ' ' << pose.native_rank << '\n';
		remark << "REMARK VINA ARCHIVE SUPPORT: " << pose.archive_votes << ' ' << pose.archive_eligible
		       << ' ' << (m_trajectory_archive == archive_mode::support) << '\n';
		remark << "REMARK VINA ARCHIVE SOURCE: " << (pose.archive_from_accepted ? "accepted" : "historical") << '\n';
		remark << "REMARK VINA ARCHIVE MODE: " << trajectory_archive_name() << '\n';
		remark << std::fixed << std::showpoint;
	} else if (pose.polar_radial >= 0) {
		remark << std::defaultfloat << std::noshowpoint
		       << std::setprecision(std::numeric_limits<fl>::max_digits10);
		remark << "REMARK VINA POLAR RESULT: " << polar_rank_score(pose) << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA NATIVE: " << pose.e << ' ' << pose.native_rank << '\n';
		remark << "REMARK VINA POLAR SCORES: " << pose.polar_distance_score << ' ' << pose.polar_directional_score << '\n';
		remark << "REMARK VINA POLAR TERMS: " << pose.polar_radial << ' ' << pose.polar_angular << '\n';
		remark << "REMARK VINA POLAR MODE: " << polar_geometry_name() << '\n';
		remark << std::fixed << std::showpoint;
	} else if (pose.template_overlap >= 0) {
		remark << std::defaultfloat << std::noshowpoint
		       << std::setprecision(std::numeric_limits<double>::max_digits10);
		remark << "REMARK VINA TEMPLATE RESULT: " << pose.template_rank_score << ' ' << lb << ' ' << ub << '\n';
		remark << "REMARK VINA NATIVE: " << pose.e << ' ' << pose.native_rank << '\n';
		remark << "REMARK VINA TEMPLATE OVERLAP: " << pose.template_overlap << '\n';
		remark << "REMARK VINA TEMPLATE STAGE: " << pose_template_stage_name() << '\n';
		remark << std::fixed << std::showpoint;
	} else {
		remark << (pose.native_rank ? "REMARK VINA PRO RESULT: " : "REMARK VINA RESULT: ")
		       << std::setw(9) << std::setprecision(3) << pose.e
		       << "  " << std::setw(9) << std::setprecision(3) << lb
		       << "  " << std::setw(9) << std::setprecision(3) << ub << '\n';
		if (pose.native_rank) {
			remark << "REMARK VINA NATIVE: " << pose.native_energy << "  " << pose.native_rank << '\n';
			remark << "REMARK VINA PRO GAUSS2 SCALE: " << std::defaultfloat << std::noshowpoint
			       << std::setprecision(std::numeric_limits<double>::max_digits10) << m_rerank_gauss2_scale << '\n';
			remark << std::fixed << std::showpoint;
		}
	}

	remark << "REMARK INTER + INTRA:    " << std::setw(12) << std::setprecision(3) << pose.total << "\n";
	remark << "REMARK INTER:            " << std::setw(12) << std::setprecision(3) << pose.inter << "\n";
	remark << "REMARK INTRA:            " << std::setw(12) << std::setprecision(3) << pose.intra << "\n";
	if (m_sf_choice == SF_AD42)
		remark << "REMARK CONF_INDEPENDENT: " << std::setw(12) << std::setprecision(3) << pose.conf_independent << "\n";
	remark << "REMARK UNBOUND:          " << std::setw(12) << std::setprecision(3) << pose.unbound << "\n";

	return remark.str();
}

namespace {
void validate_pose_output(int how_many, double energy_range) {
	if (how_many < 0)
		throw vina_runtime_error("number of poses written must be greater than zero.");
	if (energy_range < 0)
		throw vina_runtime_error("energy range must be greater than zero.");
}
}

std::string Vina::get_poses(int how_many, double energy_range) {
	validate_pose_output(how_many, energy_range);
	std::ostringstream out;
	write_pose_models(out, how_many, energy_range);
	return out.str();
}

void Vina::write_pose_models(std::ostream& out, int how_many, double energy_range) {
	int n = 0;
	std::string remarks;

	if (!m_poses.empty()) {
		// Get energy from the best conf
		const double best_energy = output_energy_reference();
		const bool range_excludes_first = uses_native_energy_window(m_poses[0]) &&
		                                 m_poses[0].e > best_energy + energy_range;
		boost::optional<model> range_reference;
		boost::optional<contact_rotors> rotor_report;
		if (m_contact_rotor_score)
			rotor_report.emplace(m_model);

		const auto selected = ranked_output_indices(m_poses, how_many, m_output_min_rmsd, best_energy + energy_range);
		for (sz i : selected) {
			// Push the current pose to model
			m_model.set(m_poses[i].c);

			// Write conf
			fl lb = m_poses[i].lb;
			fl ub = m_poses[i].ub;
			if (range_excludes_first) {
				if (!range_reference)
					range_reference = m_model;
				lb = m_model.rmsd_lower_bound(*range_reference);
				ub = m_model.rmsd_upper_bound(*range_reference);
			}
			remarks = vina_remarks(m_poses[i], lb, ub);
			if (rotor_report)
				remarks += rotor_report->remark(m_model, m_poses[i].total - m_poses[i].unbound);
			// Write into the final buffer without a separate string per model.
			const bool clean_score_remarks = m_poses[i].native_rank || m_late_pose_selection || m_armijo_refinement || m_contact_rotor_score;
			m_model.write_model(out, n + 1, remarks, clean_score_remarks ? is_vina_score_remark : nullptr);

			n++;
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);

	} else {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
	}

}

void Vina::write_poses(const std::string& output_name, int how_many, double energy_range) {
	if (!m_search_completed && m_poses.empty()) {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
		return;
	}
	// A completed search can legitimately have no in-domain output. Replace
	// a previous run's file with an empty result rather than leaving stale poses.
	validate_pose_output(how_many, energy_range);
	ofile out(make_path(output_name));
	write_pose_models(out, how_many, energy_range);
}

void Vina::write_pose(const std::string& output_name, const std::string& remark) {
	std::ostringstream format_remark;
	format_remark.setf(std::ios::fixed, std::ios::floatfield);
	format_remark.setf(std::ios::showpoint);

	// Add REMARK keyword to be PDB valid
	if(!remark.empty()){
		format_remark << "REMARK " << remark << " \n";
	}

	ofile f(make_path(output_name));
	m_model.write_structure(f, format_remark.str());
}

void Vina::randomize(const int max_steps) {
	// Randomize ligand/flex residues conformation
	// Check the box was defined
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do ligand randomization. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do ligand randomization. Affinity maps were not initialized.");
	}

	conf c;
	int seed = generate_seed();
	double penalty = 0;
	double best_clash_penalty = 0;
	std::stringstream sstm;
	rng generator(static_cast<rng::result_type>(seed));

	// It's okay to take the initial conf since we will randomize it
	conf init_conf = m_model.get_initial_conf();
	conf best_conf = init_conf;

	sstm << "Randomize conformation (random seed: " << seed << ")";
	doing(sstm.str(), m_verbosity, 0);
	VINA_FOR(i, max_steps) {
		c = init_conf;
		c.randomize(m_grid.corner1(), m_grid.corner2(), generator);
		m_model.set(c);
		penalty = m_model.clash_penalty();

		if (i == 0 || penalty < best_clash_penalty) {
			best_conf = c;
			best_clash_penalty = penalty;
		}
	}
	done(m_verbosity, 0);

	m_model.set(best_conf);

	if (m_verbosity > 1) {
		std::cout << "Clash penalty: " << best_clash_penalty << "\n";
	}
}

void Vina::show_score(const std::vector<double> energies) {
	if (m_contact_rotor_score) {
		double raw = -energies.at(7);
		for (sz i = 1; i <= 5; ++i) raw += energies[i];
		std::cout << contact_rotors(m_model).remark(m_model, raw);
	}
	std::cout << "Estimated Free Energy of Binding   : " << std::fixed << std::setprecision(3) << energies[0] << " (kcal/mol) [=(1)+(2)+(3)-(4)]\n";
	std::cout << "(1) Final Intermolecular Energy    : " << std::fixed << std::setprecision(3) << energies[1] + energies[2] << " (kcal/mol)\n";
	std::cout << "    Ligand - Receptor              : " << std::fixed << std::setprecision(3) << energies[1] << " (kcal/mol)\n";
	std::cout << "    Ligand - Flex side chains      : " << std::fixed << std::setprecision(3) << energies[2] << " (kcal/mol)\n";
	std::cout << "(2) Final Total Internal Energy    : " << std::fixed << std::setprecision(3) << energies[3] + energies[4] + energies[5] << " (kcal/mol)\n";
	std::cout << "    Ligand                         : " << std::fixed << std::setprecision(3) << energies[5] << " (kcal/mol)\n";
	std::cout << "    Flex   - Receptor              : " << std::fixed << std::setprecision(3) << energies[3] << " (kcal/mol)\n";
	std::cout << "    Flex   - Flex side chains      : " << std::fixed << std::setprecision(3) << energies[4] << " (kcal/mol)\n";
	std::cout << "(3) Torsional Free Energy          : " << std::fixed << std::setprecision(3) << energies[6] << " (kcal/mol)\n";
	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		std::cout << "(4) Unbound System's Energy        : " << std::fixed << std::setprecision(3) << energies[7] << " (kcal/mol)\n";
	} else {
		std::cout << "(4) Unbound System's Energy [=(2)] : " << std::fixed << std::setprecision(3) << energies[7] << " (kcal/mol)\n";
	}
}

const igrid& Vina::vina_receptor_interactions() {
	if (m_use_explicit_receptor)
		return m_non_cache;
	ensure_vina_maps();
	return m_grid;
}

std::vector<double> Vina::score(double intramolecular_energy) {
	// Score the current conf in the model
	double total = 0;
	double inter = 0;
	double intra = 0;
	double all_grids = 0; // ligand & flex
	double lig_grids = 0;
	double flex_grids = 0;
	double lig_intra = 0;
	double conf_independent = 0;
	double inter_pairs = 0;
	double intra_pairs = 0;
	const vec authentic_v(1000, 1000, 1000);
	std::vector<double> energies;

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		const igrid& receptor = vina_receptor_interactions();
		// Inter
		all_grids = receptor.eval(m_model, authentic_v[1]); // [1] ligand & flex -- receptor
		inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
		// Intra
		flex_grids = receptor.eval_intra(m_model, authentic_v[1]); // [1] flex -- receptor
		intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] flex_i -- flex_i and flex_i -- flex_j
		lig_grids = all_grids - flex_grids;
		inter = lig_grids + inter_pairs;
		lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2] ligand_i -- ligand_i
		intra = flex_grids + intra_pairs + lig_intra;
		// Total
		total = m_scoring_function->conf_independent(m_model, inter + intra - intramolecular_energy); // we pass intermolecular energy from the best pose
		// Torsion, we want to know how much torsion penalty was added to the total energy
		conf_independent = total - (inter + intra - intramolecular_energy);
	} else {
		// Inter
		lig_grids = m_ad4grid.eval(m_model, authentic_v[1]); // [1] ligand -- grid
		inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
		inter = lig_grids + inter_pairs;
		// Intra
		flex_grids = m_ad4grid.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
		intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // [1] flex_i -- flex_i and flex_i -- flex_j
		lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // [2] ligand_i -- ligand_i
		intra = flex_grids + intra_pairs + lig_intra;
		// Torsion
		conf_independent = m_scoring_function->conf_independent(m_model, 0); // [3] we can pass e=0 because we do not modify the energy like in vina
		// Total
		total = inter + conf_independent; // (+ intra - intra)
	}

	energies.push_back(total);
	energies.push_back(lig_grids);
	energies.push_back(inter_pairs);
	energies.push_back(flex_grids);
	energies.push_back(intra_pairs);
	energies.push_back(lig_intra);
	energies.push_back(conf_independent);

	if (m_sf_choice == SF_VINA  || m_sf_choice == SF_VINARDO) {
		energies.push_back(intramolecular_energy);
	} else {
		energies.push_back(intra);
	}

	return energies;
}

bool Vina::current_pose_in_grid() const {
	return m_sf_choice == SF_AD42 ? m_ad4grid.is_in_grid(m_model) : m_grid.is_in_grid(m_model);
}

std::vector<double> Vina::score() {
	// Score the current conf in the model
	// Check if ff and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot score the pose. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot score the pose. Affinity maps were not initialized.");
	} else if (!current_pose_in_grid()) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);

	if(m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, vina_receptor_interactions(), authentic_v);
	}

	std::vector<double> energies = score(intramolecular_energy);
	return energies;
}

std::vector<double> Vina::optimize(output_type& out, int max_steps, int max_evals, bool rigid_metric) {
	// Local optimization of the ligand conf
	const vec authentic_v(1000, 1000, 1000);
	std::vector<double> energies_before_opt;
	std::vector<double> energies_after_opt;
	int evalcount = 0;

	// Define the number minimization steps based on the number moving atoms
	if (max_steps == 0) {
		max_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
		if (m_verbosity > 1)
			std::cout << "Number of local optimization steps: " << max_steps << "\n";
	}
	quasi_newton quasi_newton_par(m_model.get_size(), max_steps);
	quasi_newton_par.strict_armijo = m_armijo_refinement;
	quasi_newton_par.geometry_metric = m_geometry_refinement;
	quasi_newton_par.rigid_metric = rigid_metric;

	if (m_verbosity > 1) {
		std::cout << "Before local optimization:\n";
		energies_before_opt = score();
		show_score(energies_before_opt);
	}

	doing("Performing local search", m_verbosity, 0);
	// Share the gradient budget across all box-recovery attempts.
	for (sz attempt = 0; attempt < 5 && (max_evals == 0 || evalcount < max_evals); ++attempt) {
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			quasi_newton_par(m_model, m_precalculated_byatom, vina_receptor_interactions(), out, authentic_v, evalcount, nullptr, max_evals);
		} else {
			quasi_newton_par(m_model, m_precalculated_byatom, m_ad4grid, out, authentic_v, evalcount);
		}
		if (current_pose_in_grid()) break;
	}
	done(m_verbosity, 0);

	energies_after_opt = score();
	if (m_armijo_refinement) {
		std::ostringstream report;
		const char* optimizer = rigid_metric ? "rigid_metric" : (m_geometry_refinement ? "geometry" : "armijo");
		report << "Local refinement: optimizer=" << optimizer
		       << " gradient_evals=" << evalcount << " objective_kcal_mol="
		       << std::setprecision(std::numeric_limits<fl>::max_digits10) << out.e << '\n';
		std::cout << report.str();
	}

	return energies_after_opt;
}

std::vector<double> Vina::optimize(int max_steps) {
	return optimize(max_steps, 0);
}

std::vector<double> Vina::optimize(int max_steps, int max_evals, bool rigid_metric) {
	// Local optimization of the ligand conf
	// Check if ff, box and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Affinity maps were not initialized.");
	} else if (!current_pose_in_grid()) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

	if (m_armijo_refinement && (!m_use_explicit_receptor || m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("Armijo refinement requires one ligand and a rigid explicit receptor");
	if (max_evals < 0 || (max_evals > 0 && !m_armijo_refinement))
		throw vina_runtime_error("A local gradient budget requires Armijo refinement and non-negative max_evals");

	if (rigid_metric && (!m_armijo_refinement || m_geometry_refinement))
		throw vina_runtime_error("Rigid metric requires Armijo without Cartesian preconditioning");

	double e = 0;
	conf c;

	if (!m_poses.empty()) {
		// if m_poses is not empty, it means that we did a docking before
		// But it is really that useful to minimize after docking?
		e = m_poses[0].e;
		c = m_poses[0].c;
	} else {
		c = m_model.get_initial_conf();
	}

	output_type out(c, e);

	return optimize(out, max_steps, max_evals, rigid_metric);
}

output_container Vina::remove_redundant(const output_container &in, fl min_rmsd) {
	output_container tmp;
	VINA_FOR_IN(i, in)
	add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

void Vina::rerank_poses(output_container& poses, const vec& authentic_v) {
	const precalculate_byatom ligand_scores(*m_rerank_scoring_function, m_model);
	precalculate receptor_scores(m_rerank_scoring_function->get_cutoff());
	m_non_cache.prepare_for_model(m_model, *m_rerank_scoring_function, receptor_scores);
	VINA_FOR_IN(i, poses) {
		output_type& pose = poses[i];
		pose.native_energy = pose.e;
		pose.native_rank = i + 1;
		m_model.set(pose.c);
		pose.inter = m_non_cache.eval(m_model, authentic_v[1], receptor_scores);
		pose.intra = m_model.evali(ligand_scores, authentic_v);
		pose.total = pose.inter + pose.intra;
	}

	// Apply Vina's common unbound-energy convention to the new scoring function:
	// every candidate subtracts the internal energy of the lowest-total-energy pose.
	const fl unbound = std::min_element(poses.begin(), poses.end(), [](const output_type& a, const output_type& b) {
		return a.total < b.total;
	})->intra;
	for (output_type& pose : poses) {
		const fl corrected = pose.total - unbound;
		pose.e = m_rerank_scoring_function->conf_independent(m_model, corrected);
		pose.conf_independent = pose.e - corrected;
		pose.unbound = unbound;
	}
	poses.sort([](const output_type& a, const output_type& b) {
		return a.e == b.e ? a.native_rank < b.native_rank : a.e < b.e;
	});
}

void Vina::set_search_anchor(const std::string& filename, const std::string& trajectories) {
	anchor_trajectories selected;
	if (trajectories == "all") selected = anchor_trajectories::all;
	else if (trajectories == "even") selected = anchor_trajectories::even;
	else if (trajectories == "odd") selected = anchor_trajectories::odd;
	else throw vina_runtime_error("search_anchor_trajectories must be all, even, or odd");
	if (filename.empty() && selected != anchor_trajectories::all)
		throw vina_runtime_error("search_anchor_trajectories requires search_anchor");
	m_search_anchor_path = filename;
	m_anchor_trajectories = selected;
}

void Vina::set_late_pose_selection(bool enabled) {
	if (enabled && (m_sf_choice != SF_VINA || m_no_refine))
		throw vina_runtime_error("late pose selection requires Vina scoring with explicit receptor refinement");
	if (enabled != m_late_pose_selection && m_map_initialized)
		throw vina_runtime_error("set late pose selection before initializing search maps");
	m_late_pose_selection = enabled;
}

void Vina::set_armijo_refinement(bool enabled, bool geometry_metric) {
	if (geometry_metric && !enabled)
		throw vina_runtime_error("Geometry preconditioning requires Armijo refinement");
	if (enabled && (m_sf_choice != SF_VINA || m_no_refine))
		throw vina_runtime_error("Armijo refinement requires Vina scoring with explicit receptor refinement");
	if ((enabled != m_armijo_refinement || geometry_metric != m_geometry_refinement) && m_map_initialized)
		throw vina_runtime_error("set Armijo refinement before initializing search maps");
	m_armijo_refinement = enabled;
	m_geometry_refinement = geometry_metric;
}

void Vina::set_rigid_reorientation(bool mixed, bool centered) {
	if ((mixed || centered) && (m_sf_choice != SF_VINA || m_no_refine))
		throw vina_runtime_error("rigid_reorientation requires Vina scoring with explicit receptor refinement");
	m_mixed_reorientation = mixed;
	m_centered_rotation = centered;
}

void Vina::set_separate_hunt_chain(bool enabled) {
	if (enabled && (m_sf_choice != SF_VINA || m_no_refine))
		throw vina_runtime_error("separate_hunt_chain requires Vina scoring with explicit receptor refinement");
	m_separate_hunt_chain = enabled;
}

void Vina::set_pocket_start_radius(double radius) {
	if (!std::isfinite(radius) || radius <= 0)
		throw vina_runtime_error("pocket_start_radius must be finite and positive");
	if (m_sf_choice == SF_AD42)
		throw vina_runtime_error("pocket_start_radius requires Vina or Vinardo scoring");
	m_pocket_start_radius = radius;
}

void Vina::set_pose_template(const std::string& filename, const std::string& stage) {
	if (m_map_initialized)
		throw vina_runtime_error("set pose_template before initializing search maps");
	if (filename.empty())
		throw vina_runtime_error("pose_template must name a heavy-atom coordinate file");
	if (stage == "rank")
		m_template_stage = template_stage::rank;
	else if (stage == "contact_rank")
		m_template_stage = template_stage::contact_rank;
	else if (stage == "select")
		m_template_stage = template_stage::select;
	else if (stage == "refine")
		m_template_stage = template_stage::refine;
	else
		throw vina_runtime_error("template_stage must be rank, select or refine for element overlap, or contact_rank for receptor contacts");
	m_pose_template_path = filename;
}

void Vina::set_polar_geometry(const std::string& mode) {
	if (m_map_initialized)
		throw vina_runtime_error("set polar_geometry before initializing search maps");
	if (mode == "distance")
		m_polar_geometry = polar_mode::distance;
	else if (mode == "directional")
		m_polar_geometry = polar_mode::directional;
	else if (mode == "penalty")
		m_polar_geometry = polar_mode::penalty;
	else
		throw vina_runtime_error("polar_geometry must be distance, directional or penalty");
}

void Vina::set_polar_burial(const std::string& mode) {
	if (m_map_initialized)
		throw vina_runtime_error("set polar_burial before initializing search maps");
	if (mode == "burial")
		m_polar_burial = burial_mode::burial;
	else if (mode == "geometry")
		m_polar_burial = burial_mode::geometry;
	else if (mode == "compensated")
		m_polar_burial = burial_mode::compensated;
	else
		throw vina_runtime_error("polar_burial must be burial, geometry or compensated");
}

const char* Vina::polar_burial_name() const {
	switch (m_polar_burial) {
		case burial_mode::burial: return "burial";
		case burial_mode::geometry: return "geometry";
		case burial_mode::compensated: return "compensated";
		case burial_mode::off: return "";
	}
	return "";
}

bool Vina::supports_polar_burial() const {
	return m_sf_choice == SF_VINA && !m_no_refine && m_receptor_initialized &&
	       m_model.num_ligands() == 1 && m_model.num_flex() == 0 &&
	       m_rerank_gauss2_scale == 1 && !m_armijo_refinement && m_pose_template_path.empty() &&
	       m_polar_geometry == polar_mode::off && m_trajectory_archive == archive_mode::off;
}

std::array<fl, 3> Vina::polar_burial_deltas(fl burial, fl geometry, fl uncompensated) const {
	// Standard Vina's only conformation-independent transform divides energy
	// by a positive torsional normalization fixed for the ligand.
	constexpr fl strength = 0.6;
	const fl normalization = m_scoring_function->conf_independent(m_model, strength);
	return {{normalization * burial, normalization * geometry, normalization * uncompensated}};
}

non_cache Vina::fixed_pose_receptor() {
	const auto box = grid_dimensions_from_ligand(1);
	grid_dims dimensions;
	VINA_FOR_IN(i, dimensions) {
		dimensions[i].begin = box[i] - box[i + 3] / 2;
		dimensions[i].end = box[i] + box[i + 3] / 2;
		dimensions[i].n_voxels = 1;
	}
	non_cache receptor(m_receptor, dimensions, m_scoring_function, 0);
	receptor.prepare_for_model(m_model);
	return receptor;
}

void Vina::show_polar_burial() {
	if (!supports_polar_burial())
		throw vina_runtime_error("polar_burial score_only requires Vina with one ligand and an explicit rigid receptor");
	const polar_geometry geometry(m_model);
	const auto terms = geometry.evaluate_burial(m_model);
	const auto deltas = polar_burial_deltas(terms.burial, terms.geometry, terms.uncompensated);

	// Reuse the native explicit scorer with an index enclosing the supplied
	// pose. This gives the rigid inter and intra terms without map population
	// or box penalties. Use a zero common unbound reference for diagnostics.
	const non_cache receptor = fixed_pose_receptor();
	const vec authentic_v(1000, 1000, 1000);
	const fl inter = receptor.eval(m_model, authentic_v[1]);
	const fl intra = m_model.evali(m_precalculated_byatom, authentic_v);
	const fl native_zero = m_scoring_function->conf_independent(m_model, inter + intra);

	std::ostringstream report;
	report << std::fixed << std::setprecision(12)
	       << "Polar burial terms: " << terms.burial << ' ' << terms.geometry << ' ' << terms.uncompensated << '\n'
	       << "Polar burial deltas: " << deltas[0] << ' ' << deltas[1] << ' ' << deltas[2] << '\n'
	       << "Polar burial native: " << inter << ' ' << intra << ' ' << native_zero << '\n';
	std::cout << report.str();
}

void Vina::rank_polar_burial_poses(output_container& poses) {
	const polar_geometry geometry(m_model);
	const sz mode = static_cast<sz>(m_polar_burial) - 1;
	VINA_FOR_IN(i, poses) {
		output_type& pose = poses[i];
		m_model.set(pose.c);
		const auto terms = geometry.evaluate_burial(m_model);
		pose.native_rank = i + 1;
		pose.polar_burial = terms.burial;
		pose.polar_geometry_defect = terms.geometry;
		pose.polar_uncompensated = terms.uncompensated;
		pose.polar_burial_deltas = polar_burial_deltas(terms.burial, terms.geometry, terms.uncompensated);
		pose.polar_burial_rank_score = pose.e + pose.polar_burial_deltas[mode];
	}
	poses.sort([](const output_type& a, const output_type& b) {
		return a.polar_burial_rank_score == b.polar_burial_rank_score ? a.native_rank < b.native_rank
		                                                            : a.polar_burial_rank_score < b.polar_burial_rank_score;
	});
}

void Vina::set_net_polar(const std::string& mode, const std::string& waters) {
	if (m_map_initialized)
		throw vina_runtime_error("set net_polar before initializing search maps");
	net_polar_mode selection;
	if (mode == "native")
		selection = net_polar_mode::native;
	else if (mode == "direction")
		selection = net_polar_mode::direction;
	else if (mode == "burial")
		selection = net_polar_mode::burial;
	else if (mode == "combined")
		selection = net_polar_mode::combined;
	else if (mode == "water")
		selection = net_polar_mode::water;
	else if (mode == "satisfied")
		selection = net_polar_mode::satisfied;
	else if (mode == "satisfied_water")
		selection = net_polar_mode::satisfied_water;
	else
		throw vina_runtime_error("net_polar must be native, direction, burial, combined, water, satisfied or satisfied_water");

	vecv water_points;
	if (!waters.empty()) {
		std::ifstream input(waters);
		if (!input)
			throw vina_runtime_error("cannot open net_polar_waters file: " + waters);
		std::string line;
		sz line_number = 0;
		while (std::getline(input, line)) {
			++line_number;
			std::istringstream record(line.substr(0, line.find('#')));
			record >> std::ws;
			if (record.eof())
				continue;
			vec point;
			std::string extra;
			if (!(record >> point[0] >> point[1] >> point[2]) || record >> extra ||
			    !std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2]))
				throw vina_runtime_error("net_polar_waters expects three finite coordinates on line " +
				                         std::to_string(line_number));
			water_points.push_back(point);
		}
		if (input.bad() || water_points.empty())
			throw vina_runtime_error("net_polar_waters requires a readable file with at least one aligned water oxygen");
	}
	if ((selection == net_polar_mode::water || selection == net_polar_mode::satisfied_water) && water_points.empty())
		throw vina_runtime_error("net_polar " + mode + " requires net_polar_waters from independent aligned structural evidence");
	m_net_polar_waters = std::move(water_points);
	m_net_polar = selection;
}

const char* Vina::net_polar_name() const {
	switch (m_net_polar) {
		case net_polar_mode::native: return "native";
		case net_polar_mode::direction: return "direction";
		case net_polar_mode::burial: return "burial";
		case net_polar_mode::combined: return "combined";
		case net_polar_mode::water: return "water";
		case net_polar_mode::satisfied: return "satisfied";
		case net_polar_mode::satisfied_water: return "satisfied_water";
		case net_polar_mode::off: return "";
	}
	return "";
}

bool Vina::uses_net_polar_satisfaction() const {
	return m_net_polar == net_polar_mode::satisfied || m_net_polar == net_polar_mode::satisfied_water;
}

std::array<fl, 4> Vina::net_polar_deltas(const net_polar_terms& terms, fl normalization) const {
	const fl direction = terms.direction_delta * normalization;
	const fl burial = net_polar::area_cost() * terms.buried_area * normalization;
	const fl water_credit = net_polar::area_cost() * terms.water_area * normalization;
	return {{direction, burial, direction + burial, direction + burial - water_credit}};
}

static net_polar_satisfaction_score satisfaction_score(const net_polar_terms& terms, fl normalization,
                                                       fl direction) {
	const fl area_cost = net_polar::area_cost() * normalization;
	const fl dry_area = terms.buried_area - terms.direct_credit_area;
	// Roundoff in the summed credits must not turn a residual cost into a reward.
	const fl wet_area = std::max(fl(0), dry_area - terms.additional_water_credit_area);
	return {area_cost * terms.direct_credit_area, area_cost * terms.additional_water_credit_area,
	        direction + area_cost * dry_area, direction + area_cost * wet_area};
}

void Vina::show_net_polar() {
	if (!supports_polar_burial() || m_polar_burial != burial_mode::off)
		throw vina_runtime_error("net_polar score_only requires Vina with one ligand and an explicit rigid receptor");
	const non_cache receptor = fixed_pose_receptor();
	const vec authentic_v(1000, 1000, 1000);
	const net_polar polar(m_model, m_net_polar_waters);
	net_polar_terms terms;
	const fl inter = receptor.directional_eval(m_model, authentic_v[1], *m_scoring_function, polar, terms);
	const fl intra = m_model.evali(m_precalculated_byatom, authentic_v);
	const fl native_zero = m_scoring_function->conf_independent(m_model, inter + intra);
	const fl normalization = m_scoring_function->conf_independent(m_model, 1);
	const bool satisfaction = uses_net_polar_satisfaction();
	polar.add_surface(m_model, terms, satisfaction);
	const auto deltas = net_polar_deltas(terms, normalization);

	std::ostringstream report;
	report << std::fixed << std::setprecision(12)
	       << "Net polar terms: " << terms.native_hbond << ' ' << terms.directed_hbond << ' '
	       << terms.direction_delta << ' ' << terms.buried_area << ' ' << terms.water_area << '\n'
	       << "Net polar deltas: " << deltas[0] << ' ' << deltas[1] << ' ' << deltas[2] << ' ' << deltas[3] << '\n'
	       << "Net polar native: " << inter << ' ' << intra << ' ' << native_zero << '\n'
	       << "Net polar normalization: " << normalization << '\n'
	       << "Net polar water points: " << m_net_polar_waters.size() << '\n';
	if (satisfaction)
		write_net_polar_satisfaction(report, "Net polar satisfaction:", terms,
		                             satisfaction_score(terms, normalization, deltas[0]));
	std::cout << report.str();
}

void Vina::rank_net_polar_poses(output_container& poses) {
	const net_polar polar(m_model, m_net_polar_waters);
	const fl normalization = m_scoring_function->conf_independent(m_model, 1);
	const bool native = m_net_polar == net_polar_mode::native;
	const bool satisfaction = uses_net_polar_satisfaction();
	VINA_FOR_IN(i, poses) {
		output_type& pose = poses[i];
		m_model.set(pose.c);
		net_polar_terms terms;
		m_non_cache.directional_eval(m_model, 1000, *m_scoring_function, polar, terms);
		polar.add_surface(m_model, terms, satisfaction);
		const auto deltas = net_polar_deltas(terms, normalization);
		pose.native_rank = i + 1;
		net_polar_pose_score score{terms, deltas, pose.e, boost::none};
		if (satisfaction) {
			score.satisfaction = satisfaction_score(terms, normalization, deltas[0]);
			score.rank_score += m_net_polar == net_polar_mode::satisfied ?
			                   score.satisfaction->dry_delta : score.satisfaction->water_delta;
		} else if (!native) {
			const auto mode = static_cast<sz>(m_net_polar) - static_cast<sz>(net_polar_mode::direction);
			score.rank_score += deltas[mode];
		}
		pose.net_polar_score = score;
	}
	if (!native) {
		poses.sort([](const output_type& a, const output_type& b) {
			const fl first = a.net_polar_score->rank_score;
			const fl second = b.net_polar_score->rank_score;
			return first == second ? a.native_rank < b.native_rank : first < second;
		});
	}
}

const char* Vina::polar_geometry_name() const {
	switch (m_polar_geometry) {
		case polar_mode::distance: return "distance";
		case polar_mode::directional: return "directional";
		case polar_mode::penalty: return "penalty";
		case polar_mode::off: return "";
	}
	return "";
}

void Vina::set_trajectory_archive(const std::string& mode) {
	if (m_map_initialized)
		throw vina_runtime_error("set trajectory_archive before initializing search maps");
	if (mode == "accepted")
		m_trajectory_archive = archive_mode::accepted;
	else if (mode == "support")
		m_trajectory_archive = archive_mode::support;
	else
		throw vina_runtime_error("trajectory_archive must be accepted or support");
}

const char* Vina::trajectory_archive_name() const {
	return m_trajectory_archive == archive_mode::support ? "support" : "accepted";
}

fl Vina::polar_rank_score(const output_type& pose) const {
	switch (m_polar_geometry) {
		case polar_mode::distance: return pose.polar_distance_score;
		case polar_mode::directional: return pose.polar_directional_score;
		case polar_mode::penalty: return pose.polar_penalty_score;
		case polar_mode::off: return pose.e;
	}
	return pose.e;
}

void Vina::rank_polar_poses(output_container& poses) {
	const polar_geometry geometry(m_model);
	const fl zero = m_scoring_function->conf_independent(m_model, 0);
	bool has_reward = false;
	VINA_FOR_IN(i, poses) {
		output_type& pose = poses[i];
		m_model.set(pose.c);
		const auto terms = geometry.evaluate(m_model);
		pose.native_rank = i + 1;
		pose.polar_radial = terms.radial;
		pose.polar_angular = terms.angular;
		pose.polar_distance_score = pose.e - (m_scoring_function->conf_independent(m_model, terms.radial) - zero);
		pose.polar_directional_score = pose.e - (m_scoring_function->conf_independent(m_model, terms.angular) - zero);
		if (m_polar_geometry == polar_mode::penalty)
			pose.polar_penalty_score = pose.e + (m_scoring_function->conf_independent(m_model, terms.radial - terms.angular) - zero);
		has_reward |= polar_rank_score(pose) != pose.e;
	}
	if (has_reward) {
		poses.sort([this](const output_type& a, const output_type& b) {
			const fl first = polar_rank_score(a);
			const fl second = polar_rank_score(b);
			return first == second ? a.native_rank < b.native_rank : first < second;
		});
	}
}

const char* Vina::pose_template_stage_name() const {
	switch (m_template_stage) {
		case template_stage::rank: return "rank";
		case template_stage::contact_rank: return "contact_rank";
		case template_stage::select: return "select";
		case template_stage::refine: return "refine";
	}
	return "";
}

void Vina::rank_template_poses(output_container& poses, pose_template& reference, sz& rank_evals) {
	if (m_template_stage == template_stage::refine) {
		// The retained poses were ordered by the combined refinement score.
		// Recover their native order before giving each pose its final rank.
		poses.sort();
	}
	const fl transformed_zero = m_scoring_function->conf_independent(m_model, 0);
	VINA_FOR_IN(i, poses) {
		output_type& pose = poses[i];
		pose.native_rank = i + 1;
		if (pose.template_overlap < 0) {
			m_model.set(pose.c);
			pose.template_overlap = m_model.template_overlap(reference);
			++rank_evals;
		}
		const fl reward = m_scoring_function->conf_independent(m_model, pose_template::strength * pose.template_overlap)
		                - transformed_zero;
		pose.template_rank_score = pose.e - reward;
	}
	poses.sort([](const output_type& a, const output_type& b) {
		return a.template_rank_score == b.template_rank_score ? a.native_rank < b.native_rank
		                                                       : a.template_rank_score < b.template_rank_score;
	});
}

void Vina::global_search(const int exhaustiveness, const int n_poses, const double min_rmsd, const int max_evals,
                         boost::optional<int> refine_candidates) {
	// Vina search (Monte-carlo and local optimization)
	// Check if ff, box and ligand were initialized
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the global search. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the global search. Affinity maps were not initialized.");
	} else if (exhaustiveness < 1) {
		throw vina_runtime_error("Exhaustiveness must be 1 or greater");
	}
	if (m_armijo_refinement && (!m_receptor_initialized || m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("Armijo refinement requires one ligand and a rigid receptor");
	if ((m_mixed_reorientation || m_centered_rotation) &&
	    (!m_use_explicit_receptor || m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("rigid_reorientation requires one ligand and a rigid explicit receptor");
	if (m_separate_hunt_chain && (!m_use_explicit_receptor || m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("separate_hunt_chain requires one ligand and a rigid explicit receptor");
	if (m_pocket_start_radius > 0) {
		if (!m_receptor_initialized || m_model.num_ligands() != 1 || m_model.num_flex() != 0)
			throw vina_runtime_error("pocket_start_radius requires one ligand and a rigid receptor");
		const vec width = m_grid.corner2() - m_grid.corner1();
		if (m_pocket_start_radius > 0.5 * std::min({width[0], width[1], width[2]}))
			throw vina_runtime_error("pocket_start_radius exceeds the smallest search-box half-width");
	}
	const bool late_selection = m_late_pose_selection;
	fl output_radius = 0;
	if (m_use_min_rmsd_for_merge &&
	    (late_selection || m_sf_choice != SF_VINA || !m_use_explicit_receptor ||
	     m_model.num_ligands() != 1 || m_model.num_flex() != 0 ||
	     !std::isfinite(min_rmsd) || min_rmsd < 0))
		throw vina_runtime_error("use_min_rmsd_for_merge requires rigid single-ligand native Vina docking and finite non-negative min_rmsd");
	if (late_selection && (!m_receptor_initialized || m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("late pose selection requires one ligand and a rigid receptor");
	if (late_selection && (n_poses <= 0 || !std::isfinite(min_rmsd) || min_rmsd < 0))
		throw vina_runtime_error("late pose selection requires positive num_modes and finite non-negative min_rmsd");
	const bool rerank = m_rerank_gauss2_scale != 1.0;
	if (rerank && (m_model.num_ligands() != 1 || m_model.num_flex() != 0))
		throw vina_runtime_error("rerank_gauss2_scale below 1 requires one ligand and a rigid receptor");
	const bool polar_rank = m_polar_geometry != polar_mode::off;
	const bool burial_rank = m_polar_burial != burial_mode::off;
	const bool net_polar_rank = m_net_polar != net_polar_mode::off;
	if (net_polar_rank && (!late_selection || burial_rank || !supports_polar_burial()))
		throw vina_runtime_error("net_polar requires rigid single-ligand Late Vina docking with native refinement and no other experimental scoring");
	if (burial_rank && (!late_selection || net_polar_rank || !supports_polar_burial()))
		throw vina_runtime_error("polar_burial requires rigid single-ligand Late Vina docking with native refinement and no other experimental scoring");
	if (polar_rank && (!late_selection || rerank || !m_pose_template_path.empty()))
		throw vina_runtime_error("polar_geometry requires rigid single-ligand Late Vina docking without template or Gaussian rescaling");
	if (refine_candidates) {
		if (*refine_candidates <= 0 || *refine_candidates < n_poses)
			throw vina_runtime_error("refine_candidates must be positive and at least num_modes");
		if (!late_selection || rerank || m_armijo_refinement || m_sf_choice != SF_VINA || m_no_refine ||
		    polar_rank || !m_pose_template_path.empty())
			throw vina_runtime_error("refine_candidates requires Late Vina docking with native refinement and no template, polar geometry, or Gaussian rescaling");
	}
	const bool archive_enabled = m_trajectory_archive != archive_mode::off;
	const bool archive_support = m_trajectory_archive == archive_mode::support;
	if (m_mixed_cutoff_search && (!late_selection || exhaustiveness < 2 || !m_use_explicit_receptor ||
	    m_sf_choice != SF_VINA || m_scoring_function->smooth_cutoff() || m_scoring_function->get_cutoff() <= 1 ||
	    rerank || m_armijo_refinement || polar_rank || burial_rank || net_polar_rank ||
	    !m_pose_template_path.empty() || archive_enabled))
		throw vina_runtime_error("mixed_cutoff_search requires native Late Vina refinement and ranking, a rigid receptor, exhaustiveness >= 2, and cutoff > 1 Angstrom");
	if (m_cubic_search_grid && (!late_selection || !m_use_explicit_receptor || m_sf_choice != SF_VINA ||
	                            m_soft_receptor_search || m_mixed_cutoff_search))
		throw vina_runtime_error("cubic_search_grid requires rigid Late Vina with a single receptor search grid");
	if (m_soft_receptor_search && (!late_selection || exhaustiveness < 2 || !m_use_explicit_receptor ||
	    m_sf_choice != SF_VINA || m_no_refine || m_mixed_cutoff_search || m_scoring_function->smooth_cutoff() ||
	    !m_search_anchor_path.empty() || !m_search_amide_path.empty()))
		throw vina_runtime_error("soft_receptor_search requires native Late Vina with no other search grid, a rigid receptor and at least two trajectories");
	if (m_prepared_torsion_starts && exhaustiveness < 2)
		throw vina_runtime_error("prepared_torsion_starts requires at least two trajectories");
	const bool split_anchor = m_anchor_trajectories != anchor_trajectories::all;
	if (split_anchor && exhaustiveness < 2)
		throw vina_runtime_error("search_anchor_trajectories even/odd requires exhaustiveness >= 2");
	if (!m_search_anchor_path.empty() || !m_search_amide_path.empty() || !m_rank_amide_path.empty()) {
		if (!late_selection || !m_use_explicit_receptor || m_sf_choice != SF_VINA ||
		    m_no_refine || m_mixed_cutoff_search || m_scoring_function->smooth_cutoff() ||
		    rerank || m_armijo_refinement || polar_rank || burial_rank || net_polar_rank ||
		    !m_pose_template_path.empty() || archive_enabled)
			throw vina_runtime_error("search_anchor, search_amide and rank_amide require rigid single-ligand Late Vina with native refinement and ranking");
	}
	std::vector<amide_search_site> rank_amides;
	if (!m_rank_amide_path.empty()) {
		ifile file(m_rank_amide_path);
		rank_amides = read_search_amides(file, m_model);
		// Keep a zero-barrier control byte-for-byte identical to native Late output.
		rank_amides.erase(std::remove_if(rank_amides.begin(), rank_amides.end(),
			[](const amide_search_site& site) { return site.barrier == 0; }), rank_amides.end());
	}
	const cache::cubic_view cubic_grid(m_grid);
	const igrid& base_grid = m_cubic_search_grid ? static_cast<const igrid&>(cubic_grid) : static_cast<const igrid&>(m_grid);
	boost::optional<amide_search_grid> amide_grid;
	if (!m_search_amide_path.empty()) {
		ifile file(m_search_amide_path);
		auto sites = read_search_amides(file, m_model);
		if (m_verbosity > 0 || m_search_counters)
			std::cout << "Search amide: sites=" << sites.size() << " search=all_trajectories final=native_late\n";
		amide_grid.emplace(base_grid, std::move(sites));
	}
	// Anchor trajectory selection affects only the anchor. Amides guide both subsets.
	const igrid& unanchored_grid = amide_grid ? static_cast<const igrid&>(*amide_grid) : base_grid;
	boost::optional<anchored_search_grid> anchor_grid;
	if (!m_search_anchor_path.empty()) {
		ifile file(m_search_anchor_path);
		auto sites = read_search_anchors(file, m_model, m_grid.corner1(), m_grid.corner2());
		if (m_verbosity > 0 || m_search_counters)
			std::cout << "Search anchor: sites=" << sites.front().size() << " assignments=" << sites.size()
			          << (m_anchor_trajectories == anchor_trajectories::even ? " search=even_1based" :
			              m_anchor_trajectories == anchor_trajectories::odd ? " search=odd_1based" : " search=all_trajectories")
			          << " final=native_late\n";
		anchor_grid.emplace(unanchored_grid, std::move(sites));
	}
	if (archive_enabled) {
		if (!late_selection || rerank || m_armijo_refinement || m_sf_choice != SF_VINA || m_no_refine ||
		    polar_rank || !m_pose_template_path.empty())
			throw vina_runtime_error("trajectory_archive requires rigid single-ligand Late Vina docking with native refinement and no template, polar geometry, or Gaussian rescaling");
		VINA_FOR(i, m_model.num_atoms()) {
			const atom a = m_model.get_atom(i);
			if (a.el == EL_TYPE_Dummy || (a.ad >= AD_TYPE_CG0 && a.ad <= AD_TYPE_CG3))
				throw vina_runtime_error("trajectory_archive does not support ligand dummy atoms or macrocycle closure types");
		}
	}
	std::unique_ptr<pose_template> template_reference;
	if (!m_pose_template_path.empty()) {
		const bool contact = m_template_stage == template_stage::contact_rank;
		if (contact) {
			if (m_sf_choice != SF_VINA || !m_use_explicit_receptor || m_model.num_ligands() != 1 ||
			    m_model.num_flex() != 0 || rerank || n_poses <= 0 || !std::isfinite(min_rmsd) || min_rmsd < 0)
				throw vina_runtime_error("contact_rank requires one rigid Vina receptor, one ligand, explicit refinement, positive num_modes, finite nonnegative min_rmsd, and no Gaussian rescaling");
		} else if (!late_selection || rerank || m_armijo_refinement || m_sf_choice != SF_VINA || m_no_refine)
			throw vina_runtime_error("pose_template requires Late Vina docking with native refinement and no Gaussian rescaling");
		ifile input(m_pose_template_path);
		const atomv* receptor = contact ? &m_model.grid_atoms() : nullptr;
		template_reference = std::make_unique<pose_template>(input, m_model.chemical_heavy_atom_elements(), receptor);
	}
	const bool template_selection = template_reference && (m_template_stage == template_stage::select || m_template_stage == template_stage::refine);
	const bool template_refinement = template_reference && m_template_stage == template_stage::refine;
	sz template_selection_evals = 0;
	sz template_refinement_evals = 0;
	sz template_rank_evals = 0;

	if (exhaustiveness < m_cpu) {
		std::cerr << "WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.\n";
	}

	double e = 0;
	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);
	model best_model;
	output_container poses;
	std::vector<trajectory_archive> archives;
	std::vector<search_trajectory_counters> trajectory_counters;
	auto* counters = m_search_counters ? &trajectory_counters : nullptr;
	sz explicit_inputs = 0;
	sz explicit_calls = 0;
	int explicit_evals = 0;
	std::stringstream sstm;
	rng generator(static_cast<rng::result_type>(m_seed));

	// Setup Monte-Carlo search
	parallel_mc parallelmc;
	sz heuristic = m_model.num_movable_atoms() + 10 * m_model.get_size().num_degrees_of_freedom();
	parallelmc.mc.global_steps = unsigned(70 * 3 * (50 + heuristic) / 2); // 2 * 70 -> 8 * 20 // FIXME
	parallelmc.mc.local_steps = unsigned((25 + m_model.num_movable_atoms()) / 3);
	parallelmc.mc.max_evals = max_evals;
	parallelmc.mc.min_rmsd = min_rmsd;
	parallelmc.mc.num_saved_mins = n_poses;
	parallelmc.mc.hunt_cap = vec(10, 10, 10);
	parallelmc.mc.mixed_reorientation = m_mixed_reorientation;
	parallelmc.mc.separate_hunt_chain = m_separate_hunt_chain;
	parallelmc.mc.centered_rotation = m_centered_rotation;
	parallelmc.mc.balanced_torsion = m_balanced_torsion;
	parallelmc.prepared_torsion_starts = m_prepared_torsion_starts;
	parallelmc.num_tasks = exhaustiveness;
	parallelmc.num_threads = m_cpu;
	parallelmc.display_progress = (m_verbosity > 0);
	parallelmc.keep_trajectory_minima = late_selection;
	if (m_sf_choice == SF_AD42 || !m_use_explicit_receptor)
		parallelmc.final_map_domain = m_sf_choice == SF_AD42 ? m_ad4grid.gd() : m_grid.gd();
	parallelmc.use_min_rmsd_for_merge = m_use_min_rmsd_for_merge;
	parallelmc.pocket_start_radius = m_pocket_start_radius;

	// Docking search
	ensure_vina_maps();
	sstm << "Performing docking (random seed: " << m_seed << ")";
	doing(sstm.str(), m_verbosity, 0);
	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		boost::optional<precalculate_byatom> smooth_pairs;
		boost::optional<cache> soft_receptor_grid;
		boost::optional<search_landscape> alternate;
		const igrid& guided_grid = anchor_grid ? static_cast<const igrid&>(*anchor_grid) : unanchored_grid;
		const bool guide_even_only = m_anchor_trajectories == anchor_trajectories::even;
		if (m_mixed_cutoff_search) {
			ScoringFunction smooth(SF_VINA, m_weights, m_scoring_function->get_cutoff(), true);
			if (!m_mixed_search_grid) m_mixed_search_grid.emplace(m_grid.gd());
			const auto missing = m_mixed_search_grid->missing_map_types(m_model.get_movable_atom_types(atom_type::XS));
			if (missing.any()) m_non_cache.populate_maps_from(*m_mixed_search_grid, m_receptor, missing, smooth);
			smooth_pairs.emplace(smooth, m_model);
			alternate.emplace(search_landscape{*smooth_pairs, *m_mixed_search_grid});
			if (m_verbosity > 0 || m_search_counters)
				std::cout << "\nMixed cutoff search: native=" << (exhaustiveness + 1) / 2
				          << " smooth=" << exhaustiveness / 2 << " smooth_on_even_1based=true final=native_late\n";
		} else if (m_soft_receptor_search) {
			// Relax receptor contacts only. Intraligand sterics and the final
			// explicit receptor energy must retain the original potential.
			flv weights = m_weights;
			weights[2] *= 0.5;
			ScoringFunction soft(SF_VINA, weights, m_scoring_function->get_cutoff());
			soft_receptor_grid.emplace(m_grid.gd());
			const auto types = soft_receptor_grid->missing_map_types(m_model.get_movable_atom_types(atom_type::XS));
			m_non_cache.populate_maps_from(*soft_receptor_grid, m_receptor, types, soft);
			alternate.emplace(search_landscape{m_precalculated_byatom, *soft_receptor_grid});
			if (m_verbosity > 0 || m_search_counters)
				std::cout << "\nSoft receptor search: even_1based=true repulsion_scale=0.5 stage=search_grid_only\n";
		} else if (split_anchor) {
			alternate.emplace(search_landscape{m_precalculated_byatom, guide_even_only ? guided_grid : unanchored_grid});
		}
		const igrid& search_grid = guide_even_only ? unanchored_grid : guided_grid;
		parallelmc(m_model, poses, m_precalculated_byatom, search_grid, m_grid.corner1(), m_grid.corner2(), generator,
		           m_progress_callback, archive_enabled ? &archives : nullptr, archive_support, counters,
		           alternate.get_ptr());
	} else {
		parallelmc(m_model, poses, m_precalculated_byatom, m_ad4grid, m_ad4grid.corner1(), m_ad4grid.corner2(), generator, m_progress_callback, nullptr, false, counters);
	}
	done(m_verbosity, 1);

	// Docking post-processing and rescoring
	const sz saved_minima = poses.size();
	const sz accepted_minima = archive_enabled ? append_accepted_minima(poses, archives) : 0;
	const sz prescore_inputs = poses.size();
	if (late_selection) {
		// Use the saved conformation with chemical atoms only for selection geometry.
		// Score with native weights and the same box penalty used for final refinement ranking.
		m_non_cache.slope = 1e6;
		for (auto& pose : poses) {
			m_model.set(pose.c);
			pose.coords = m_model.get_chemical_heavy_atom_movable_coords();
			pose.e = m_non_cache.eval(m_model, authentic_v[1])
			       + m_model.evali(m_precalculated_byatom, authentic_v);
			if (template_selection) {
				pose.e -= pose_template::strength * m_model.template_overlap(*template_reference);
				++template_selection_evals;
			}
		}
		const sz refinement_limit = refine_candidates.value_or(n_poses);
		poses = select_distinct_poses(poses, min_rmsd, refinement_limit);
	} else {
		poses = remove_redundant(poses, min_rmsd);
	}

	if (!poses.empty()) {
		// For the Vina scoring function, we take the intramolecular energy from the best pose
		// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			// Refine poses if no_refine is false and got receptor
			if (m_use_explicit_receptor) {
				const sz refinement_inputs = poses.size();
				const sz accepted_refinement_inputs = archive_enabled
					? std::count_if(poses.begin(), poses.end(), [](const output_type& pose) {
						return pose.archive_from_accepted;
					})
					: 0;
				sz refinement_calls = 0;
				quasi_newton quasi_newton_par(m_model.get_size(),
				                             unsigned((25 + m_model.num_movable_atoms()) / 3));
				quasi_newton_par.strict_armijo = m_armijo_refinement;
				quasi_newton_par.geometry_metric = m_geometry_refinement;
				int evalcount = 0;
				const fl slope = 1e6;
				m_non_cache.slope = slope;

				VINA_FOR_IN(i, poses){
					VINA_FOR(p, 5){
						m_non_cache.slope = 100 * std::pow(10.0, 2.0*p);
						quasi_newton_par(m_model, m_precalculated_byatom, m_non_cache, poses[i], authentic_v, evalcount,
						                 template_refinement ? template_reference.get() : nullptr);
						++refinement_calls;
						if(m_non_cache.within(m_model))
							break;
					}
					poses[i].coords = late_selection ? m_model.get_chemical_heavy_atom_movable_coords()
					                                : m_model.get_heavy_atom_movable_coords();
					if (late_selection && m_model.num_ligands() == 1 && m_model.num_flex() == 0) {
						const vec shift = printable_box_shift(m_model.get_grid_atom_movable_coords(),
															  m_grid.corner1(), m_grid.corner2());
						if (sqr(shift) > 0) {
							poses[i].c.ligands[0].rigid.position += shift;
							m_model.set(poses[i].c);
							poses[i].coords = m_model.get_chemical_heavy_atom_movable_coords();
						}
					}
					m_non_cache.slope = slope;
					// rescoring in case a ligand or flex sidechain atom is outside box
					// ensuring poses will be sorted with the same slope (a.k.a. out of
					// box penalty) that will be used to calculate final energies.
					double all_grids = m_non_cache.eval(m_model, authentic_v[1]);
					double inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // ligand -- flex
					double intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // flex_i -- flex_i and flex_i -- flex_j
					double lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // ligand_i -- ligand_i
					poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra;
					if (template_refinement) {
						output_type& pose = poses[i];
						pose.native_energy = pose.e;
						pose.template_overlap = m_model.template_overlap(*template_reference);
						pose.e -= pose_template::strength * pose.template_overlap;
						++template_rank_evals;
					}
				}
				explicit_inputs = refinement_inputs;
				explicit_calls = refinement_calls;
				explicit_evals = evalcount;
				if (template_refinement)
					template_refinement_evals = evalcount;
				if (late_selection) {
					poses = select_distinct_poses(poses, min_rmsd, n_poses);
					if (m_verbosity > 0)
						std::cout << "Late pose selection: saved=" << saved_minima
						          << " refine=" << refinement_inputs << " refinement_calls=" << refinement_calls
						          << " retained=" << poses.size() << '\n';
					if (archive_enabled && m_verbosity > 0) {
						sz search_evals = 0;
						sz search_refinements = 0;
						for (const auto& archive : archives) {
							search_evals += archive.native_evaluations;
							search_refinements += archive.native_refinements;
						}
						std::cout << "Trajectory archive: mode=" << trajectory_archive_name()
						          << " trajectories=" << archives.size() << " native=" << saved_minima
						          << " extra=" << accepted_minima << " prescore=" << prescore_inputs
						          << " extra_prescore=" << prescore_inputs - saved_minima
						          << " refine=" << refinement_inputs << " extra_refine=" << accepted_refinement_inputs
						          << " refinement_calls=" << refinement_calls << " retained=" << poses.size()
						          << " native_search_evals=" << search_evals
						          << " native_search_refinements=" << search_refinements << '\n';
					}
				}
			}

			if (!late_selection)
				poses.sort(); // late selection already stably sorted on refined energy
			const output_type& native_best = template_refinement
				? *std::min_element(poses.begin(), poses.end(), [](const output_type& a, const output_type& b) {
					return a.native_energy < b.native_energy;
				})
				: poses[0];
			m_model.set(native_best.c);
			intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, vina_receptor_interactions(), authentic_v);
		}

		VINA_FOR_IN(i, poses) {
			if (m_verbosity > 1)
				std::cout << "ENERGY FROM SEARCH: " << poses[i].e << "\n";

			m_model.set(poses[i].c);

			// For AD42 intramolecular_energy is equal to 0
			std::vector<double> energies = score(intramolecular_energy);
			// Store energy components in current pose
			poses[i].e = energies[0]; // specific to each scoring function
			poses[i].inter = energies[1] + energies[2];
			poses[i].intra = energies[3] + energies[4] + energies[5];
			poses[i].total = poses[i].inter + poses[i].intra; // cost function for optimization
			poses[i].conf_independent = energies[6]; // "torsion"
			poses[i].unbound = energies[7]; // specific to each scoring function

			if (m_verbosity > 1) {
				std::cout << (rerank ? "NATIVE FINAL ENERGY: \n" : "FINAL ENERGY: \n");
				show_score(energies);
			}
		}

		// In AD4, the unbound energy is intra for each pose, so order may have changed
		// The order does not change in Vina because unbound is intra of the 1st pose
		if (m_sf_choice == SF_AD42)
			poses.sort();
		if (rerank)
			rerank_poses(poses, authentic_v);
		if (!rank_amides.empty()) {
			VINA_FOR_IN(i, poses) {
				auto& pose = poses[i];
				m_model.set(pose.c);
				pose.native_rank = i + 1;
				pose.amide_penalty = amide_search_grid::penalty(m_model, rank_amides);
			}
			poses.sort([](const output_type& a, const output_type& b) {
				return a.amide_rank_score() == b.amide_rank_score() ? a.native_rank < b.native_rank
				                                                  : a.amide_rank_score() < b.amide_rank_score();
			});
		}
		if (polar_rank) {
			rank_polar_poses(poses);
			if (m_verbosity > 0)
				std::cout << "Polar geometry scoring: mode=" << polar_geometry_name() << " poses=" << poses.size() << '\n';
		}
		if (burial_rank) {
			rank_polar_burial_poses(poses);
			if (m_verbosity > 0)
				std::cout << "Polar burial scoring: mode=" << polar_burial_name() << " poses=" << poses.size() << '\n';
		}
		if (net_polar_rank) {
			rank_net_polar_poses(poses);
			if (m_verbosity > 0)
				std::cout << "Net polar scoring: mode=" << net_polar_name() << " poses=" << poses.size() << '\n';
		}
		if (template_reference) {
			// Freeze the same native-order output cap before scoring. Native
			// merge and Late selection already impose this cap on both paths.
			if (m_template_stage == template_stage::contact_rank && poses.size() > sz(n_poses))
				poses.erase(poses.begin() + n_poses, poses.end());
			rank_template_poses(poses, *template_reference, template_rank_evals);
			if (m_verbosity > 0)
				std::cout << "Template pose scoring: stage=" << pose_template_stage_name()
				          << " selection_evals=" << template_selection_evals
				          << " refinement_evals=" << template_refinement_evals
				          << " rank_evals=" << template_rank_evals << '\n';
		}
		if (archive_enabled) {
			// Late refinement already has stable tie order. Rank the final
			// reported native energies explicitly at their full precision.
			VINA_FOR_IN(i, poses)
				poses[i].native_rank = i + 1;
			poses.sort([](const output_type& a, const output_type& b) {
				return a.e == b.e ? a.native_rank < b.native_rank : a.e < b.e;
			});
			const fl eligible_energy = poses.front().e + 1.0;
			VINA_FOR_IN(i, poses) {
				output_type& pose = poses[i];
				pose.native_rank = i + 1;
				pose.archive_eligible = pose.e <= eligible_energy;
				pose.archive_votes = archive_support ? archived_trajectory_support(pose.coords, archives) : 0;
			}
			if (archive_support) {
				poses.sort([](const output_type& a, const output_type& b) {
					if (a.archive_eligible != b.archive_eligible)
						return a.archive_eligible;
					if (a.archive_eligible && a.archive_votes != b.archive_votes)
						return a.archive_votes > b.archive_votes;
					return a.native_rank < b.native_rank;
				});
			}
		}

		// Native explicit optimization moves the representatives after the grid-stage merge.
		// Native energy order makes every energy window a prefix, so prune once in place.
		// A separately ranked native window can skip earlier poses: defer its distance test
		// until the caller's window is known, so excluded poses never suppress eligible ones.
		if (!late_selection && explicit_inputs && min_rmsd > 0) {
			if (uses_native_energy_window(poses.front())) {
				output_radius = min_rmsd;
			} else {
				const auto keep = ranked_output_indices(poses, poses.size(), min_rmsd, boost::none);
				if (keep.size() != poses.size()) {
					sz remaining = keep.size();
					for (sz i = poses.size(); i > 0; --i) {
						if (remaining && keep[remaining - 1] == i - 1) --remaining;
						else poses.erase(poses.begin() + (i - 1));
					}
				}
			}
		}

		// Now compute RMSD from the best model
		// Necessary to do it in two pass for AD4 scoring function
		m_model.set(poses[0].c);
		best_model = m_model;

		if (m_verbosity > 0) {
			std::cout << '\n';
			if (rerank)
				std::cout << "Pose ranking: FastDock Vina Pro\n";
			if (!rank_amides.empty() || template_reference || polar_rank || burial_rank || net_polar_rank) {
				std::cout << "mode | " << (!rank_amides.empty() ? "amide score" : net_polar_rank ? "net polar score" : burial_rank ? "polar burial score" : polar_rank ? "polar score" : "template score") << " | dist from best mode\n";
				std::cout << "     | (rank score)   | rmsd l.b.| rmsd u.b.\n";
			} else {
				std::cout << "mode |   affinity | dist from best mode\n";
				std::cout << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
			}
			std::cout << "-----+------------+----------+----------\n";
		}

		const auto visible = m_verbosity > 0
		    ? ranked_output_indices(poses, poses.size(), output_radius, boost::none) : pose_indices();
		sz reported = 0;
		VINA_FOR_IN(i, poses) {
			m_model.set(poses[i].c);

			// Get RMSD between current pose and best_model
			poses[i].lb = m_model.rmsd_lower_bound(best_model);
			poses[i].ub = m_model.rmsd_upper_bound(best_model);

			if (reported < visible.size() && i == visible[reported]) {
				++reported;
				const fl reported_score = !rank_amides.empty() ? poses[i].amide_rank_score()
				                        : net_polar_rank ? poses[i].net_polar_score->rank_score
				                        : burial_rank ? poses[i].polar_burial_rank_score
				                        : polar_rank ? polar_rank_score(poses[i])
				                        : template_reference ? poses[i].template_rank_score : poses[i].e;
				std::cout << std::setw(4) << reported << "    " << std::setw(9) << std::setprecision(4) << reported_score;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].lb;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].ub << "\n";
			}
		}

		// Clean up by putting back the best pose in model
		m_model.set(poses[0].c);
	} else {
		if (parallelmc.final_map_domain)
			std::cerr << "WARNING: No saved poses lie inside the affinity map domain.\n";
		else
			std::cerr << "WARNING: Zero poses in output container after global search.\n";
	}

	if (m_search_counters) {
		if (!late_selection)
			std::cout << "Native merge: rmsd=" << (m_use_min_rmsd_for_merge ? min_rmsd : 2)
			          << " candidates=" << saved_minima << '\n';
		VINA_FOR_IN(i, trajectory_counters) {
			const auto& counts = trajectory_counters[i];
			std::cout << "Search trajectory: index=" << i + 1 << " steps=" << counts.hunt_calls
			          << " hunt_calls=" << counts.hunt_calls << " hunt_evals=" << counts.hunt_evals
			          << " second_calls=" << counts.second_calls << " second_evals=" << counts.second_evals
			          << " rigid_reorientations=" << counts.rigid_reorientations << '\n';
			std::cout << "Hunt chain: index=" << i + 1 << " separate=" << int(m_separate_hunt_chain)
			          << " accepted=" << counts.accepted_steps
			          << " mixed_accept_checks=" << counts.mixed_accept_checks << '\n';
			if (counts.initial_centroid) {
				const auto& center = *counts.initial_centroid;
				std::ostringstream line;
				line << "Pocket start: index=" << i + 1 << " anchored=" << int(counts.pocket_anchored)
				     << std::fixed << std::setprecision(12)
				     << " cx=" << center[0] << " cy=" << center[1] << " cz=" << center[2] << '\n';
				std::cout << line.str();
			}
		}
		const char* reorientation = m_centered_rotation ? "centered" : (m_mixed_reorientation ? "mixed" : "native");
		std::cout << "Search counters: rigid_reorientation=" << reorientation
		          << " explicit_inputs=" << explicit_inputs << " explicit_calls=" << explicit_calls
		          << " explicit_evals=" << explicit_evals << " late_prescore=" << (late_selection ? prescore_inputs : 0)
		          << '\n';
	}

	// Store results in Vina object
	m_poses.swap(poses);
	m_search_completed = true;
	m_output_min_rmsd = output_radius;
}

Vina::~Vina() = default;
