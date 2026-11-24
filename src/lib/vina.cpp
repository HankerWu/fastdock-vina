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
   to avoid deep copies without changing docking behavior.
   See docs/PROVENANCE.md for upstream and experimental source identities. */

#include <utility>

#include "vina.h"
#include "scoring_function.h"
#include "precalculate.h"


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
	m_receptor_initialized = true;
	// If we are reading another receptor we should not consider the ligand and the map as initialized anymore
	m_ligand_initialized = false;
	m_map_initialized = false;
	if (m_batch_map_scoring && m_no_refine)
		m_non_cache = non_cache();
	m_batch_map_scoring.reset();
}

void Vina::reset_ligand_model(atom_type::t atom_typing) {
	if (m_receptor_initialized)
		m_model = m_receptor;
	else
		m_model = model(atom_typing);
}

void Vina::finish_ligand_setup(atom_type::t atom_typing) {
	precalculate_byatom precalculated_byatom(*m_scoring_function, m_model);
	if (m_map_initialized) {
		const szv atom_types = m_model.get_movable_atom_types(atom_typing);
		if (m_batch_map_scoring) {
			const auto missing = m_grid.missing_map_types(atom_types);
			if (missing.any()) {
				const auto& index = m_non_cache.receptor_index();
				if (m_no_refine) {
					const precalculate tables(*m_batch_map_scoring, index.atom_types(), missing);
					m_grid.populate(m_receptor, tables, index, missing);
				} else {
					m_grid.populate(m_receptor, m_non_cache.receptor_tables(), index, missing);
				}
			}
		} else if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			if (!m_grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		} else {
			if (!m_ad4grid.are_atom_types_grid_initialized(atom_types))
				exit(EXIT_FAILURE);
		}
	}

	m_poses.clear();
	m_precalculated_byatom = std::move(precalculated_byatom);
	m_ligand_initialized = true;
}

void Vina::set_ligand_from_string(const std::string& ligand_string) {
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand file. Ligand string is empty.\n";
		exit(EXIT_FAILURE);
	}
	const atom_type::t atom_typing = m_scoring_function->get_atom_typing();
	reset_ligand_model(atom_typing);
	m_model.append(parse_ligand_pdbqt_from_string(ligand_string, atom_typing));
	finish_ligand_setup(atom_typing);
}

void Vina::set_ligand_from_string(const std::vector<std::string>& ligand_string) {
	if (ligand_string.empty()) {
		std::cerr << "ERROR: Cannot read ligand list. Ligands list is empty.\n";
		exit(EXIT_FAILURE);
	}
	const atom_type::t atom_typing = m_scoring_function->get_atom_typing();
	reset_ligand_model(atom_typing);
	for (const auto& ligand : ligand_string)
		m_model.append(parse_ligand_pdbqt_from_string(ligand, atom_typing));
	finish_ligand_setup(atom_typing);
}

void Vina::set_ligand_from_file(const std::string& ligand_name) {
	set_ligand_from_string(get_file_contents(ligand_name));
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
    m_scoring_function = std::make_shared<ScoringFunction>(m_sf_choice, m_weights, cutoff);
}

// FastDock: bind the cutoff before any ligand or map lookup tables are created.
void Vina::set_interaction_cutoff(double cutoff) {
	if (m_sf_choice == SF_AD42)
		throw vina_runtime_error("interaction_cutoff is supported by Vina and Vinardo");
	if (m_ligand_initialized || m_map_initialized)
		throw vina_runtime_error("set interaction_cutoff before loading ligands or maps");
	m_scoring_function = std::make_shared<ScoringFunction>(m_sf_choice, m_weights, cutoff);
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

void Vina::compute_vina_maps(double center_x, double center_y, double center_z, double size_x, double size_y, double size_z, double granularity, bool force_even_voxels) {
	build_vina_maps(vec(center_x, center_y, center_z), vec(size_x, size_y, size_z),
	                granularity, force_even_voxels, false);
}

void Vina::initialize_vina_batch_maps(double center_x, double center_y, double center_z,
                                     double size_x, double size_y, double size_z, double granularity) {
	build_vina_maps(vec(center_x, center_y, center_z), vec(size_x, size_y, size_z),
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

	if (m_sf_choice == SF_VINA)
		doing("Computing Vina grid", m_verbosity, 0);
	else
		doing("Computing Vinardo grid", m_verbosity, 0);

	cache grid(gd, slope);
	szv_grid receptor_index(m_receptor, szv_grid_dims(gd), sqr(m_scoring_function->get_cutoff()));
	std::bitset<XS_TYPE_SIZE> map_types;
	if (!defer_population || m_ligand_initialized) {
		map_types = cache::map_types(m_ligand_initialized
			? m_model.get_movable_atom_types(m_scoring_function->get_atom_typing())
			: m_scoring_function->get_atom_types());
	}
	precalculate receptor_tables;
	if (!m_no_refine || map_types.any()) {
		auto interaction_types = map_types;
		// Explicit scoring also uses raw XS types and supports replacing the ligand.
		if (!m_no_refine) interaction_types.set();
		receptor_tables = precalculate(*m_scoring_function, receptor_index.atom_types(), interaction_types);
		grid.populate(m_receptor, receptor_tables, receptor_index, map_types);
	}
	done(m_verbosity, 0);

	if (!m_no_refine || defer_population) {
		if (m_no_refine) receptor_tables = precalculate();
		m_non_cache = non_cache(std::move(receptor_index), gd, std::move(receptor_tables), slope);
	} else {
		m_non_cache = non_cache();
	}
	m_grid = std::move(grid);
	m_batch_map_scoring = defer_population ? m_scoring_function : nullptr;
	m_map_initialized = true;
}

void Vina::load_maps(std::string maps) {
	if (m_sf_choice != SF_AD42 && m_scoring_function->get_cutoff() != default_xs_interaction_cutoff)
		throw vina_runtime_error("map files require the default interaction_cutoff");
	const fl slope = 1e6; // FIXME: too large? used to be 100

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		doing("Reading Vina maps", m_verbosity, 0);
		cache grid(slope);
		grid.read(maps);
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

	// Imported maps cannot be extended by the receptor used for a previous batch.
	if (m_batch_map_scoring && m_no_refine)
		m_non_cache = non_cache();
	m_batch_map_scoring.reset();
	m_map_initialized = true;
}

void Vina::write_maps(const std::string& map_prefix, const std::string& gpf_filename,
					  const std::string& fld_filename, const std::string& receptor_filename) {
	if (m_sf_choice != SF_AD42 && m_scoring_function->get_cutoff() != default_xs_interaction_cutoff)
		throw vina_runtime_error("map files require the default interaction_cutoff");
	if (!m_map_initialized) {
		throw vina_runtime_error("Cannot write affinity maps. Affinity maps were not initialized.");
	}

	szv atom_types;
	atom_type::t atom_typing = m_scoring_function->get_atom_typing();

	if (m_ligand_initialized)
		atom_types = m_model.get_movable_atom_types(atom_typing);
	else if (m_batch_map_scoring)
		throw vina_runtime_error("Load a ligand before writing deferred affinity maps.");
	else
		atom_types = m_scoring_function->get_atom_types();

	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
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

std::vector< std::vector<double> > Vina::get_poses_coordinates(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::vector< std::vector<double> > coordinates;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			m_model.set(m_poses[i].c);
			coordinates.push_back(m_model.get_ligand_coords());

			n++;
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);
	} else {
		std::cerr << "WARNING: Could not find any pose coordinaates.\n";
	}

	return coordinates;
}

std::vector< std::vector<double> > Vina::get_poses_energies(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::vector< std::vector<double> > energies;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses asked must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			energies.push_back({m_poses[i].e,
								m_poses[i].inter, m_poses[i].intra,
								m_poses[i].conf_independent, m_poses[i].unbound});

			n++;
		}
	} else {
		std::cerr << "WARNING: Could not find any pose energies.\n";
	}

	return energies;
}

std::string Vina::vina_remarks(output_type &pose, fl lb, fl ub) {
	std::ostringstream remark;

	remark.setf(std::ios::fixed, std::ios::floatfield);
	remark.setf(std::ios::showpoint);

	remark << "REMARK VINA RESULT: "
		   << std::setw(9) << std::setprecision(3) << pose.e
		   << "  " << std::setw(9) << std::setprecision(3) << lb
		   << "  " << std::setw(9) << std::setprecision(3) << ub
		   << '\n';

	remark << "REMARK INTER + INTRA:    " << std::setw(12) << std::setprecision(3) << pose.total << "\n";
	remark << "REMARK INTER:            " << std::setw(12) << std::setprecision(3) << pose.inter << "\n";
	remark << "REMARK INTRA:            " << std::setw(12) << std::setprecision(3) << pose.intra << "\n";
	if (m_sf_choice == SF_AD42)
		remark << "REMARK CONF_INDEPENDENT: " << std::setw(12) << std::setprecision(3) << pose.conf_independent << "\n";
	remark << "REMARK UNBOUND:          " << std::setw(12) << std::setprecision(3) << pose.unbound << "\n";

	return remark.str();
}

std::string Vina::get_poses(int how_many, double energy_range) {
	int n = 0;
	double best_energy = 0;
	std::ostringstream out;
	std::string remarks;

	if (how_many < 0) {
		throw vina_runtime_error("number of poses written must be greater than zero.");
	}

	if (energy_range < 0) {
		throw vina_runtime_error("energy range must be greater than zero.");
	}

	if (!m_poses.empty()) {
		// Get energy from the best conf
		best_energy = m_poses[0].e;

		VINA_FOR_IN(i, m_poses) {
			/* Stop if:
				- We wrote the number of conf asked
				- If there is no conf to write
				- The energy of the current conf is superior than best_energy + energy_range
			*/
			if (n >= how_many || !not_max(m_poses[i].e) || m_poses[i].e > best_energy + energy_range)
				break; // check energy_range sanity FIXME

			// Push the current pose to model
			m_model.set(m_poses[i].c);

			// Write conf
			remarks = vina_remarks(m_poses[i], m_poses[i].lb, m_poses[i].ub);
			// Write into the final buffer without a separate string per model.
			m_model.write_model(out, n + 1, remarks);

			n++;
		}

		// Push back the best conf in model
		m_model.set(m_poses[0].c);

	} else {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
	}

	return out.str();
}

void Vina::write_poses(const std::string& output_name, int how_many, double energy_range) {
	if (!m_poses.empty()) {
		// Open output file
		ofile f(make_path(output_name));
		f << get_poses(how_many, energy_range);
	} else {
		std::cerr << "WARNING: Could not find any poses. No poses were written.\n";
	}
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
		// Inter
		if (m_no_refine || !m_receptor_initialized)
			all_grids = m_grid.eval(m_model, authentic_v[1]); // [1] ligand & flex -- grid
		else
			all_grids = m_non_cache.eval(m_model, authentic_v[1]); // [1] ligand & flex -- grid
		inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // [1] ligand -- flex
		// Intra
		if (m_no_refine || !m_receptor_initialized)
			flex_grids = m_grid.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
		else
			flex_grids = m_non_cache.eval_intra(m_model, authentic_v[1]); // [1] flex -- grid
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

std::vector<double> Vina::score() {
	// Score the current conf in the model
	// Check if ff and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot score the pose. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot score the pose. Affinity maps were not initialized.");
	} else if (!m_grid.is_in_grid(m_model)) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);

	if(m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);
	}

	std::vector<double> energies = score(intramolecular_energy);
	return energies;
}

std::vector<double> Vina::optimize(output_type& out, int max_steps) {
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

	if (m_verbosity > 1) {
		std::cout << "Before local optimization:\n";
		energies_before_opt = score();
		show_score(energies_before_opt);
	}

	doing("Performing local search", m_verbosity, 0);
	// Try 5 five times to optimize locally the conformation
	VINA_FOR(p, 5) {
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			quasi_newton_par(m_model, m_precalculated_byatom, m_grid,    out, authentic_v, evalcount);
			// Break if we succeed to bring (back) the ligand within the grid
			if (m_grid.is_in_grid(m_model))
				break;
		} else {
			quasi_newton_par(m_model, m_precalculated_byatom, m_ad4grid, out, authentic_v, evalcount);
			if (m_ad4grid.is_in_grid(m_model))
				break;
		}
	}
	done(m_verbosity, 0);

	energies_after_opt = score();

	return energies_after_opt;
}

std::vector<double> Vina::optimize(int max_steps) {
	// Local optimization of the ligand conf
	// Check if ff, box and ligand were initialized
	// Check if the ligand is not outside the box
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the optimization. Affinity maps were not initialized.");
	} else if (!m_grid.is_in_grid(m_model)) {
		throw vina_runtime_error("The ligand is outside the grid box. Increase the size of the grid box or center it accordingly around the ligand.");
	}

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

	std::vector<double> energies = optimize(out, max_steps);

	return energies;
}

output_container Vina::remove_redundant(const output_container &in, fl min_rmsd) {
	output_container tmp;
	VINA_FOR_IN(i, in)
	add_to_output_container(tmp, in[i], min_rmsd, in.size());
	return tmp;
}

void Vina::global_search(const int exhaustiveness, const int n_poses, const double min_rmsd, const int max_evals) {
	// Vina search (Monte-carlo and local optimization)
	// Check if ff, box and ligand were initialized
	if (!m_ligand_initialized) {
		throw vina_runtime_error("Cannot do the global search. Ligand(s) was(ere) not initialized.");
	} else if (!m_map_initialized) {
		throw vina_runtime_error("Cannot do the global search. Affinity maps were not initialized.");
	} else if (exhaustiveness < 1) {
		throw vina_runtime_error("Exhaustiveness must be 1 or greater");
	}

	if (exhaustiveness < m_cpu) {
		std::cerr << "WARNING: At low exhaustiveness, it may be impossible to utilize all CPUs.\n";
	}

	double e = 0;
	double intramolecular_energy = 0;
	const vec authentic_v(1000, 1000, 1000);
	model best_model;
	boost::optional<model> ref;
	output_container poses;
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
	parallelmc.num_tasks = exhaustiveness;
	parallelmc.num_threads = m_cpu;
	parallelmc.display_progress = (m_verbosity > 0);

	// Docking search
	sstm << "Performing docking (random seed: " << m_seed << ")";
	doing(sstm.str(), m_verbosity, 0);
	if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
		parallelmc(m_model, poses, m_precalculated_byatom,    m_grid, m_grid.corner1(), m_grid.corner2(), generator, m_progress_callback);
	} else {
		parallelmc(m_model, poses, m_precalculated_byatom, m_ad4grid, m_ad4grid.corner1(), m_ad4grid.corner2(), generator, m_progress_callback);
	}
	done(m_verbosity, 1);

	// Docking post-processing and rescoring
	poses = remove_redundant(poses, min_rmsd);

	if (!poses.empty()) {
		// For the Vina scoring function, we take the intramolecular energy from the best pose
		// the order must not change because of non-decreasing g (see paper), but we'll re-sort in case g is non strictly increasing
		if (m_sf_choice == SF_VINA || m_sf_choice == SF_VINARDO) {
			// Refine poses if no_refine is false and got receptor
			if (!m_no_refine & m_receptor_initialized) {
				quasi_newton quasi_newton_par(m_model.get_size(),
				                             unsigned((25 + m_model.num_movable_atoms()) / 3));
				int evalcount = 0;
				const fl slope = 1e6;
				m_non_cache.slope = slope;

				VINA_FOR_IN(i, poses){
					VINA_FOR(p, 5){
						m_non_cache.slope = 100 * std::pow(10.0, 2.0*p);
						quasi_newton_par(m_model, m_precalculated_byatom, m_non_cache, poses[i], authentic_v, evalcount);
						if(m_non_cache.within(m_model))
							break;
					}
					poses[i].coords = m_model.get_heavy_atom_movable_coords();
					m_non_cache.slope = slope;
					// rescoring in case a ligand or flex sidechain atom is outside box
					// ensuring poses will be sorted with the same slope (a.k.a. out of
					// box penalty) that will be used to calculate final energies.
					double all_grids = m_non_cache.eval(m_model, authentic_v[1]);
					double inter_pairs = m_model.eval_inter(m_precalculated_byatom, authentic_v); // ligand -- flex
					double intra_pairs = m_model.evalo(m_precalculated_byatom, authentic_v); // flex_i -- flex_i and flex_i -- flex_j
					double lig_intra = m_model.evali(m_precalculated_byatom, authentic_v); // ligand_i -- ligand_i
					poses[i].e = all_grids + inter_pairs + intra_pairs + lig_intra;
				}
			}

			poses.sort(); // order often changes after non_cache refinement
			m_model.set(poses[0].c);
			if (m_no_refine || !m_receptor_initialized)
				intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_grid, authentic_v);
			else
				intramolecular_energy = m_model.eval_intramolecular(m_precalculated_byatom, m_non_cache, authentic_v);
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
				std::cout << "FINAL ENERGY: \n";
				show_score(energies);
			}
		}

		// In AD4, the unbound energy is intra for each pose, so order may have changed
		// The order does not change in Vina because unbound is intra of the 1st pose
		if (m_sf_choice == SF_AD42)
			poses.sort();

		// Now compute RMSD from the best model
		// Necessary to do it in two pass for AD4 scoring function
		m_model.set(poses[0].c);
		best_model = m_model;

		if (m_verbosity > 0) {
			std::cout << '\n';
			std::cout << "mode |   affinity | dist from best mode\n";
			std::cout << "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n";
			std::cout << "-----+------------+----------+----------\n";
		}

		VINA_FOR_IN(i, poses) {
			m_model.set(poses[i].c);

			// Get RMSD between current pose and best_model
			const model &r = ref ? ref.get() : best_model;
			poses[i].lb = m_model.rmsd_lower_bound(r);
			poses[i].ub = m_model.rmsd_upper_bound(r);

			if (m_verbosity > 0) {
				std::cout << std::setw(4) << i + 1 << "    " << std::setw(9) << std::setprecision(4) << poses[i].e;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].lb;
				std::cout << "  " << std::setw(9) << std::setprecision(4) << poses[i].ub << "\n";
			}
		}

		// Clean up by putting back the best pose in model
		m_model.set(poses[0].c);
	} else {
		std::cerr << "WARNING: Zero poses in output container after global search. This should not be happening and is likely a bug.\n";
		std::cerr << "WARNING: If possible, please file a bug report with your input files and random seed on GitHub.\n";
	}

	// Store results in Vina object
	m_poses.swap(poses);
}

Vina::~Vina() = default;
