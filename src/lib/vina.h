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

// FastDock: share receptor data and defer affinity maps.

#ifndef VINA_H
#define VINA_H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <exception>
#include <vector> // ligand paths
#include <cmath> // for ceila
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp> // hardware_concurrency // FIXME rm ?
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>
//#include <openbabel/mol.h>
#include "parse_pdbqt.h"
#include "parallel_mc.h"
#include "file.h"
#include "conf.h"
#include "model.h"
#include "common.h"
#include "cache.h"
#include "non_cache.h"
#include "ad4cache.h"
#include "quasi_newton.h"
#include "coords.h" // add_to_output_container
#include "utils.h"
#include "scoring_function.h"
#include "precalculate.h"
#include <memory>


class vina_runtime_error : public std::exception {
public:
    explicit vina_runtime_error(const std::string & message)
            : m_message("\n\nVina runtime error: " + message + "\n") {}

    virtual const char* what() const throw () {
        return m_message.c_str();
    }

private:
    const std::string m_message;
};

class Vina {
public:
	// Constructor
	Vina(const std::string &sf_name="vina", int cpu=0, int seed=0, int verbosity=1, bool no_refine=false, std::function<void(double)>* progress_callback = NULL) {
		m_verbosity = verbosity;
		m_receptor_initialized = false;
		m_ligand_initialized = false;
		m_map_initialized = false;
		m_seed = generate_seed(seed);
		m_no_refine = no_refine;
		m_progress_callback = progress_callback;

		// Look for the number of cpu
		if (cpu <= 0) {
			unsigned num_cpus = boost::thread::hardware_concurrency();

			if (num_cpus > 0) {
				m_cpu = num_cpus;
			} else {
				std::cerr << "WARNING: Could not determined the number of concurrent thread supported on this machine. ";
				std::cerr << "You might need to set it manually using cpu argument or fix the issue.\n";
				exit(EXIT_FAILURE);
			}
		} else {
			m_cpu = cpu;
		}

		if (sf_name.compare("vina") == 0) {
			m_sf_choice = SF_VINA;
			set_vina_weights();
		} else if (sf_name.compare("vinardo") == 0) {
			m_sf_choice = SF_VINARDO;
			set_vinardo_weights();
		} else if (sf_name.compare("ad4") == 0) {
			m_sf_choice = SF_AD42;
			set_ad4_weights();
		} else {
			std::cerr << "ERROR: Scoring function " << sf_name << " not implemented (choices: vina, vinardo or ad4)\n";
			exit (EXIT_FAILURE);
		}
	}
	// Destructor
	~Vina();

	void cite();
	int seed() { return m_seed; }
	void set_receptor(const std::string &rigid_name=std::string(), const std::string &flex_name=std::string());
	void set_ligand_from_string(const std::string &ligand_string);
	void set_ligand_from_string(const std::vector<std::string> &ligand_string);
	void set_ligand_from_file(const std::string& ligand_name);
	void set_ligand_from_file(const std::vector<std::string>& ligand_name);
	// Fixed score only. Preserve prepared topology while replacing lab coordinates.
	void set_pose_coordinates_for_scoring(const std::string& pose_name);
	void show_rigid_window_pool(const std::string& path, unsigned samples) const;
	void set_ligand_from_file_with_bonds(const std::string& ligand_name, const std::string& bonds_name);
	// Configure the Vina/Vinardo cutoff before loading ligands or maps.
	void set_interaction_cutoff(double cutoff, bool smooth_cutoff = false);
	// Configure final Vina pose ranking before building search maps.
	void set_rerank_gauss2_scale(double scale);
	// Use explicit-receptor rescoring before merging trajectories, configured before maps.
	void set_late_pose_selection(bool enabled);
	void set_use_min_rmsd_for_merge(bool enabled) { m_use_min_rmsd_for_merge = enabled; }
	// Enable experimental template guidance before building search maps.
	void set_rank_amide(const std::string& filename) { m_rank_amide_path = filename; }
	void set_search_amide(const std::string& filename) { m_search_amide_path = filename; }
	void set_search_anchor(const std::string& filename, const std::string& trajectories = "all");
	void set_pose_template(const std::string& filename, const std::string& stage);
	// Enable experimental polar-geometry final ranking before building search maps.
	void set_polar_geometry(const std::string& mode);
	void set_polar_burial(const std::string& mode);
	void show_polar_burial();
	void set_net_polar(const std::string& mode, const std::string& waters);
	void show_net_polar();
	// Enable experimental trajectory archiving before building search maps.
	void set_trajectory_archive(const std::string& mode);
	// Use strict Armijo line search only in explicit final refinement.
	void set_armijo_refinement(bool enabled, bool geometry_metric = false);
	void set_rigid_reorientation(bool mixed, bool centered = false);
	void set_separate_hunt_chain(bool enabled);
	void set_prepared_torsion_starts(bool enabled) { m_prepared_torsion_starts = enabled; }
	void set_balanced_torsion_search(bool enabled) { m_balanced_torsion = enabled; }
	void set_cubic_search_grid(bool enabled) { m_cubic_search_grid = enabled; }
	void set_soft_receptor_search(bool enabled) { m_soft_receptor_search = enabled; }
	void set_mixed_cutoff_search(bool enabled) { m_mixed_cutoff_search = enabled; }
	void set_contact_rotor_score(bool enabled) {
		if (enabled && (m_sf_choice != SF_VINA || m_no_refine))
			throw vina_runtime_error("contact_rotor_score requires the Vina scoring function and explicit receptor refinement");
		m_contact_rotor_score = enabled;
	}
	void set_search_counters(bool enabled) { m_search_counters = enabled; }
	void set_pocket_start_radius(double radius);
	//void set_ligand(OpenBabel::OBMol* mol);
	//void set_ligand(std::vector<OpenBabel::OBMol*> mol);
	void set_vina_weights(double weight_gauss1=-0.035579, double weight_gauss2=-0.005156,
						       double weight_repulsion=0.840245, double weight_hydrophobic=-0.035069,
						       double weight_hydrogen=-0.587439, double weight_glue=50,
						       double weight_rot=0.05846);
	void set_vinardo_weights(double weight_gauss1=-0.045,
							       double weight_repulsion=0.8, double weight_hydrophobic=-0.035,
							       double weight_hydrogen=-0.600, double weight_glue=50,
							       double weight_rot=0.05846);
	void set_ad4_weights(double weight_ad4_vdw=0.1662, double weight_ad4_hb=0.1209,
						      double weight_ad4_elec=0.1406, double weight_ad4_dsolv=0.1322,
						      double weight_glue=50, double weight_ad4_rot=0.2983);
	std::vector<double> grid_dimensions_from_ligand(double buffer_size=4);
	// Deferred setup preserves the box and receptor, then builds the maps at use.
	void compute_vina_maps(double center_x, double center_y, double center_z,
								  double size_x, double size_y, double size_z,
								  double granularity=0.5, bool force_even_voxels=false,
								  bool defer_population=false);
	// Prepare one shared box for batch ligands.
	void initialize_vina_batch_maps(double center_x, double center_y, double center_z,
	                               double size_x, double size_y, double size_z, double granularity);
	void load_maps(std::string maps);
	void randomize(const int max_steps=10000);
	std::vector<double> score();
	std::vector<double> score(double intramolecular_energy);
	std::vector<double> optimize(const int max_steps=0);
	std::vector<double> optimize(int max_steps, int max_evals, bool rigid_metric = false);
	// An optional Late cap before explicit refinement leaves the search and output caps intact.
	void global_search(const int exhaustiveness=8, const int n_poses=20, const double min_rmsd=1.0, const int max_evals=0,
	                   boost::optional<int> refine_candidates=boost::none);
	std::string get_poses(int how_many=9, double energy_range=3.0);
	std::vector< std::vector<double> > get_poses_coordinates(int how_many=9, double energy_range=3.0);
	std::vector< std::vector<double> > get_poses_energies(int how_many=9, double energy_range=3.0);
	void write_pose(const std::string &output_name, const std::string &remark = std::string());
	void write_poses(const std::string &output_name, int how_many=9, double energy_range=3.0);
	void write_maps(const std::string& map_prefix="receptor", const std::string& gpf_filename="NULL",
					    const std::string& fld_filename="NULL", const std::string& receptor_filename="NULL");
	void show_score(const std::vector<double> energies);

private:
	void set_single_ligand_model(model ligand);
	// model and poses
	model m_receptor;
	model m_model;
	output_container m_poses;
	bool m_search_completed = false;
	void invalidate_poses() { m_poses.clear(); m_search_completed = false; }
	fl m_output_min_rmsd = 0; // Deferred only when eligibility depends on a reranked native energy window.
	//OpenBabel::OBMol m_mol;
	bool m_receptor_initialized;
	bool m_ligand_initialized;
	// scoring function
	scoring_function_choice m_sf_choice;
	flv m_weights;
	std::shared_ptr<ScoringFunction> m_scoring_function;
	precalculate_byatom m_precalculated_byatom;
	double m_rerank_gauss2_scale = 1.0;
	std::shared_ptr<ScoringFunction> m_rerank_scoring_function;
	// maps
	cache m_grid;
	boost::optional<cache> m_mixed_search_grid;
	ad4cache m_ad4grid;
	non_cache m_non_cache;
	bool m_deferred_map_population = false;
	bool m_map_initialized;
	bool m_use_explicit_receptor = false;
	// global search
	int m_cpu;
	int m_seed;
	bool m_late_pose_selection = false;
	bool m_use_min_rmsd_for_merge = false;
	bool m_armijo_refinement = false;
	bool m_geometry_refinement = false;
	bool m_mixed_reorientation = false;
	bool m_centered_rotation = false;
	bool m_separate_hunt_chain = false;
	bool m_search_counters = false;
	bool m_contact_rotor_score = false;
	bool m_mixed_cutoff_search = false;
	bool m_prepared_torsion_starts = false;
	bool m_balanced_torsion = false;
	bool m_cubic_search_grid = false;
	bool m_soft_receptor_search = false;
	double m_pocket_start_radius = 0;
	enum class template_stage { rank, select, refine, contact_rank };
	template_stage m_template_stage = template_stage::rank;
	std::string m_pose_template_path;
	std::string m_search_anchor_path, m_search_amide_path, m_rank_amide_path;
	enum class anchor_trajectories { all, even, odd };
	anchor_trajectories m_anchor_trajectories = anchor_trajectories::all;
	enum class polar_mode { off, distance, directional, penalty };
	polar_mode m_polar_geometry = polar_mode::off;
	enum class burial_mode { off, burial, geometry, compensated };
	burial_mode m_polar_burial = burial_mode::off;
	enum class net_polar_mode { off, native, direction, burial, combined, water, satisfied, satisfied_water };
	net_polar_mode m_net_polar = net_polar_mode::off;
	vecv m_net_polar_waters;
	enum class archive_mode { off, accepted, support };
	archive_mode m_trajectory_archive = archive_mode::off;
	// others
	int m_verbosity;
	bool m_no_refine;
	std::function<void(double)>* m_progress_callback;

	std::string vina_remarks(output_type& pose, fl lb, fl ub);
	const char* pose_template_stage_name() const;
	const char* polar_geometry_name() const;
	const char* trajectory_archive_name() const;
	fl polar_rank_score(const output_type& pose) const;
	fl output_energy_reference() const;
	void write_pose_models(std::ostream& out, int how_many, double energy_range);
	bool current_pose_in_grid() const;
	output_container remove_redundant(const output_container& in, fl min_rmsd);
	void rerank_poses(output_container& poses, const vec& authentic_v);
	void rank_template_poses(output_container& poses, pose_template& reference, sz& rank_evals);
	void rank_polar_poses(output_container& poses);
	void rank_polar_burial_poses(output_container& poses);
	std::array<fl, 3> polar_burial_deltas(fl burial, fl geometry, fl uncompensated) const;
	bool supports_polar_burial() const;
	const char* polar_burial_name() const;
	const char* net_polar_name() const;
	bool uses_net_polar_satisfaction() const;
	void rank_net_polar_poses(output_container& poses);
	std::array<fl, 4> net_polar_deltas(const net_polar_terms& terms, fl normalization) const;
	non_cache fixed_pose_receptor();

	void reset_ligand_model(atom_type::t atom_typing);
	void finish_ligand_setup(atom_type::t atom_typing);
	void ensure_vina_maps();
	void build_vina_maps(const vec& center, const vec& span, fl granularity,
	                     bool force_even_voxels, bool defer_population);
	void set_forcefield();
	const igrid& vina_receptor_interactions();
	std::vector<double> optimize(output_type& out, int max_steps, int max_evals, bool rigid_metric);
	int generate_seed(const int seed=0);
};

#endif
