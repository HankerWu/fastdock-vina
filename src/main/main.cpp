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

// FastDock: prepare affinity maps on demand and validate search box options.

#include <iostream>
#include <set>
#include <string>
#include <vector> // ligand paths
#include <exception>
#include "parse_error.h"
#include <boost/program_options.hpp>
#include "vina.h"
#include "utils.h"
#include "scoring_function.h"
#include <unordered_map>
#include <boost/filesystem.hpp>

struct usage_error : public std::runtime_error {
	usage_error(const std::string& message) : std::runtime_error(message) {}
};

int main(int argc, char* argv[]) {
	using namespace boost::program_options;
	const std::string git_version = VERSION;
	const std::string version_string = "AutoDock Vina " + git_version;
	const std::string error_message = "\n\n\
Please report bugs through the Issue Tracker on GitHub \n\
(https://github.com/ccsb-scripps/AutoDock-Vina/issues)., so\n\
that this problem can be resolved. The reproducibility of the\n\
error may be vital, so please remember to include the following in\n\
your problem report:\n\
* the EXACT error message,\n\
* your version of the program,\n\
* the type of computer system you are running it on,\n\
* all command line options,\n\
* configuration file (if used),\n\
* ligand file as PDBQT,\n\
* receptor file as PDBQT,\n\
* flexible side chains file as PDBQT (if used),\n\
* output file as PDBQT (if any),\n\
* input (if possible),\n\
* random seed the program used (this is printed when the program starts).\n\
\n\
Thank you!\n";

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

	try {
		std::string rigid_name;
		std::string flex_name;
		std::string config_name;
		std::string out_name;
		std::string out_dir;
		std::string out_maps;
		std::vector<std::string> ligand_names;
		std::vector<std::string> batch_ligand_names;
		std::string maps;
		std::string sf_name = "vina";
		std::string pose_selection = "native";
		std::string refinement_optimizer = "native";
		std::string rigid_reorientation = "native";
		bool separate_hunt_chain = false;
		bool use_min_rmsd_for_merge = false;
		std::string pose_coordinates, rigid_window_pool;
		unsigned rigid_window_samples = 2048;
		std::string ligand_bonds;
		bool search_counters = false;
		bool contact_rotor_score = false;
		double pocket_start_radius = 0;
		std::string pose_template_path;
		std::string search_anchor_path, search_amide_path, rank_amide_path;
		std::string search_anchor_trajectories = "all";
		std::string template_stage;
		std::string polar_geometry;
		std::string polar_burial;
		std::string net_polar;
		std::string net_polar_waters;
		std::string trajectory_archive;
		boost::optional<int> refine_candidates;
		double center_x;
		double center_y;
		double center_z;
		double size_x;
		double size_y;
		double size_z;
		int cpu = 0;
		int seed = 0;
		int exhaustiveness = 8;
		int max_evals = 0;
		int verbosity = 1;
		int num_modes = 9;
		double min_rmsd = 1.0;
		double energy_range = 3.0;
		double grid_spacing = 0.375;
		double interaction_cutoff = default_xs_interaction_cutoff;
		bool smooth_cutoff = false;
		bool mixed_cutoff_search = false;
		bool prepared_torsion_starts = false;
		bool balanced_torsion_search = false;
		bool cubic_search_grid = false;
		bool soft_receptor_search = false;
		double buffer_size = 4;
		double unbound_energy = NAN;

		// autodock4.2 weights
		double weight_ad4_vdw   = 0.1662;
		double weight_ad4_hb    = 0.1209;
		double weight_ad4_elec  = 0.1406;
		double weight_ad4_dsolv = 0.1322;
		double weight_ad4_rot   = 0.2983;

		// vina weights
		double weight_gauss1      = -0.035579;
		double weight_gauss2      = -0.005156;
		double rerank_gauss2_scale = 1.0;
		double weight_repulsion   =  0.840245;
		double weight_hydrophobic = -0.035069;
		double weight_hydrogen    = -0.587439;
		double weight_rot         =  0.05846;

		// vinardo weights
		double weight_vinardo_gauss1 = -0.045;
		double weight_vinardo_repulsion = 0.8;
		double weight_vinardo_hydrophobic = -0.035;
		double weight_vinardo_hydrogen = -0.600;
		double weight_vinardo_rot = 0.05846;

		// macrocycle closure
		double weight_glue        = 50.000000; // linear attraction

		bool score_only = false;
		bool local_only = false;
		bool no_refine = false;
		bool force_even_voxels = false;
		bool randomize_only = false;
		bool help = false;
		bool help_advanced = false;
		bool version = false; // FIXME
		bool autobox = false;
		variables_map vm;

		positional_options_description positional; // remains empty

		options_description inputs("Input");
		inputs.add_options()
			("receptor", value<std::string>(&rigid_name), "rigid part of the receptor (PDBQT)")
			("flex", value<std::string>(&flex_name), "flexible side chains, if any (PDBQT)")
			("ligand", value< std::vector<std::string> >(&ligand_names)->multitoken(), "ligand (PDBQT)")
			("batch", value< std::vector<std::string> >(&batch_ligand_names)->multitoken(), "batch directory or ligands (PDBQT)")
			("scoring", value<std::string>(&sf_name)->default_value(sf_name), "scoring function (ad4, vina or vinardo)")
		;
		options_description search_area("Search space");
		search_area.add_options()
			("maps", value<std::string>(&maps), "affinity maps for the autodock4.2 (ad4) or vina scoring function")
			("center_x", value<double>(&center_x), "X coordinate of the center (Angstrom)")
			("center_y", value<double>(&center_y), "Y coordinate of the center (Angstrom)")
			("center_z", value<double>(&center_z), "Z coordinate of the center (Angstrom)")
			("size_x", value<double>(&size_x), "size in the X dimension (Angstrom)")
			("size_y", value<double>(&size_y), "size in the Y dimension (Angstrom)")
			("size_z", value<double>(&size_z), "size in the Z dimension (Angstrom)")
			("autobox", bool_switch(&autobox), "set maps dimensions based on input ligand(s) (for --score_only and --local_only)")
		;
		//options_description outputs("Output prefixes (optional - by default, input names are stripped of .pdbqt\nare used as prefixes. _001.pdbqt, _002.pdbqt, etc. are appended to the prefixes to produce the output names");
		options_description outputs("Output (optional)");
		outputs.add_options()
			("out", value<std::string>(&out_name), "output models (PDBQT), the default is chosen based on the ligand file name")
			("dir", value<std::string>(&out_dir), "output directory for batch mode")
			("write_maps", value<std::string>(&out_maps), "output filename (directory + prefix name) for maps. Option --force_even_voxels may be needed to comply with .map format")
		;
		options_description advanced("Advanced options (see the manual)");
		advanced.add_options()
			("score_only",     bool_switch(&score_only),     "score only - single-system receptor input defaults to autobox")
			("ligand_bonds", value<std::string>(&ligand_bonds),
			 "trusted complete ligand covalent graph as undirected PDBQT serial pairs")
			("rigid_window_pool", value<std::string>(&rigid_window_pool),
			 "integrate rigid-body fluctuations around each MODEL of a native pose pool (score_only)")
			("rigid_window_samples", value<unsigned>(&rigid_window_samples),
			 "Sobol samples per scale for rigid_window_pool (power of two, default 2048)")
			("contact_rotor_score", bool_switch(&contact_rotor_score),
			 "report the experimental contact-dependent rotor score alongside each native affinity")
			("pose_coordinates", value<std::string>(&pose_coordinates),
			 "fixed score only: use coordinates from one pose with the topology prepared from --ligand")
			("local_only",     bool_switch(&local_only),     "do local search only")
			("unbound_energy", value<double>(&unbound_energy)->default_value(unbound_energy), "Explicitly set the Unbound System's Energy for --score_only jobs")
			("no_refine", bool_switch(&no_refine),  "when --receptor is provided, do not use explicit receptor atoms (instead of precalculated grids) for: (1) local optimization and scoring after docking, (2) --local_only jobs, and (3) --score_only jobs")
			("refinement_optimizer", value<std::string>(&refinement_optimizer)->default_value(refinement_optimizer),
			 "native, armijo, geometry (Cartesian preconditioning), or rigid_metric (centroid preconditioning, local_only)")
			("rigid_reorientation", value<std::string>(&rigid_reorientation)->default_value(rigid_reorientation),
			 "native, mixed local/global, or centered local rotation for single-ligand rigid Vina docking")
			("separate_hunt_chain", bool_switch(&separate_hunt_chain),
			 "keep the Monte Carlo state on the hunt objective and refine only the saved copy")
			("pocket_start_radius", value<double>(&pocket_start_radius),
			 "place half the trajectories' initial ligand heavy centroids within this radius (Angstrom) of the search box center")
			("search_counters", bool_switch(&search_counters),
			 "print per-trajectory gradient counts and explicit refinement work")
			("force_even_voxels", bool_switch(&force_even_voxels),  "calculated grid maps will have an even number of voxels (intervals) in each dimension (odd number of grid points)")
			("randomize_only", bool_switch(&randomize_only), "randomize input, attempting to avoid clashes")

			("weight_gauss1", value<double>(&weight_gauss1)->default_value(weight_gauss1),                "gauss_1 weight")
			("weight_gauss2", value<double>(&weight_gauss2)->default_value(weight_gauss2),                "gauss_2 weight")
			("rerank_gauss2_scale", value<double>(&rerank_gauss2_scale)->default_value(rerank_gauss2_scale), "scale gauss_2 for FastDock Vina Pro final pose ranking only (0 to 1; 1 uses normal Vina)")
			("weight_repulsion", value<double>(&weight_repulsion)->default_value(weight_repulsion),       "repulsion weight")
			("weight_hydrophobic", value<double>(&weight_hydrophobic)->default_value(weight_hydrophobic), "hydrophobic weight")
			("weight_hydrogen", value<double>(&weight_hydrogen)->default_value(weight_hydrogen),          "Hydrogen bond weight")
			("weight_rot", value<double>(&weight_rot)->default_value(weight_rot),                         "N_rot weight")

			("weight_vinardo_gauss1", value<double>(&weight_vinardo_gauss1)->default_value(weight_vinardo_gauss1), "Vinardo gauss_1 weight")
			("weight_vinardo_repulsion", value<double>(&weight_vinardo_repulsion)->default_value(weight_vinardo_repulsion), "Vinardo repulsion weight")
			("weight_vinardo_hydrophobic", value<double>(&weight_vinardo_hydrophobic)->default_value(weight_vinardo_hydrophobic), "Vinardo hydrophobic weight")
			("weight_vinardo_hydrogen", value<double>(&weight_vinardo_hydrogen)->default_value(weight_vinardo_hydrogen), "Vinardo Hydrogen bond weight")
			("weight_vinardo_rot", value<double>(&weight_vinardo_rot)->default_value(weight_vinardo_rot), "Vinardo N_rot weight")

			("weight_ad4_vdw", value<double>(&weight_ad4_vdw)->default_value(weight_ad4_vdw), "ad4_vdw weight")
			("weight_ad4_hb", value<double>(&weight_ad4_hb)->default_value(weight_ad4_hb), "ad4_hb weight")
			("weight_ad4_elec", value<double>(&weight_ad4_elec)->default_value(weight_ad4_elec), "ad4_elec weight")
			("weight_ad4_dsolv", value<double>(&weight_ad4_dsolv)->default_value(weight_ad4_dsolv), "ad4_dsolv weight")
			("weight_ad4_rot", value<double>(&weight_ad4_rot)->default_value(weight_ad4_rot), "ad4_rot weight")

			("weight_glue", value<double>(&weight_glue)->default_value(weight_glue),                      "macrocycle glue weight")
		;
		options_description misc("Misc (optional)");
		misc.add_options()
			("cpu", value<int>(&cpu)->default_value(0), "the number of CPUs to use (the default is to try to detect the number of CPUs or, failing that, use 1)")
			("seed", value<int>(&seed)->default_value(0), "explicit random seed")
			("exhaustiveness", value<int>(&exhaustiveness)->default_value(8), "exhaustiveness of the global search (roughly proportional to time): 1+")
			("max_evals", value<int>(&max_evals)->default_value(0), "evaluation budget per MC run, or hard gradient limit for local_only with an Armijo optimizer (0 uses heuristics)")
			("num_modes", value<int>(&num_modes)->default_value(9), "maximum number of binding modes to generate")
			("pose_selection", value<std::string>(&pose_selection)->default_value(pose_selection),
			 "native search merge or late selection using explicit Vina scores before refinement")
			("use_min_rmsd_for_merge", bool_switch(&use_min_rmsd_for_merge),
			 "use min_rmsd for the native cross-trajectory merge instead of 2 Angstrom")
			("refine_candidates", value<int>(),
			 "maximum Late candidates sent to explicit refinement (at least num_modes; default num_modes)")
			("search_anchor", value<std::string>(&search_anchor_path),
			 "guide grid search (rows: ligand serial x y z [shell radius]; OR separates alternative assignments); native Late refinement and scoring")
			("search_amide", value<std::string>(&search_amide_path),
			 "guide all search trajectories with amide planarity (one O C N X [X2] barrier row per active C-N; native Late refinement and scoring)")
			("rank_amide", value<std::string>(&rank_amide_path),
			 "rank refined Late poses by native affinity plus amide strain (O C N X [X2] barrier); energy_range uses native affinity")
			("search_anchor_trajectories", value<std::string>(&search_anchor_trajectories),
			 "trajectories receiving search_anchor: all (default), even or odd (1-based)")
			("pose_template", value<std::string>(&pose_template_path),
			 "aligned external template with element, x, y and z columns (Angstrom)")
			("template_stage", value<std::string>(&template_stage),
			 "template guidance stage: rank, select, refine (element overlap), or contact_rank (receptor contact overlap)")
			("polar_geometry", value<std::string>(&polar_geometry),
			 "final Late Vina ranking by polar distance, directional geometry or angular penalty (default off)")
			("polar_burial", value<std::string>(&polar_burial),
			 "final Late Vina ranking by burial, geometry or compensated penalty, or fixed coordinates with score_only (default off)")
			("net_polar", value<std::string>(&net_polar),
			 "final Late Vina ranking or fixed-coordinate diagnostics (native|direction|burial|combined|water|satisfied|satisfied_water; default off)")
			("net_polar_waters", value<std::string>(&net_polar_waters),
			 "independent aligned water oxygens, one x y z row per water (Angstrom)")
			("trajectory_archive", value<std::string>(&trajectory_archive),
			 "add accepted Monte Carlo minima to Late Vina, ranked by native score or trajectory support (accepted|support, default off)")
			("min_rmsd", value<double>(&min_rmsd)->default_value(1.0), "minimum RMSD between output poses")
			("energy_range", value<double>(&energy_range)->default_value(3.0), "maximum energy difference between the best binding mode and the worst one displayed (kcal/mol)")
			("spacing", value<double>(&grid_spacing)->default_value(0.375), "grid spacing (Angstrom)")
			("prepared_torsion_starts", bool_switch(&prepared_torsion_starts),
			 "start even trajectories at input ligand torsions with a random rigid placement")
			("balanced_torsion_search", bool_switch(&balanced_torsion_search),
			 "distribute each ligand torsion proposal between both sides to minimize heavy-atom motion about the bond")
			("cubic_search_grid", bool_switch(&cubic_search_grid),
			 "use tricubic receptor-map interpolation in Late Vina search")
			("soft_receptor_search", bool_switch(&soft_receptor_search),
			 "halve receptor repulsion on even search trajectories only; keep configured Late postprocessing")
			("mixed_cutoff_search", bool_switch(&mixed_cutoff_search),
			 "search alternating trajectories with the smooth cutoff and select pooled poses with native Late Vina")
			("smooth_cutoff", bool_switch(&smooth_cutoff),
			 "taper ordinary Vina pair energies over the final 1 Angstrom of the cutoff")
			("interaction_cutoff", value<double>(&interaction_cutoff)->default_value(default_xs_interaction_cutoff),
			 "Vina/Vinardo interaction cutoff (Angstrom); custom values require direct receptor input")
			("verbosity", value<int>(&verbosity)->default_value(1), "verbosity (0=no output, 1=normal, 2=verbose)")
		;
		options_description config("Configuration file (optional)");
		config.add_options()
			("config", value<std::string>(&config_name), "the above options can be put here")
		;
		options_description info("Information (optional)");
		info.add_options()
			("help",          bool_switch(&help), "display usage summary")
			("help_advanced", bool_switch(&help_advanced), "display usage summary with advanced options")
			("version",       bool_switch(&version), "display program version")
		;
		options_description desc, desc_config, desc_simple;
		desc       .add(inputs).add(search_area).add(outputs).add(advanced).add(misc).add(config).add(info);
		desc_config.add(inputs).add(search_area).add(outputs).add(advanced).add(misc);
		desc_simple.add(inputs).add(search_area).add(outputs).add(misc).add(config).add(info);

		std::cout << version_string << '\n';
		try {
			//store(parse_command_line(argc, argv, desc, command_line_style::default_style ^ command_line_style::allow_guessing), vm);
			store(command_line_parser(argc, argv)
				.options(desc)
				.style(command_line_style::default_style ^ command_line_style::allow_guessing)
				.positional(positional)
				.run(),
				vm);
			notify(vm);
		} catch(boost::program_options::error& e) {
			std::cerr << "Command line parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
			return 1;
		}

		if (vm.count("config")) {
			try {
				path name = make_path(config_name);
				ifile config_stream(name);
				store(parse_config_file(config_stream, desc_config), vm);
				notify(vm);
			}
			catch(boost::program_options::error& e) {
				std::cerr << "Configuration file parse error: " << e.what() << '\n' << "\nCorrect usage:\n" << desc_simple << '\n';
				return 1;
			}
		}

		if (help) {
			std::cout << desc_simple << '\n';
			return 0;
		}

		if (help_advanced) {
			std::cout << desc << '\n';
			return 0;
		}

		if (version) {
			return 0;
		}

		if (verbosity > 0) {
			std::cout << cite_message << '\n';
		}

		const bool regular_rigid_vina = sf_name == "vina" && !rigid_name.empty() &&
		    !vm.count("flex") && !vm.count("maps") && !vm.count("write_maps") &&
		    ligand_names.size() <= 1 && !no_refine && !score_only && !local_only &&
		    !randomize_only && !autobox && vm["unbound_energy"].defaulted();
		if (rerank_gauss2_scale >= 0 && rerank_gauss2_scale < 1 && !regular_rigid_vina)
			throw usage_error("rerank_gauss2_scale below 1 requires regular Vina docking with a rigid receptor, one ligand per run, and explicit receptor refinement");

		if ((cubic_search_grid || mixed_cutoff_search || soft_receptor_search || vm.count("search_anchor") || vm.count("search_amide") || vm.count("rank_amide")) && vm["pose_selection"].defaulted()) pose_selection = "late";
		if (vm.count("search_anchor_trajectories") && !vm.count("search_anchor"))
			throw usage_error("search_anchor_trajectories requires search_anchor");
		if (vm.count("search_anchor") && (search_anchor_path.empty() || !regular_rigid_vina ||
		    ligand_names.size() != 1 || vm.count("batch") || pose_selection != "late"))
			throw usage_error("search_anchor requires regular rigid single-ligand Late Vina docking");
		if (vm.count("search_amide") && (search_amide_path.empty() || !regular_rigid_vina ||
		    ligand_names.size() != 1 || vm.count("batch") || pose_selection != "late"))
			throw usage_error("search_amide requires regular rigid single-ligand Late Vina docking");
		if (vm.count("rank_amide") && (rank_amide_path.empty() || !regular_rigid_vina ||
		    ligand_names.size() != 1 || vm.count("batch") || pose_selection != "late"))
			throw usage_error("rank_amide requires regular rigid single-ligand Late Vina docking");
		if (cubic_search_grid && (!regular_rigid_vina || pose_selection != "late" || soft_receptor_search || mixed_cutoff_search))
			throw usage_error("cubic_search_grid requires rigid Late Vina with a single receptor search grid");
		if (soft_receptor_search && (!regular_rigid_vina || ligand_names.size() != 1 || vm.count("batch") ||
		    pose_selection != "late" || exhaustiveness < 2))
			throw usage_error("soft_receptor_search requires rigid single-ligand Late Vina docking and at least two trajectories");
		if (mixed_cutoff_search && (!regular_rigid_vina || pose_selection != "late" || smooth_cutoff || exhaustiveness < 2))
			throw usage_error("mixed_cutoff_search requires regular rigid Vina docking, native Late scoring and exhaustiveness >= 2");
		if (pose_selection != "native" && pose_selection != "late")
			throw usage_error("pose_selection must be native or late");
		if (pose_selection == "late" && !regular_rigid_vina)
			throw usage_error("late pose selection requires regular Vina docking with a rigid receptor, one ligand per search, and explicit receptor refinement");
		if (use_min_rmsd_for_merge && (!regular_rigid_vina || pose_selection != "native" ||
		    !std::isfinite(min_rmsd) || min_rmsd < 0))
			throw usage_error("use_min_rmsd_for_merge requires rigid single-ligand native Vina docking and finite non-negative min_rmsd");
		if (rigid_reorientation != "native" && rigid_reorientation != "mixed" && rigid_reorientation != "centered")
			throw usage_error("rigid_reorientation must be native, mixed or centered");
		if (rigid_reorientation != "native" && !regular_rigid_vina)
			throw usage_error(rigid_reorientation + " rigid_reorientation requires regular Vina docking with a rigid receptor, one ligand per search, and explicit receptor refinement");
		if (separate_hunt_chain && !regular_rigid_vina)
			throw usage_error("separate_hunt_chain requires regular Vina docking with a rigid receptor, one ligand per search, and explicit receptor refinement");
		if (vm.count("pocket_start_radius")) {
			if (!std::isfinite(pocket_start_radius) || pocket_start_radius <= 0)
				throw usage_error("pocket_start_radius must be finite and positive");
			if (score_only || local_only || randomize_only || ligand_names.size() != 1 ||
			    rigid_name.empty() || vm.count("batch") || vm.count("flex") || vm.count("maps") ||
			    vm.count("write_maps") || (sf_name != "vina" && sf_name != "vinardo"))
				throw usage_error("pocket_start_radius requires global Vina or Vinardo docking with one ligand and a rigid receptor");
		}
		if (prepared_torsion_starts && (randomize_only || score_only || local_only || exhaustiveness < 2))
			throw usage_error("prepared_torsion_starts requires global docking with at least two trajectories");
		if (balanced_torsion_search && (randomize_only || score_only || local_only))
			throw usage_error("balanced_torsion_search requires global docking");
		if (search_counters && (randomize_only || score_only || local_only))
			throw usage_error("search_counters requires global docking");
		const bool rigid_metric = refinement_optimizer == "rigid_metric";
		if (refinement_optimizer != "native" && refinement_optimizer != "armijo" &&
		    refinement_optimizer != "geometry" && !rigid_metric)
			throw usage_error("refinement_optimizer must be native, armijo, geometry, or rigid_metric");
		const bool rigid_local_vina = local_only && !score_only && !randomize_only && sf_name == "vina" &&
		    !rigid_name.empty() && ligand_names.size() == 1 && !no_refine && !vm.count("batch") &&
		    !vm.count("flex") && !vm.count("maps") && !vm.count("write_maps");
		if (rigid_metric && !rigid_local_vina)
			throw usage_error("rigid_metric requires local_only with one ligand and an explicit rigid Vina receptor");
		if (refinement_optimizer != "native" && !regular_rigid_vina && !rigid_local_vina)
			throw usage_error("Armijo refinement requires Vina scoring, one ligand per search, and an explicit rigid receptor");
		if (rigid_local_vina && refinement_optimizer != "native" && max_evals < 0)
			throw usage_error("Local refinement requires non-negative max_evals");
		if (vm.count("pose_template") != vm.count("template_stage"))
			throw usage_error("Specify pose_template and template_stage together");
		if (vm.count("pose_template")) {
			const bool contact = template_stage == "contact_rank";
			if (!regular_rigid_vina || ligand_names.size() != 1 || vm.count("batch") || rerank_gauss2_scale != 1 ||
			    (!contact && (pose_selection != "late" || refinement_optimizer != "native")))
				throw usage_error(contact ? "contact_rank requires one rigid Vina receptor, one ligand, explicit refinement, and no Gaussian rescaling" :
				                  "pose_template requires Late Vina docking with one rigid receptor, one ligand, native refinement, and no Gaussian rescaling");
		}
		if (vm.count("polar_geometry") && (!regular_rigid_vina || ligand_names.size() != 1 || vm.count("batch") ||
		    pose_selection != "late" || rerank_gauss2_scale != 1 || vm.count("pose_template")))
			throw usage_error("polar_geometry requires rigid single-ligand Late Vina docking without template or Gaussian rescaling");
		if (vm.count("ligand_bonds") &&
		    (ligand_names.size() != 1 || rigid_name.empty() || vm.count("batch") || vm.count("flex") ||
		     vm.count("maps") || vm.count("write_maps") || (sf_name != "vina" && sf_name != "vinardo")))
			throw usage_error("ligand_bonds requires one prepared Vina or Vinardo ligand and a rigid receptor");
		if (contact_rotor_score && (sf_name != "vina" || rigid_name.empty() || ligand_names.size() != 1 ||
		    vm.count("flex") || vm.count("maps") || vm.count("batch") || vm.count("write_maps") ||
		    vm.count("rigid_window_pool") || randomize_only || local_only || no_refine))
			throw usage_error("contact_rotor_score requires one rigid-receptor Vina docking or explicit score_only calculation");
		if (vm.count("pose_coordinates") &&
		    (!score_only || local_only || randomize_only || ligand_names.size() != 1 ||
		     rigid_name.empty() || vm.count("batch") || vm.count("flex") || vm.count("maps") ||
		     vm.count("write_maps") || no_refine || (sf_name != "vina" && sf_name != "vinardo")))
			throw usage_error("pose_coordinates requires score_only with one prepared Vina or Vinardo ligand and an explicit rigid receptor");
		if (vm.count("rigid_window_samples") && !vm.count("rigid_window_pool"))
			throw usage_error("rigid_window_samples requires rigid_window_pool");
		if (vm.count("rigid_window_pool") && (!score_only || local_only || randomize_only || no_refine ||
		    ligand_names.size()!=1 || rigid_name.empty() || sf_name!="vina" || vm.count("maps") ||
		    vm.count("batch") || vm.count("flex") || vm.count("write_maps") || (vm.count("unbound_energy") && !vm["unbound_energy"].defaulted()) ||
		    vm.count("pose_coordinates") || vm.count("net_polar") || vm.count("polar_burial")))
			throw usage_error("rigid_window_pool requires score_only, one prepared ligand, and an explicit rigid Vina receptor");
		const bool polar_burial_score_only = vm.count("polar_burial") && score_only;
		const bool net_polar_score_only = vm.count("net_polar") && score_only;
		if (vm.count("net_polar_waters") && (!vm.count("net_polar") || net_polar_waters.empty()))
			throw usage_error("net_polar_waters requires net_polar and a nonempty filename");
		if (vm.count("polar_burial") && vm.count("net_polar"))
			throw usage_error("net_polar and polar_burial select separate final scoring modes");
		if (vm.count("polar_burial") || vm.count("net_polar")) {
			const bool single_rigid_vina = sf_name == "vina" && !rigid_name.empty() && ligand_names.size() == 1 &&
			    !vm.count("batch") && !vm.count("flex") && !vm.count("maps") && !vm.count("write_maps") &&
			    !no_refine && !local_only && !randomize_only && vm["unbound_energy"].defaulted();
			if (!single_rigid_vina || (!score_only && (!regular_rigid_vina || pose_selection != "late")) ||
			    (score_only && vm.count("refine_candidates")) ||
			    rerank_gauss2_scale != 1 || refinement_optimizer != "native" || vm.count("pose_template") ||
			    vm.count("polar_geometry") || vm.count("trajectory_archive"))
				throw usage_error(std::string(vm.count("net_polar") ? "net_polar" : "polar_burial") +
			                  " requires one ligand, an explicit rigid receptor, and either Late Vina docking with native refinement or score_only");
		}
		if (vm.count("trajectory_archive") && (trajectory_archive != "accepted" && trajectory_archive != "support"))
			throw usage_error("trajectory_archive must be accepted or support");
		if (vm.count("trajectory_archive") && (!regular_rigid_vina || ligand_names.size() != 1 || vm.count("batch") ||
		    pose_selection != "late" || rerank_gauss2_scale != 1 || refinement_optimizer != "native" ||
		    vm.count("pose_template") || vm.count("polar_geometry")))
			throw usage_error("trajectory_archive requires rigid single-ligand Late Vina docking with native refinement and no template, polar geometry, or Gaussian rescaling");
		if (vm.count("refine_candidates")) {
			refine_candidates = vm["refine_candidates"].as<int>();
			if (*refine_candidates <= 0 || *refine_candidates < num_modes)
				throw usage_error("refine_candidates must be positive and at least num_modes");
			if (!regular_rigid_vina || pose_selection != "late" || rerank_gauss2_scale != 1 ||
			    refinement_optimizer != "native" || vm.count("pose_template") || vm.count("polar_geometry"))
				throw usage_error("refine_candidates requires Late Vina docking with native refinement and no template, polar geometry, or Gaussian rescaling");
		}

		if (vm.count("receptor") && vm.count("maps")) {
			std::cerr << "ERROR: Cannot specify both receptor and affinity maps at the same time, --flex argument is allowed with receptor or maps.\n";
			exit(EXIT_FAILURE);
		}

		if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
			if (!vm.count("receptor") && !vm.count("maps")) {
				std::cerr << desc_simple << "ERROR: The receptor or affinity maps must be specified.\n";
				exit(EXIT_FAILURE);
			}
		} else if (sf_name.compare("ad4") == 0) {
			if (vm.count("receptor")) {
				std::cerr << "ERROR: No receptor allowed, only --flex argument with the AD4 scoring function.\n";
				exit(EXIT_FAILURE);
			}
			if (!vm.count("maps")) {
				std::cerr << desc_simple << "\n\nERROR: Affinity maps are missing.\n";
				exit(EXIT_FAILURE);
			}
		} else {
			std::cerr << desc_simple << "Scoring function " << sf_name << " unknown.\n";
			exit(EXIT_FAILURE);
		}

		if (!vm.count("ligand") && !vm.count("batch")) {
			std::cerr << desc_simple << "\n\nERROR: Missing ligand(s).\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("ligand") && vm.count("batch")) {
			std::cerr << desc_simple << "\n\nERROR: Can't use both --ligand and --batch arguments simultaneously.\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("batch") && !vm.count("dir")) {
			std::cerr << desc_simple << "\n\nERROR: Need to specify an output directory for batch mode.\n";
			exit(EXIT_FAILURE);
		} else if (vm.count("dir")) {
			if (!is_directory(out_dir)) {
				std::cerr << "ERROR: Directory " << out_dir << " does not exist.\n";
				exit(EXIT_FAILURE);
			}
		} else if (vm.count("ligand") && vm.count("dir")) {
			std::cerr << "WARNING: In ligand mode, --dir argument is ignored.\n";
		}

		const bool ligand_autobox = vm.count("ligand") && (score_only || local_only);
		if (!vm.count("maps")) {
			const auto box_fields = vm.count("center_x") + vm.count("center_y") + vm.count("center_z")
			                      + vm.count("size_x") + vm.count("size_y") + vm.count("size_z");
			if (box_fields == 0 && ligand_autobox && score_only && !randomize_only)
				autobox = true;
			if (!(ligand_autobox && autobox) && box_fields != 6)
				throw usage_error("Provide all three center and all three size values for the search box");
		}

		if (!score_only) {
			if (!vm.count("out") && ligand_names.size() == 1) {
				out_name = default_output(ligand_names[0]);
				std::cout << "Output will be " << out_name << '\n';
			} else if (!vm.count("out") && ligand_names.size() >= 1) {
				std::cerr << desc_simple << "\n\nERROR: Output name must be defined when docking simultaneously multiple ligands.\n";
				exit(EXIT_FAILURE);
			}
		}

		if (verbosity > 0) {
			std::cout << "Scoring function : " << sf_name << "\n";
			if (vm.count("receptor"))
				std::cout << "Rigid receptor: " << rigid_name << "\n";
			if (vm.count("flex"))
				std::cout << "Flex receptor: " << flex_name << "\n";
			if (ligand_names.size() == 1) {
				std::cout << "Ligand: " << ligand_names[0] << "\n";
			} else if (ligand_names.size() > 1) {
				std::cout << "Ligands:\n";
				VINA_RANGE(i, 0, ligand_names.size()) {
					std::cout << "  - " << ligand_names[i] << "\n";
				}
			} else if (batch_ligand_names.size() > 1) {
				std::cout << "Ligands (batch mode): " << batch_ligand_names.size() << " molecules\n";
			}
			if (!vm.count("maps") && !autobox) {
				std::cout << "Grid center: X " << center_x << " Y " << center_y << " Z " << center_z << "\n";
				std::cout << "Grid size  : X " << size_x << " Y " << size_y << " Z " << size_z << "\n";
				std::cout << "Grid space : " << grid_spacing << "\n";
			} else if (autobox) {
				std::cout << "Grid center: ligand center (autobox)\n";
				std::cout << "Grid size  : ligand size + " << buffer_size << " A in each dimension (autobox)\n";
				std::cout << "Grid space : " << grid_spacing << "\n";
			}
			std::cout << "Exhaustiveness: " << exhaustiveness << "\n";
			std::cout << "CPU: " << cpu << "\n";
			if (!vm.count("seed"))
				std::cout << "Seed: " << seed << "\n";
			std::cout << "Verbosity: " << verbosity << "\n";
			std::cout << "\n";
		}

		Vina v(sf_name, cpu, seed, verbosity, no_refine);
		v.set_armijo_refinement(refinement_optimizer != "native", refinement_optimizer == "geometry");
		v.set_rigid_reorientation(rigid_reorientation == "mixed", rigid_reorientation == "centered");
		v.set_separate_hunt_chain(separate_hunt_chain);
		v.set_mixed_cutoff_search(mixed_cutoff_search);
		v.set_prepared_torsion_starts(prepared_torsion_starts);
		v.set_balanced_torsion_search(balanced_torsion_search);
		v.set_cubic_search_grid(cubic_search_grid);
		v.set_soft_receptor_search(soft_receptor_search);
		v.set_search_counters(search_counters);
		v.set_contact_rotor_score(contact_rotor_score);
		if (vm.count("pocket_start_radius"))
			v.set_pocket_start_radius(pocket_start_radius);
		v.set_late_pose_selection(pose_selection == "late");
		v.set_use_min_rmsd_for_merge(use_min_rmsd_for_merge);
		v.set_rerank_gauss2_scale(rerank_gauss2_scale);
		if (vm.count("search_anchor"))
			v.set_search_anchor(search_anchor_path, search_anchor_trajectories);
		if (vm.count("rank_amide"))
			v.set_rank_amide(rank_amide_path);
		if (vm.count("search_amide"))
			v.set_search_amide(search_amide_path);
		if (vm.count("pose_template"))
			v.set_pose_template(pose_template_path, template_stage);
		if (vm.count("polar_geometry"))
			v.set_polar_geometry(polar_geometry);
		if (vm.count("polar_burial"))
			v.set_polar_burial(polar_burial);
		if (vm.count("net_polar"))
			v.set_net_polar(net_polar, net_polar_waters);
		if (vm.count("trajectory_archive"))
			v.set_trajectory_archive(trajectory_archive);
		if (smooth_cutoff || !vm["interaction_cutoff"].defaulted()) {
			if (smooth_cutoff && (vm.count("maps") || vm.count("write_maps")))
				throw usage_error("map files cannot encode smooth_cutoff");
			if (interaction_cutoff != default_xs_interaction_cutoff && (vm.count("maps") || vm.count("write_maps")))
				throw usage_error("map files require the default interaction_cutoff");
			v.set_interaction_cutoff(interaction_cutoff, smooth_cutoff);
			if (smooth_cutoff && verbosity > 0)
				std::cout << "Vina smooth cutoff: " << interaction_cutoff - 1 << " to " << interaction_cutoff << " Angstrom\n";
		}

		// rigid_name variable can be ignored for AD4
		if (vm.count("receptor") || vm.count("flex"))
			v.set_receptor(rigid_name, flex_name);

		// Technically we don't have to initialize weights,
		// because they are initialized during the Vina object creation with the default weights
		// but we still do it in case the user decided to change them
		if (sf_name.compare("vina") == 0) {
			v.set_vina_weights(weight_gauss1, weight_gauss2, weight_repulsion,
							   weight_hydrophobic, weight_hydrogen, weight_glue, weight_rot);
		} else if (sf_name.compare("vinardo") == 0) {
			v.set_vinardo_weights(weight_vinardo_gauss1, weight_vinardo_repulsion,
								  weight_vinardo_hydrophobic, weight_vinardo_hydrogen, weight_glue, weight_vinardo_rot);
		} else {
			v.set_ad4_weights(weight_ad4_vdw, weight_ad4_hb, weight_ad4_elec,
							  weight_ad4_dsolv, weight_glue, weight_ad4_rot);
			v.load_maps(maps);

			// It works, but why would you do this?!
			if (vm.count("write_maps"))
				v.write_maps(out_maps);
		}

		const auto report_fixed_score = [&] {
			v.show_score(vm["unbound_energy"].defaulted() ? v.score() : v.score(unbound_energy));
		};
		if (vm.count("ligand")) {
			if (vm.count("ligand_bonds"))
				v.set_ligand_from_file_with_bonds(ligand_names.front(), ligand_bonds);
			else
				v.set_ligand_from_file(ligand_names);
			if (vm.count("rigid_window_pool")) {
				v.show_rigid_window_pool(rigid_window_pool,rigid_window_samples);
				return 0;
			}
			if (vm.count("pose_coordinates"))
				v.set_pose_coordinates_for_scoring(pose_coordinates);
			if (net_polar_score_only) {
				v.show_net_polar();
				return 0;
			}
			if (polar_burial_score_only) {
				v.show_polar_burial();
				return 0;
			}

			if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
				if (vm.count("maps")) {
					v.load_maps(maps);
				} else {
					const bool defer_maps = (score_only || randomize_only) && !vm.count("write_maps");

					if (ligand_autobox && autobox) {
						std::vector<double> dim = v.grid_dimensions_from_ligand(buffer_size);
						v.compute_vina_maps(dim[0], dim[1], dim[2], dim[3], dim[4], dim[5],
						                    grid_spacing, force_even_voxels, defer_maps);
					} else {
						v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z,
						                    grid_spacing, force_even_voxels, defer_maps);
					}

					if (vm.count("write_maps"))
						v.write_maps(out_maps);
				}
			}

			if (randomize_only) {
				v.randomize();
				v.write_pose(out_name);
			} else if (score_only) {
				report_fixed_score();
			} else if (local_only) {
				const auto energies = refinement_optimizer == "native" ? v.optimize() : v.optimize(0, max_evals, rigid_metric);
				v.write_pose(out_name);
				v.show_score(energies);
			} else {
				v.global_search(exhaustiveness, num_modes, min_rmsd, max_evals, refine_candidates);
				v.write_poses(out_name, num_modes, energy_range);
			}
		} else if (vm.count("batch")) {
			if (sf_name.compare("vina") == 0 || sf_name.compare("vinardo") == 0) {
				if (vm.count("maps")) {
					v.load_maps(maps);
				} else if (vm.count("write_maps")) {
					v.compute_vina_maps(center_x, center_y, center_z, size_x, size_y, size_z, grid_spacing);
					v.write_maps(out_maps);
				} else {
					v.initialize_vina_batch_maps(center_x, center_y, center_z, size_x, size_y, size_z, grid_spacing);
				}
			}

			// Iter all .pdbqt files in the batch input directory
			if (batch_ligand_names.size() == 1 && is_directory(batch_ligand_names[0])) {
				std::string in_dir = batch_ligand_names[0];
				batch_ligand_names.clear();
				for (const auto& entry : boost::filesystem::directory_iterator(in_dir)) {
					if (entry.path().extension() == ".pdbqt") {
						batch_ligand_names.push_back(entry.path().string());
					}
				}
			}

			std::set<std::string> repeated_names;
			std::set<std::string> raw_names;
			std::string name;
			VINA_RANGE(i, 0, batch_ligand_names.size()) {
				name = get_filename(batch_ligand_names[i]);
				if (raw_names.count(name)) {
					repeated_names.insert(name);
				}
				raw_names.insert(name);
			}
			std::unordered_map<std::string, int> instance_counter;

			std::size_t failed_ligand_parsing = 0;
			std::size_t failed_scores = 0;
			VINA_RANGE(i, 0, batch_ligand_names.size()) {
				name = get_filename(batch_ligand_names[i]);
				if (repeated_names.count(name)) {
					if (instance_counter.count(name)) {
						instance_counter[name] += 1;
					} else {
						instance_counter[name] = 1;
					}
					out_name = default_output(name, out_dir, instance_counter[name]);
				} else {
					out_name = default_output(name, out_dir);
				}

				try {
					v.set_ligand_from_file(batch_ligand_names[i]);
				}
				catch(pdbqt_parse_error& e) {
					std::cerr << e.what();
					std::cout << "Failed parsing " << batch_ligand_names[i] << ". Skipping it.\n";
					failed_ligand_parsing++;
					continue;
				}

				if (randomize_only) {
					v.randomize();
					v.write_pose(out_name);
				} else if (score_only) {
					std::cout << "Scoring ligand: " << batch_ligand_names[i] << '\n';
					try {
						report_fixed_score();
					} catch (const vina_runtime_error& error) {
						std::cerr << "Failed scoring " << batch_ligand_names[i] << ": " << error.what();
						++failed_scores;
					}
				} else if (local_only) {
					v.optimize();
					v.write_pose(out_name);
				} else {
					v.global_search(exhaustiveness, num_modes, min_rmsd, max_evals, refine_candidates);
					v.write_poses(out_name, num_modes, energy_range);
				}
			}
			if (repeated_names.size()) {
				std::cout << "Found " << repeated_names.size() << " repeated filenames in the input batch.\n";
				std::cout << "The corresponding output filenames are suffixed with _instance<n>_out.pdbqt\n";
			}
			if (failed_ligand_parsing) {
				std::cout << "Failed to parse " << failed_ligand_parsing << " ligands.\n";
			}
			if (failed_scores)
				std::cerr << "Failed to score " << failed_scores << " ligands.\n";
			if (score_only && (failed_scores || failed_ligand_parsing))
				return 1;
		}
	}

	catch(pdbqt_parse_error& e) {
		std::cerr << e.what();
		return 1;
	}
	catch(file_error& e) {
		std::cerr << "\n\nError: could not open \"" << e.name.string() << "\" for " << (e.in ? "reading" : "writing") << ".\n";
		return 1;
	}
	catch(boost::filesystem::filesystem_error& e) {
		std::cerr << "\n\nFile system error: " << e.what() << '\n';
		return 1;
	}
	catch(usage_error& e) {
		std::cerr << "\n\nUsage error: " << e.what() << ".\n";
		return 1;
	}
	catch(std::bad_alloc&) {
		std::cerr << "\n\nError: insufficient memory!\n";
		return 1;
	}

	// Errors that shouldn't happen:
	catch(std::exception& e) {
		std::cerr << "\n\nAn error occurred: " << e.what() << ". " << error_message;
		return 1;
	}
	catch(internal_error& e) {
		std::cerr << "\n\nAn internal error occurred in " << e.file << "(" << e.line << "). " << error_message;
		return 1;
	}
	catch(...) {
		std::cerr << "\n\nAn unknown error occurred. " << error_message;
		return 1;
	}
}
