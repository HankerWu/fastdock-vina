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

// FastDock: restore only mutable search state when reusing a worker's model.

#ifndef VINA_MODEL_H
#define VINA_MODEL_H

#include <memory>
#include <boost/optional.hpp> // for context

#include "file.h"
#include "tree.h"
#include "matrix.h"
#include "igrid.h"
#include "grid_dim.h"


struct interacting_pair {
	// Index the score table by XS type pair or, for charge-dependent AD4, by atom pair.
	sz table_index;
	sz a;
	sz b;
	interacting_pair(sz table_index_, sz a_, sz b_) : table_index(table_index_), a(a_), b(b_) {}
};

typedef std::vector<interacting_pair> interacting_pairs;

typedef std::pair<std::string, boost::optional<sz> > parsed_line;
typedef std::vector<parsed_line> context;

struct ligand : public flexible_body, atom_range {
	unsigned active_heavy_torsions = 0; // Two-sided PDBQT axes used in the XS score divisor.
	unsigned degrees_of_freedom; // can be different from the apparent number of rotatable bonds, because of the disabled torsions
	interacting_pairs pairs;
	context cont;
	ligand(const flexible_body& f, unsigned degrees_of_freedom_) : flexible_body(f), atom_range(0, 0), degrees_of_freedom(degrees_of_freedom_) {}
	void set_range();
};

struct residue : public main_branch {
	residue(const main_branch& m) : main_branch(m) {}
};

enum distance_type {DISTANCE_FIXED, DISTANCE_ROTOR, DISTANCE_VARIABLE};
typedef strictly_triangular_matrix<distance_type> distance_type_matrix;

struct cache; // forward declaration
struct szv_grid; // forward declaration
struct pdbqt_initializer; // forward declaration - only declared in parse_pdbqt.cpp
struct precalculate_byatom; // forward declaration
class pose_template;

fl eval_interacting_pairs(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, const bool closure_pairs=false);
fl eval_interacting_pairs_deriv(const precalculate_byatom& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces, const bool closure_pairs=false);

struct model {
	using context_line_filter = bool (*)(const std::string&);

public:
	// Had to move it from private to public to make it work. 
	// So we might have to fix that later
	model() : m_num_movable_atoms(0), m_atom_typing_used(atom_type::XS) {}
	model(atom_type::t atype) : m_num_movable_atoms(0), m_atom_typing_used(atype) {}

    atom get_atom(sz i) const { return atoms[i]; }
    ligand get_ligand(sz i) const { return ligands[i]; }
	vec get_coords(sz i) const { return coords[i]; }
	interacting_pairs get_other_pairs() const { return other_pairs; }
	interacting_pairs get_inter_pairs() const { return inter_pairs; }
	interacting_pairs get_glue_pairs() const { return glue_pairs; }

	void append(const model& m);
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }

	bool is_atom_in_ligand(sz a) const;
	bool is_movable_atom(sz a) const;
	std::vector<double> center() const;
	sz find_ligand(sz a) const;
	sz num_atoms() const { return atoms.size(); }
	sz num_movable_atoms() const { return m_num_movable_atoms; }
	sz num_internal_pairs() const;
	sz num_other_pairs() const { return other_pairs.size(); }
	sz num_ligands() const { return ligands.size(); }
	sz num_flex() const { return flex.size(); }
	sz ligand_degrees_of_freedom(sz ligand_number) const { return ligands[ligand_number].degrees_of_freedom; }
	sz ligand_longest_branch(sz ligand_number) const;
	sz ligand_length(sz ligand_number) const;

	szv get_movable_atom_types(atom_type::t atom_typing_used_) const;
	vecv get_ligand_coords() const;
	std::vector<double> get_ligand_coords();
	bool within_grid(const grid_dims& gd, fl margin) const;
	vecv get_grid_atom_movable_coords() const;
	vecv get_heavy_atom_movable_coords() const;
	vecv get_chemical_heavy_atom_movable_coords() const;
	// Set c and return the regularized inverse diagonal of its Cartesian RMS metric.
	fl max_heavy_atom_speed(const conf& base, const flv& direction);
	flv geometry_inverse_metric(const conf& c);
	std::vector<std::pair<sz, std::string>> chemical_heavy_atom_elements() const;
	fl template_overlap(pose_template& reference) const;
	conf_size get_size() const;
	conf get_initial_conf() const; // torsions = 0, orientations = identity, ligand positions = current

	void write_flex  (                  const path& name, const std::string& remark) const { write_context(flex_context, name, remark); }
	void write_ligand(sz ligand_number, const path& name, const std::string& remark) const { VINA_CHECK(ligand_number < ligands.size()); write_context(ligands[ligand_number].cont, name, remark); }
	void write_structure(std::ostream& out) const {
		write_structure_context(out, nullptr);
	}
	void write_structure(std::ostream& out, const std::string& remark) const {
		out << remark;
		write_structure(out);
	}
	void write_structure(ofile& out, std::vector<std::string>& remarks) const {
		VINA_FOR_IN(i, remarks) out << remarks[i];
		write_structure(out);
	}

	void write_structure(const path& name) const { ofile out(name); write_structure(out); }
	void write_model(std::ostream& out, sz model_number, const std::string& remark, context_line_filter omit_line = nullptr) const {
		out << "MODEL " << model_number << '\n';
		out << remark;
		write_structure_context(out, omit_line);
		out << "ENDMDL\n";
	}

	void set (const conf& c);

	void set_balanced_torsion(ligand_conf& c, sz ligand_number, sz torsion, fl angle);
	fl gyration_radius(sz ligand_number) const; // uses coords
	vec ligand_heavy_centroid(sz ligand_number) const;
	vec ligand_heavy_centroid(const ligand_conf& conf) const;

	const atomv& grid_atoms() const;
	// Indexed by movable receptor atom, before ligands are appended.
	std::vector<szv> rigid_bond_exclusions() const;
	const atom_base& movable_atom  (sz i) const { assert(i < m_num_movable_atoms); return  atoms[i]; }
	const vec&       movable_coords(sz i) const { assert(i < m_num_movable_atoms); return coords[i]; }

	const vec& atom_coords(const atom_index& i) const;
	fl distance_sqr_between(const atom_index& a, const atom_index& b) const;
	bool atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const; // there is an atom closer to both a and b then they are to each other and immobile relative to them

	distance_type distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const;

	// clean up
	fl evalo     (const precalculate_byatom& p,                  const vec& v           ) const;
	fl evali     (const precalculate_byatom& p,                  const vec& v           ) const;
	fl eval_inter(const precalculate_byatom& p,                  const vec& v           ) const;
	// Gradient uses the same layout as conf::set_incremented.
	fl eval_deriv(const precalculate_byatom& p, const igrid& ig, const vec& v, flv& gradient,
	              pose_template* reference = nullptr);
	fl eval_intramolecular(const precalculate_byatom& p, const igrid& ig, const vec& v);

	fl rmsd_lower_bound(const model& m) const; // uses coords
	fl rmsd_upper_bound(const model& m) const; // uses coords
	fl rmsd_ligands_upper_bound(const model& m) const; // uses coords

	void verify_bond_lengths() const;
	void about() const;
	void check_ligand_internal_pairs() const;
	void show_pairs() const;
	void show_atoms() const;
	void show_forces() const;
	void print_stuff(bool show_coords=true, bool show_internal=true, bool show_atoms=true, bool show_grid=true, bool show_about=true) const; // FIXME rm

	fl clash_penalty() const;

	const atom& get_atom(const atom_index& i) const { return i.in_grid ? (*m_grid_atoms)[i.i] : atoms[i.i]; }

private:
	friend struct parallel_mc;
	friend struct precalculate_byatom;
	friend struct cache;
	friend class anchored_search_grid;
	friend class amide_search_grid;
	friend class rotor_partition;
	friend struct non_cache;
	friend struct ad4cache;
	friend struct szv_grid;
	friend struct appender_info;
	friend struct pdbqt_initializer;
	friend class polar_geometry;
	friend class net_polar;

	void reset_search_state(const model& initial) {
		coords = initial.coords;
		minus_forces = initial.minus_forces;
		flex = initial.flex;
		VINA_FOR_IN(i, ligands)
			static_cast<flexible_body&>(ligands[i]) = initial.ligands[i];
	}

	void write_structure_context(std::ostream& out, context_line_filter omit_line) const {
		VINA_FOR_IN(i, ligands)
			write_context(ligands[i].cont, out, omit_line);
		if (num_flex() > 0)
			write_context(flex_context, out, omit_line);
	}
	void write_context(const context& c, std::ostream& out, context_line_filter omit_line = nullptr) const;
	void write_context(const context& c, ofile& out, const std::string& remark) const {
		out << remark;
	}
	void write_context(const context& c, const path& name) const {
		ofile out(name);
		write_context(c, out);
	}
	void write_context(const context& c, const path& name, const std::string& remark) const {
		ofile out(name);
		write_context(c, out, remark);
	}
	fl rmsd_lower_bound_asymmetric(const model& x, const model& y) const; // actually static
	
	atom_index sz_to_atom_index(sz i) const; // grid_atoms, atoms
	bool bonded_to_HD(const atom& a) const;
	bool bonded_to_heteroatom(const atom& a) const;
	szv bonded_to(sz a, sz n) const;
    bool is_unmatched_closure_dummy(sz i, sz j) const;
	bool is_glue_pair(sz i, sz j) const;

	atom& initializing_atom(const atom_index& i, atomv& rigid_atoms) { return i.in_grid ? rigid_atoms[i.i] : atoms[i.i]; }
	void assign_bonds(const distance_type_matrix& mobility, atomv& rigid_atoms);
	void assign_types(atomv& rigid_atoms);
	void initialize_pairs(const distance_type_matrix& mobility);
	vec ligand_heavy_centroid(sz ligand_number, const vecv& positions, const vec& origin) const;
	void initialize(const distance_type_matrix& mobility, atomv rigid_atoms);
	fl clash_penalty_aux(const interacting_pairs& pairs) const;

	vecv coords;
	vecv minus_forces;

	std::shared_ptr<const atomv> m_grid_atoms;
	atomv atoms; // movable, inflex
	vector_mutable<ligand> ligands;
	vector_mutable<residue> flex;
	context flex_context;
	interacting_pairs other_pairs; // INTRAmolecular interactions: flex_i - flex_j and flex_i - flex_i
	interacting_pairs inter_pairs; // INTERmolecular interactions: ligand - flex and ligand_i - ligand_j
	interacting_pairs glue_pairs; // INTRAmolecular interactions: glue_i - glue_i

	sz m_num_movable_atoms;
	atom_type::t m_atom_typing_used;
};

#endif
