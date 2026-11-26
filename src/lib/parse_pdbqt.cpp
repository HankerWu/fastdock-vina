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

#include <algorithm>
#include <string>
#include <fstream> // for getline ?
#include <sstream> // in parse_two_unsigneds
#include <cctype> // isspace
#include <exception>
#include <map>
#include <set>
#include <boost/utility.hpp> // for noncopyable 
#include <boost/optional.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
//#include <openbabel/mol.h>
//#include <openbabel/obconversion.h>
#include "model.h"
#include "atom_constants.h"
#include "file.h"
#include "convert_substring.h"
#include "utils.h"
#include "parse_pdbqt.h"
#include "parse_error.h"


struct parsed_atom : public atom {
    unsigned number; 
    parsed_atom(sz ad_, fl charge_, const vec& coords_, unsigned number_) : number(number_) {
        ad = ad_;
        charge = charge_;
        coords = coords_;
    }
};

void add_context(context& c, std::string& str) {
    c.push_back(parsed_line(str, boost::optional<sz>()));
}

std::string omit_whitespace(const std::string& str, sz i, sz j) {
    if(i < 1) i = 1;
    if(j < i-1) j = i-1; // i >= 1
    if(j < str.size()) j = str.size();

    // omit leading whitespace
    while(i <= j && std::isspace(str[i-1]))
        ++i;

    // omit trailing whitespace
    while(i <= j && std::isspace(str[j-1]))
        --j;

    VINA_CHECK(i-1 < str.size());
    VINA_CHECK(j-i+1 < str.size());

    return str.substr(i-1, j-i+1);
}

template<typename T>
T checked_convert_substring(const std::string& str, sz i, sz j, const std::string& dest_nature) {
    VINA_CHECK(i >= 1);
    VINA_CHECK(i <= j+1);

    if(j > str.size()) 
        throw pdbqt_parse_error("This line is too short.", str);

    // omit leading whitespace
    while(i <= j && std::isspace(str[i-1]))
        ++i;

    const std::string substr = str.substr(i-1, j-i+1);
    try {
        return boost::lexical_cast<T>(substr);
    }
    catch(...) {
        throw pdbqt_parse_error(dest_nature + std::string(" \"") + substr + "\" is not valid.", str);
    }
}

parsed_atom parse_pdbqt_atom_string(const std::string& str) {
    unsigned number = checked_convert_substring<unsigned>(str, 7, 11, "Atom number");
    vec coords(checked_convert_substring<fl>(str, 31, 38, "Coordinate"),
               checked_convert_substring<fl>(str, 39, 46, "Coordinate"),
               checked_convert_substring<fl>(str, 47, 54, "Coordinate"));
    fl charge = 0;
    if(!substring_is_blank(str, 69, 76))
        charge = checked_convert_substring<fl>(str, 69, 76, "Charge");
    std::string name = omit_whitespace(str, 78, 79);
    sz ad = string_to_ad_type(name);
    parsed_atom tmp(ad, charge, coords, number);

    if(is_non_ad_metal_name(name))
        tmp.xs = XS_TYPE_Met_D;
    if(tmp.acceptable_type()) 
        return tmp;
    else
        throw pdbqt_parse_error("Atom type " + name + " is not a valid AutoDock type (atom types are case-sensitive).", str);
}

struct atom_reference {
    sz index;
    bool inflex;
    atom_reference(sz index_, bool inflex_) : index(index_), inflex(inflex_) {}
};

struct movable_atom : public atom {
    vec relative_coords;
    movable_atom(const atom& a, const vec& relative_coords_) : atom(a) {
        relative_coords = relative_coords_;
    }
};

struct rigid {
    atomv atoms;
};

typedef std::vector<movable_atom> mav;

struct non_rigid_parsed {
    vector_mutable<ligand> ligands;
    vector_mutable<residue> flex;

    mav atoms;
    atomv inflex;

    distance_type_matrix atoms_atoms_bonds;
    matrix<distance_type> atoms_inflex_bonds;
    distance_type_matrix inflex_inflex_bonds;

    distance_type_matrix mobility_matrix() const {
        distance_type_matrix tmp(atoms_atoms_bonds);
        tmp.append(atoms_inflex_bonds, inflex_inflex_bonds);
        return tmp;
    }
};

struct parsing_struct {
    // start reading after this class
    template<typename T> // T == parsing_struct
    struct node_t {
        sz context_index;
        parsed_atom a;
        std::vector<T> ps;
        node_t(const parsed_atom& a_, sz context_index_) : context_index(context_index_), a(a_) {}

        // inflex atom insertion
        void insert_inflex(non_rigid_parsed& nr) {
            VINA_FOR_IN(i, ps)
                ps[i].axis_begin = atom_reference(nr.inflex.size(), true);
            nr.inflex.push_back(a);
        }
        void insert_immobiles_inflex(non_rigid_parsed& nr) {
            VINA_FOR_IN(i, ps)
                ps[i].insert_immobile_inflex(nr);
        }

        // insertion into non_rigid_parsed
        void insert(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
            VINA_FOR_IN(i, ps)
                ps[i].axis_begin = atom_reference(nr.atoms.size(), false);
            vec relative_coords; relative_coords = a.coords - frame_origin;
            c[context_index].second = nr.atoms.size();
            nr.atoms.push_back(movable_atom(a, relative_coords));
        }
        void insert_immobiles(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
            VINA_FOR_IN(i, ps)
                ps[i].insert_immobile(nr, c, frame_origin);
        }
    };

    typedef node_t<parsing_struct> node;
    boost::optional<sz> immobile_atom; // which of `atoms' is immobile, if any
    boost::optional<atom_reference> axis_begin; // the index (in non_rigid_parsed::atoms) of the parent bound to immobile atom (if already known)
    boost::optional<atom_reference> axis_end; // if immobile atom has been pushed into non_rigid_parsed::atoms, this is its index there
    std::vector<node> atoms;

    void add(const parsed_atom& a, const context& c) { 
        VINA_CHECK(c.size() > 0);
        atoms.push_back(node(a, c.size()-1)); 
    }
    const vec& immobile_atom_coords() const {
        VINA_CHECK(immobile_atom);
        VINA_CHECK(immobile_atom.get() < atoms.size());
        return atoms[immobile_atom.get()].a.coords;
    }
    // inflex insertion
    void insert_immobile_inflex(non_rigid_parsed& nr) {
        if(!atoms.empty()) {
            VINA_CHECK(immobile_atom);
            VINA_CHECK(immobile_atom.get() < atoms.size());
            axis_end = atom_reference(nr.inflex.size(), true);
            atoms[immobile_atom.get()].insert_inflex(nr);
        }
    }

    // insertion into non_rigid_parsed
    void insert_immobile(non_rigid_parsed& nr, context& c, const vec& frame_origin) {
        if(!atoms.empty()) {
            VINA_CHECK(immobile_atom);
            VINA_CHECK(immobile_atom.get() < atoms.size());
            axis_end = atom_reference(nr.atoms.size(), false);
            atoms[immobile_atom.get()].insert(nr, c, frame_origin);
        }
    }

    bool essentially_empty() const { // no sub-branches besides immobile atom, including sub-sub-branches, etc
        VINA_FOR_IN(i, atoms) {
            if(immobile_atom && immobile_atom.get() != i)
                return false;
            const node& nd = atoms[i];
            if(!nd.ps.empty())
                return false; // FIXME : iffy
        }
        return true;
    }
};

unsigned parse_one_unsigned(const std::string& str, const std::string& start) {
    std::istringstream in_str(str.substr(start.size()));
    int tmp;
    in_str >> tmp;

    if(!in_str || tmp < 0) 
        throw pdbqt_parse_error("Syntax error.", str);

    return unsigned(tmp);
}

void parse_two_unsigneds(const std::string& str, const std::string& start, unsigned& first, unsigned& second) {
    std::istringstream in_str(str.substr(start.size()));
    int tmp1, tmp2;
    in_str >> tmp1;
    in_str >> tmp2;

    if(!in_str || tmp1 < 0 || tmp2 < 0) 
        throw pdbqt_parse_error("Syntax error.", str);

    first = unsigned(tmp1);
    second = unsigned(tmp2);
}

void parse_pdbqt_rigid(const path& name, rigid& r) {
    ifile in(name);
    std::string str;

    while(std::getline(in, str)) {
        if(str.empty()) {} // ignore ""
        else if(starts_with(str, "TER")) {} // ignore 
        else if(starts_with(str, "END")) {} // ignore
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM"))
            r.atoms.push_back(parse_pdbqt_atom_string(str));
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in rigid receptor. "
                                    "Only one model can be used for the rigid receptor.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in rigid receptor.", str);
    }
}

void parse_pdbqt_root_aux(std::istream& in, parsing_struct& p, context& c) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // ignore ""
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM"))
            p.add(parse_pdbqt_atom_string(str), c);
        else if(starts_with(str, "ENDROOT"))
            return;
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

void parse_pdbqt_root(std::istream& in, parsing_struct& p, context& c) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // ignore
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "ROOT")) {
            parse_pdbqt_root_aux(in, p, c);
            break;
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

void parse_pdbqt_branch(std::istream& in, parsing_struct& p, context& c, unsigned from, unsigned to); // forward declaration

void parse_pdbqt_branch_aux(std::istream& in, const std::string& str, parsing_struct& p, context& c) {
    unsigned first, second;
    parse_two_unsigneds(str, "BRANCH", first, second);
    sz i = 0;

    for(; i < p.atoms.size(); ++i)
        if(p.atoms[i].a.number == first) {
            p.atoms[i].ps.push_back(parsing_struct());
            parse_pdbqt_branch(in, p.atoms[i].ps.back(), c, first, second);
            break;
        }

    if(i == p.atoms.size())
        throw pdbqt_parse_error("Atom number " + std::to_string(first) + " is missing in this branch.", str);
}

void parse_pdbqt_aux(std::istream& in, parsing_struct& p, context& c, boost::optional<unsigned>& torsdof, bool residue) {
    parse_pdbqt_root(in, p, c);

    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // ignore ""
        if(str[0] == '\0') {} // ignore a different kind of emptiness (potential issues on Windows)
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, str, p, c);
        else if(!residue && starts_with(str, "TORSDOF")) {
            if(torsdof)
                throw pdbqt_parse_error("TORSDOF keyword can be defined only once.");
            torsdof = parse_one_unsigned(str, "TORSDOF");
        }
        else if(residue && starts_with(str, "END_RES"))
            return; 
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}

void add_bonds(non_rigid_parsed& nr, boost::optional<atom_reference> atm, const atom_range& r) {
    if(atm)
        VINA_RANGE(i, r.begin, r.end) {
            atom_reference& ar = atm.get();
            if(ar.inflex) 
                nr.atoms_inflex_bonds(i, ar.index) = DISTANCE_FIXED; //(max_unsigned); // first index - atoms, second index - inflex
            else
                nr.atoms_atoms_bonds(ar.index, i) = DISTANCE_FIXED; // (max_unsigned);
        }
}

void set_rotor(non_rigid_parsed& nr, boost::optional<atom_reference> axis_begin, boost::optional<atom_reference> axis_end) {
    if(axis_begin && axis_end) {
        atom_reference& r1 = axis_begin.get();
        atom_reference& r2 = axis_end  .get();
        if(r2.inflex) {
            VINA_CHECK(r1.inflex); // no atom-inflex rotors
            nr.inflex_inflex_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
        }
        else
            if(r1.inflex)
                nr.atoms_inflex_bonds(r2.index, r1.index) = DISTANCE_ROTOR; // (atoms, inflex)
            else
                nr.atoms_atoms_bonds(r1.index, r2.index) = DISTANCE_ROTOR;
    }
}

typedef std::pair<sz, sz> axis_numbers;
typedef boost::optional<axis_numbers> axis_numbers_option;

void nr_update_matrixes(non_rigid_parsed& nr) {
    // atoms with indexes p.axis_begin and p.axis_end can not move relative to [b.node.begin, b.node.end)

    nr.atoms_atoms_bonds.resize(nr.atoms.size(), DISTANCE_VARIABLE);  
    nr.atoms_inflex_bonds.resize(nr.atoms.size(), nr.inflex.size(), DISTANCE_VARIABLE); // first index - inflex, second index - atoms
    nr.inflex_inflex_bonds.resize(nr.inflex.size(), DISTANCE_FIXED); // FIXME?
}

template<typename B> // B == branch or main_branch or flexible_body 
void postprocess_branch(non_rigid_parsed& nr, parsing_struct& p, context& c, B& b) {
    b.node.begin = nr.atoms.size();
    VINA_FOR_IN(i, p.atoms) {  // postprocess atoms into 'b.node'
        parsing_struct::node& p_node = p.atoms[i];
        if(p.immobile_atom && i == p.immobile_atom.get()) {} // skip immobile_atom - it's already inserted in "THERE"
        else p_node.insert(nr, c, b.node.get_origin());
        p_node.insert_immobiles(nr, c, b.node.get_origin());
    }
    b.node.end = nr.atoms.size();

    nr_update_matrixes(nr);
    add_bonds(nr, p.axis_begin, b.node); // b.node is used as atom_range
    add_bonds(nr, p.axis_end  , b.node); // b.node is used as atom_range
    set_rotor(nr, p.axis_begin, p.axis_end);

    VINA_RANGE(i, b.node.begin, b.node.end)
        VINA_RANGE(j, i+1, b.node.end)
            nr.atoms_atoms_bonds(i, j) = DISTANCE_FIXED; // FIXME


    VINA_FOR_IN(i, p.atoms) {   // postprocess children
        parsing_struct::node& p_node = p.atoms[i];
        VINA_FOR_IN(j, p_node.ps) {
            parsing_struct& ps = p_node.ps[j];
            if(!ps.essentially_empty()) { // immobile already inserted // FIXME ?!
                b.children.push_back(segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords, b.node)); // postprocess_branch will assign begin and end
                postprocess_branch(nr, ps, c, b.children.back());
            }
        }
    }
    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

void postprocess_ligand(non_rigid_parsed& nr, parsing_struct& p, context& c, unsigned torsdof) {
    VINA_CHECK(!p.atoms.empty());
    nr.ligands.push_back(ligand(flexible_body(rigid_body(p.atoms[0].a.coords, 0, 0)), torsdof)); // postprocess_branch will assign begin and end
    postprocess_branch(nr, p, c, nr.ligands.back());
    nr_update_matrixes(nr); // FIXME ?
}

void postprocess_residue(non_rigid_parsed& nr, parsing_struct& p, context& c) {
    VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
        parsing_struct::node& p_node = p.atoms[i];
        p_node.insert_inflex(nr);
        p_node.insert_immobiles_inflex(nr);
    }
    VINA_FOR_IN(i, p.atoms) { // iterate over "root" of a "residue"
        parsing_struct::node& p_node = p.atoms[i];
        VINA_FOR_IN(j, p_node.ps) {
            parsing_struct& ps = p_node.ps[j];
            if(!ps.essentially_empty()) { // immobile atom already inserted // FIXME ?!
                nr.flex.push_back(main_branch(first_segment(ps.immobile_atom_coords(), 0, 0, p_node.a.coords))); // postprocess_branch will assign begin and end
                postprocess_branch(nr, ps, c, nr.flex.back());
            }
        }
    }
    nr_update_matrixes(nr); // FIXME ?
    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_1() == nr.atoms.size());
    VINA_CHECK(nr.atoms_inflex_bonds.dim_2() == nr.inflex.size());
}

//dkoes, stream version
void parse_pdbqt_ligand(std::istream& in, non_rigid_parsed& nr, context& c) {
    parsing_struct p;
    boost::optional<unsigned> torsdof;

    parse_pdbqt_aux(in, p, c, torsdof, false);

    if(p.atoms.empty()) 
        throw pdbqt_parse_error("No atoms in this ligand.");
    if(!torsdof)
        throw pdbqt_parse_error("Missing TORSDOF keyword.");

    postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_ligand(const path& name, non_rigid_parsed& nr, context& c) {
    ifile in(name);
    parsing_struct p;
    boost::optional<unsigned> torsdof;

    parse_pdbqt_aux(in, p, c, torsdof, false);

    if(p.atoms.empty()) 
        throw pdbqt_parse_error("No atoms in this ligand.");
    if(!torsdof)
        throw pdbqt_parse_error("Missing TORSDOF keyword in this ligand.");

    postprocess_ligand(nr, p, c, unsigned(torsdof.get())); // bizarre size_t -> unsigned compiler complaint

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_residue(std::istream& in, parsing_struct& p, context& c) { 
    boost::optional<unsigned> dummy;
    parse_pdbqt_aux(in, p, c, dummy, true);
}

void parse_pdbqt_flex(const path& name, non_rigid_parsed& nr, context& c) {
    ifile in(name);
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} // ignore ""
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "BEGIN_RES")) {
            parsing_struct p;
            parse_pdbqt_residue(in, p, c);
            postprocess_residue(nr, p, c);
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue PDBQT file. "
                                    "Use \"vina_split\" to split flex residues in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue.", str);
    }

    VINA_CHECK(nr.atoms_atoms_bonds.dim() == nr.atoms.size());
}

void parse_pdbqt_branch(std::istream& in, parsing_struct& p, context& c, unsigned from, unsigned to) {
    std::string str;

    while(std::getline(in, str)) {
        add_context(c, str);

        if(str.empty()) {} //ignore ""
        else if(starts_with(str, "WARNING")) {} // ignore - AutoDockTools bug workaround
        else if(starts_with(str, "REMARK")) {} // ignore
        else if(starts_with(str, "BRANCH")) parse_pdbqt_branch_aux(in, str, p, c);
        else if(starts_with(str, "ENDBRANCH")) {
            unsigned first, second;
            parse_two_unsigneds(str, "ENDBRANCH", first, second);
            if(first != from || second != to) 
                throw pdbqt_parse_error("Inconsistent branch numbers.");
            if(!p.immobile_atom) 
                throw pdbqt_parse_error("Atom " + boost::lexical_cast<std::string>(to) + " has not been found in this branch.");
            return;
        }
        else if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
            parsed_atom a = parse_pdbqt_atom_string(str);
            if(a.number == to)
                p.immobile_atom = p.atoms.size();
            p.add(a, c);
        }
        else if(starts_with(str, "MODEL"))
            throw pdbqt_parse_error("Unexpected multi-MODEL tag found in flex residue or ligand PDBQT file. "
                                    "Use \"vina_split\" to split flex residues or ligands in multiple PDBQT files.");
        else 
            throw pdbqt_parse_error("Unknown or inappropriate tag found in flex residue or ligand.", str);
    }
}


//////////// new stuff //////////////////
bool is_vina_score_remark(const std::string& line) {
	// A ligand split from a previous docking output retains these diagnostic rows.
	for (const char* prefix : {"REMARK VINA RESULT:", "REMARK VINA PRO RESULT:", "REMARK VINA NATIVE:",
	                          "REMARK VINA AMIDE RESULT:", "REMARK VINA AMIDE PENALTY:", "REMARK VINA CONTACT ROTOR SCORE:",
	                          "REMARK VINA PRO GAUSS2 SCALE:", "REMARK VINA TEMPLATE RESULT:",
	                          "REMARK VINA TEMPLATE STAGE:", "REMARK VINA TEMPLATE OVERLAP:",
	                          "REMARK VINA POLAR RESULT:", "REMARK VINA POLAR SCORES:",
	                          "REMARK VINA POLAR TERMS:", "REMARK VINA POLAR MODE:",
	                          "REMARK VINA POLAR BURIAL RESULT:", "REMARK VINA POLAR BURIAL TERMS:",
	                          "REMARK VINA POLAR BURIAL DELTAS:", "REMARK VINA POLAR BURIAL MODE:",
	                          "REMARK VINA POLAR BURIAL NATIVE RANK:",
	                          "REMARK VINA NET POLAR RESULT:", "REMARK VINA NET POLAR TERMS:",
	                          "REMARK VINA NET POLAR DELTAS:", "REMARK VINA NET POLAR MODE:",
	                          "REMARK VINA NET POLAR WATER POINTS:", "REMARK VINA NET POLAR SATISFACTION:",
	                          "REMARK VINA ARCHIVE RESULT:", "REMARK VINA ARCHIVE SUPPORT:",
	                          "REMARK VINA ARCHIVE SOURCE:", "REMARK VINA ARCHIVE MODE:",
	                          "REMARK INTER + INTRA:", "REMARK INTER:",
	                          "REMARK INTRA:", "REMARK CONF_INDEPENDENT:", "REMARK UNBOUND:"}) {
		if (starts_with(line, prefix))
			return true;
	}
	return false;
}

struct pdbqt_initializer {

    atom_type::t atom_typing_used;
    model m;
    //pdbqt_initializer(): atom_typing_used(atom_type::XS), m(atom_type::XS) {}
    pdbqt_initializer(atom_type::t atype): atom_typing_used(atype), m(atype) {}
    void initialize_from_nrp(const non_rigid_parsed& nrp, const context& c, bool is_ligand) { // static really
        VINA_CHECK(m.ligands.empty());
        VINA_CHECK(m.flex   .empty());

        m.ligands = nrp.ligands;
        m.flex    = nrp.flex;

        VINA_CHECK(m.atoms.empty());

        sz n = nrp.atoms.size() + nrp.inflex.size();
        m.atoms.reserve(n);
        m.coords.reserve(n);

        VINA_FOR_IN(i, nrp.atoms) {
            const movable_atom& a = nrp.atoms[i];
            atom b = static_cast<atom>(a);
            b.coords = a.relative_coords;
            m.atoms.push_back(b);
            m.coords.push_back(a.coords);
        }
        VINA_FOR_IN(i, nrp.inflex) {
            const atom& a = nrp.inflex[i];
            atom b = a;
            b.coords = zero_vec; // to avoid any confusion; presumably these will never be looked at
            m.atoms.push_back(b);
            m.coords.push_back(a.coords);
        }
        VINA_CHECK(m.coords.size() == n);

        m.minus_forces = m.coords;
        m.m_num_movable_atoms = nrp.atoms.size();

        if(is_ligand) {
            VINA_CHECK(m.ligands.size() == 1);
            auto& stored = m.ligands.front().cont;
            stored = c;
            // Split amide-ranked or contact-annotated poses carry scores for old coordinates.
            // Discard that score block on input, including when ordinary docking or
            // single-pose local refinement writes the molecule back out.
            if (std::any_of(stored.begin(), stored.end(), [](const parsed_line& l) {
                return starts_with(l.first, "REMARK VINA AMIDE RESULT:") ||
                       starts_with(l.first, "REMARK VINA CONTACT ROTOR SCORE:");
            })) {
                stored.erase(std::remove_if(stored.begin(), stored.end(), [](const parsed_line& l) {
                    return is_vina_score_remark(l.first);
                }), stored.end());
            }
        }
        else
            m.flex_context = c;

    }
    void initialize(const distance_type_matrix& mobility, atomv rigid_atoms = {}) {
        m.initialize(mobility, std::move(rigid_atoms));
    }
    static void replace_fixed_pose_coordinates(model& prepared, const non_rigid_parsed& pose,
                                               const context& pose_context);
    void initialize_explicit_ligand_bonds(const distance_type_matrix& mobility, std::istream& input);
};

namespace {

unsigned fixed_pose_atom_serial(const std::string& line) {
    auto token = line.substr(6, 5);
    const auto first = token.find_first_not_of(" ");
    const auto last = token.find_last_not_of(" ");
    if (first == std::string::npos)
        throw pdbqt_parse_error("PDBQT atom serial must be a positive decimal integer.", line);
    token = token.substr(first, last - first + 1);
    if (token.find_first_not_of("0123456789") != std::string::npos)
        throw pdbqt_parse_error("PDBQT atom serial must be a positive decimal integer.", line);
    const auto serial = checked_convert_substring<unsigned>(token, 1, token.size(), "Atom number");
    if (!serial)
        throw pdbqt_parse_error("PDBQT atom serial must be a positive decimal integer.", line);
    return serial;
}

// Normalize layout tokens but keep their order and the full chemical metadata.
std::string fixed_pose_layout_line(const parsed_line& line) {
    if (line.second) {
        // Keep the atom identity (including serial), charge and AD type.
        const auto& text = line.first;
        return text.substr(0, 30) + text.substr(68, 11);
    }
    std::istringstream input(line.first);
    std::string tag, layout, token;
    input >> tag;
    if (tag == "REMARK") {
        input >> token;
        if (token != "SMILES" && token != "H" && token != "INDEX")
            return {};
        layout = tag + " " + token;
    } else {
        if (tag != "ROOT" && tag != "ENDROOT" && tag != "BRANCH" && tag != "ENDBRANCH" && tag != "TORSDOF")
            return {};
        layout = tag;
    }
    while (input >> token)
        layout += " " + token;
    return layout;
}

using fixed_pose_layout = std::vector<std::pair<std::string, const parsed_line*>>;

fixed_pose_layout fixed_pose_layout_of(const context& lines) {
    fixed_pose_layout result;
    for (const auto& line : lines) {
        auto signature = fixed_pose_layout_line(line);
        if (!signature.empty())
            result.emplace_back(std::move(signature), &line);
    }
    return result;
}

std::string unwrap_fixed_pose(const std::string& input) {
    std::istringstream stream(input);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        lines.push_back(line);
    }
    const auto content = [](const std::string& s) { return s.find_first_not_of(" \t") != std::string::npos; };
    const auto first = std::find_if(lines.begin(), lines.end(), content);
    if (first == lines.end())
        throw pdbqt_parse_error("No atoms in the supplied pose coordinates.");
    if (!starts_with(*first, "MODEL"))
        return input;
    const auto last = std::find_if(lines.rbegin(), lines.rend(), content).base() - 1;
    std::istringstream header(*first), footer(*last);
    std::string tag, extra;
    unsigned model_number = 0;
    if (!(header >> tag >> model_number) || tag != "MODEL" || !model_number || (header >> extra) ||
        !(footer >> tag) || tag != "ENDMDL" || (footer >> extra) || first == last)
        throw pdbqt_parse_error("Pose coordinates require exactly one complete MODEL or bare ligand PDBQT.");
    std::ostringstream body;
    for (auto it = first + 1; it != last; ++it) {
        if (starts_with(*it, "MODEL") || starts_with(*it, "ENDMDL"))
            throw pdbqt_parse_error("Pose coordinates require exactly one complete MODEL or bare ligand PDBQT.");
        body << *it << '\n';
    }
    return body.str();
}

} // namespace

std::map<unsigned, sz> ligand_heavy_indices_by_serial(const model& m) {
    if (m.num_ligands() != 1)
        throw std::invalid_argument("search_anchor requires one prepared ligand");
    std::map<unsigned, sz> heavy;
    std::set<unsigned> used;
    const auto ligand = m.get_ligand(0);
    for (const auto& line : ligand.cont) {
        if (!line.second) continue;
        const auto serial = fixed_pose_atom_serial(line.first);
        if (!used.insert(serial).second)
            throw std::invalid_argument("search_anchor atom serial is ambiguous");
        const sz index = *line.second;
        if (index >= m.num_movable_atoms()) continue;
        const auto& atom = m.movable_atom(index);
        if (atom.el != EL_TYPE_H && atom.el != EL_TYPE_Dummy)
            heavy.emplace(serial, index);
    }
    return heavy;
}

void pdbqt_initializer::replace_fixed_pose_coordinates(model& prepared, const non_rigid_parsed& pose,
                                                      const context& pose_context) {
    const context& prepared_context = prepared.ligands.front().cont;
    const auto original = fixed_pose_layout_of(prepared_context);
    const auto candidate = fixed_pose_layout_of(pose_context);
    if (original.size() != candidate.size() || pose.atoms.size() != prepared.num_movable_atoms())
        throw pdbqt_parse_error("Pose coordinates differ from the prepared atom layout or chemical metadata.");

    std::set<unsigned> serials;
    vecv coordinates = prepared.coords;
    VINA_FOR_IN(i, original) {
        const auto& reference = *original[i].second;
        const auto& current = *candidate[i].second;
        if (original[i].first != candidate[i].first || bool(reference.second) != bool(current.second))
            throw pdbqt_parse_error("Pose coordinates differ from the prepared atom layout or chemical metadata.",
                                    current.first);
        if (!reference.second)
            continue;
        const auto serial = fixed_pose_atom_serial(current.first);
        if (!serials.insert(serial).second)
            throw pdbqt_parse_error("Pose coordinates require unique positive atom serials.", current.first);
        const vec& xyz = pose.atoms[*current.second].coords;
        if (!std::isfinite(xyz[0]) || !std::isfinite(xyz[1]) || !std::isfinite(xyz[2]))
            throw pdbqt_parse_error("Pose coordinates must be finite.", current.first);
        coordinates[*reference.second] = xyz;
    }
    if (serials.size() != pose.atoms.size())
        throw pdbqt_parse_error("Pose coordinates do not map every prepared atom.");
    prepared.coords.swap(coordinates);
}

void replace_fixed_pose_coordinates(model& prepared, const std::string& pose) {
    VINA_CHECK(prepared.num_ligands() == 1 && prepared.num_flex() == 0);
    std::istringstream input(unwrap_fixed_pose(pose));
    non_rigid_parsed parsed;
    context lines;
    // Parse/validate the pose tree without calling initialize(), which assigns
    // bonds, XS types and interaction pairs from the pose coordinates.
    parse_pdbqt_ligand(input, parsed, lines);
    pdbqt_initializer::replace_fixed_pose_coordinates(prepared, parsed, lines);
}


void pdbqt_initializer::initialize_explicit_ligand_bonds(const distance_type_matrix& mobility,
                                                          std::istream& input) {
    std::map<unsigned, sz> indices;
    const auto& ligand = m.ligands.front();
    for (const auto& line : ligand.cont) {
        if (!line.second)
            continue;
        const auto serial = fixed_pose_atom_serial(line.first);
        const sz index = *line.second;
        if (!indices.emplace(serial, index).second)
            throw pdbqt_parse_error("ligand_bonds requires unique positive PDBQT atom serials.", line.first);
        const sz type = m.atoms[index].ad;
        if (type >= AD_TYPE_G0 && type <= AD_TYPE_W)
            throw pdbqt_parse_error("ligand_bonds does not support macrocycle or hydrated-ligand atom types.", line.first);
        const vec& xyz = m.coords[index];
        if (!std::isfinite(xyz[0]) || !std::isfinite(xyz[1]) || !std::isfinite(xyz[2]))
            throw pdbqt_parse_error("ligand_bonds requires finite prepared coordinates.", line.first);
    }
    if (indices.size() != m.atoms.size())
        throw pdbqt_parse_error("ligand_bonds requires a serial for every prepared atom.");

    std::set<std::pair<sz, sz>> edges;
    std::string line;
    while (std::getline(input, line)) {
        std::istringstream fields(line.substr(0, line.find('#')));
        std::string first_token, second_token, extra;
        if (!(fields >> first_token))
            continue;
        if (!(fields >> second_token) || (fields >> extra))
            throw pdbqt_parse_error("ligand_bonds expects two atom serials on each edge line.", line);
        const auto serial = [&](const std::string& token) {
            if (token.find_first_not_of("0123456789") != std::string::npos)
                throw pdbqt_parse_error("ligand_bonds atom serials must be positive integers.", line);
            const auto value = checked_convert_substring<unsigned>(token, 1, token.size(), "Atom number");
            const auto found = indices.find(value);
            if (found == indices.end())
                throw pdbqt_parse_error("ligand_bonds edge references an unknown prepared atom serial.", line);
            return found->second;
        };
        const sz first = serial(first_token), second = serial(second_token);
        if (first == second)
            throw pdbqt_parse_error("ligand_bonds edge cannot connect an atom to itself.", line);
        const auto edge = std::make_pair(std::min(first, second), std::max(first, second));
        if (mobility(edge.first, edge.second) == DISTANCE_VARIABLE)
            throw pdbqt_parse_error("ligand_bonds edge crosses a variable-distance torsion-tree pair.", line);
        if (!edges.insert(edge).second)
            throw pdbqt_parse_error("Duplicate ligand_bonds edge.", line);
    }
    VINA_FOR_IN(first, m.atoms) {
        VINA_RANGE(second, first + 1, m.atoms.size()) {
            if (mobility(first, second) == DISTANCE_ROTOR && !edges.count({first, second}))
                throw pdbqt_parse_error("ligand_bonds is missing a BRANCH rotor edge.");
        }
    }
    for (const auto& edge : edges) {
        const sz a = edge.first, b = edge.second;
        const fl length = std::sqrt(vec_distance_sqr(m.coords[a], m.coords[b]));
        if (!(length > epsilon_fl) || !std::isfinite(length))
            throw pdbqt_parse_error("ligand_bonds bond endpoints have invalid coordinates.");
        const bool rotatable = mobility(a, b) == DISTANCE_ROTOR;
        m.atoms[a].bonds.emplace_back(atom_index(b, false), length, rotatable);
        m.atoms[b].bonds.emplace_back(atom_index(a, false), length, rotatable);
    }
    for (const auto& atom : m.atoms) {
        if ((atom.ad == AD_TYPE_H || atom.ad == AD_TYPE_HD) && atom.bonds.size() != 1)
            throw pdbqt_parse_error("ligand_bonds requires one neighbor for each explicit hydrogen.");
    }

    std::vector<bool> visited(m.atoms.size(), false);
    szv connected{0};
    visited[0] = true;
    for (sz i = 0; i < connected.size(); ++i) {
        for (const auto& bond : m.atoms[connected[i]].bonds) {
            const sz next = bond.connected_atom_index.i;
            if (!visited[next]) {
                visited[next] = true;
                connected.push_back(next);
            }
        }
    }
    if (connected.size() != m.atoms.size())
        throw pdbqt_parse_error("ligand_bonds must connect every atom of one prepared ligand.");

    m.ligands.front().set_range();
    atomv rigid_atoms;
    m.assign_types(rigid_atoms);
    m.initialize_pairs(mobility);
}

model parse_ligand_pdbqt_with_bonds(const std::string& ligand_name, atom_type::t atype,
                                  const std::string& bonds_name) {
    non_rigid_parsed parsed;
    context lines;
    parse_pdbqt_ligand(make_path(ligand_name), parsed, lines);
    pdbqt_initializer prepared(atype);
    prepared.initialize_from_nrp(parsed, lines, true);
    ifile bonds(make_path(bonds_name));
    prepared.initialize_explicit_ligand_bonds(parsed.mobility_matrix(), bonds);
    return prepared.m;
}

model parse_ligand_pdbqt_from_file(const std::string& name, atom_type::t atype) { // can throw parse_error
    non_rigid_parsed nrp;
    context c;
    parse_pdbqt_ligand(make_path(name), nrp, c);

    pdbqt_initializer tmp(atype);
    tmp.initialize_from_nrp(nrp, c, true);
    tmp.initialize(nrp.mobility_matrix());
    return std::move(tmp.m);
}

model parse_ligand_pdbqt_from_string(const std::string& string_name, atom_type::t atype) { // can throw parse_error
    non_rigid_parsed nrp;
    context c;

    std::stringstream molstream(string_name);
    parse_pdbqt_ligand(molstream, nrp, c);

    pdbqt_initializer tmp(atype);
    tmp.initialize_from_nrp(nrp, c, true);
    tmp.initialize(nrp.mobility_matrix());
    return std::move(tmp.m);
}

model parse_receptor_pdbqt(const std::string& rigid_name, const std::string& flex_name, atom_type::t atype) { 
    // Parse PDBQT receptor with flex residues
    //if (rigid_name.empty() && flex_name.empty()) {
    //    // CONDITION 1
    //    std::cerr << "ERROR: No (rigid) receptor or flexible residues were specified.\n";
    //    exit(EXIT_FAILURE);
    //}

    rigid r;
    non_rigid_parsed nrp;
    context c;
    pdbqt_initializer tmp(atype);

    if (!rigid_name.empty()) {
        parse_pdbqt_rigid(make_path(rigid_name), r);
    }
    
    if (!flex_name.empty()) {
        parse_pdbqt_flex(make_path(flex_name), nrp, c);
    }

    if (!flex_name.empty()) {
        tmp.initialize_from_nrp(nrp, c, false);
        tmp.initialize(nrp.mobility_matrix(), std::move(r.atoms));
    } else if (!rigid_name.empty())
        tmp.initialize(distance_type_matrix(), std::move(r.atoms));

    return std::move(tmp.m);
}
