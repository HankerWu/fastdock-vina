#include "parse_pdbqt.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

// Compare the model passed to optimization, in PDBQT atom serial coordinates.
using Edge = std::pair<int, int>;
struct Topology {
    std::map<int, std::tuple<sz, sz, sz, double>> types;
    std::set<Edge> bonds, rotors;
    std::map<Edge, double> lengths;
    std::vector<std::tuple<int, int, int, sz>> pairs;
    std::vector<std::string> tree;
    unsigned dof;
};

Topology parsed(const std::string& path) {
    const model m = parse_ligand_pdbqt_from_file(path, atom_type::XS);
    if (m.num_ligands() != 1) throw std::runtime_error("Expected a single ligand: " + path);
    const auto lig = m.get_ligand(0);
    Topology t;
    t.dof = lig.degrees_of_freedom;
    std::vector<int> serial(m.num_atoms(), -1);
    for (const auto& line : lig.cont) {
        if (line.second) {
            const auto id = std::stoi(line.first.substr(6, 5));
            const auto index = *line.second;
            serial.at(index) = id;
            const auto a = m.get_atom(index);
            if (!t.types.emplace(id, std::make_tuple(a.el, a.ad, a.xs, a.charge)).second)
                throw std::runtime_error("Duplicate atom serial: " + path);
        } else if (line.first.compare(0, 6, "BRANCH") == 0 ||
                   line.first.compare(0, 9, "ENDBRANCH") == 0) {
            t.tree.push_back(line.first);
        }
    }
    if (std::find(serial.begin(), serial.end(), -1) != serial.end())
        throw std::runtime_error("Missing atom serial: " + path);
    for (sz i = 0; i < m.num_atoms(); ++i) {
        for (const auto& b : m.get_atom(i).bonds) {
            if (b.connected_atom_index.in_grid) throw std::runtime_error("Unexpected grid bond");
            const int other = serial.at(b.connected_atom_index.i);
            const Edge e = std::minmax(serial[i], other);
            t.bonds.insert(e);
            t.lengths[e] = b.length;
            if (b.rotatable) t.rotors.insert(e);
        }
    }
    auto collect = [&](const interacting_pairs& pairs, int category) {
        for (const auto& p : pairs) {
            const Edge e = std::minmax(serial.at(p.a), serial.at(p.b));
            t.pairs.emplace_back(category, e.first, e.second, p.table_index);
        }
    };
    collect(lig.pairs, 0);
    collect(m.get_glue_pairs(), 1);
    collect(m.get_inter_pairs(), 2);
    collect(m.get_other_pairs(), 3);
    std::sort(t.pairs.begin(), t.pairs.end());
    return t;
}

int main() {
    std::map<std::string, Topology> originals;
    std::string line;
    std::cout.precision(15);
    while (std::getline(std::cin, line)) {
        const auto split = line.find('\t');
        if (split == std::string::npos) throw std::runtime_error("Expected tab-separated paths");
        const auto original = line.substr(0, split), path = line.substr(split + 1);
        auto it = originals.find(original);
        if (it == originals.end()) it = originals.emplace(original, parsed(original)).first;
        const auto& reference = it->second;
        const auto pose = parsed(path);
        double delta = 0.0;
        for (const auto& b : pose.lengths) {
            const auto ref = reference.lengths.find(b.first);
            if (ref != reference.lengths.end())
                delta = std::max(delta, std::abs(b.second - ref->second));
        }
        std::cout << (pose.types != reference.types) << ' '
                  << (pose.bonds != reference.bonds) << ' '
                  << (pose.rotors != reference.rotors) << ' '
                  << (pose.dof != reference.dof || pose.tree != reference.tree) << ' '
                  << (pose.pairs != reference.pairs) << ' '
                  << delta << '\n';
    }
}
