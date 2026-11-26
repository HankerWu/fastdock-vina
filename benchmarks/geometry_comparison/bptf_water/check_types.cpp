// Check parser-derived XS types and connectivity for appended water sites.
#include "parse_pdbqt.h"
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto dry = parse_receptor_pdbqt(argv[1]).grid_atoms();
    const auto steric = parse_receptor_pdbqt(argv[2]).grid_atoms();
    const auto polar = parse_receptor_pdbqt(argv[3]).grid_atoms();
    const sz sites = steric.size() - dry.size();
    if (!sites || polar.size() != dry.size() + 3 * sites)
        throw std::runtime_error("Unexpected water count");
    for (sz i = 0; i < dry.size(); ++i)
        if (steric[i].xs != dry[i].xs || polar[i].xs != dry[i].xs ||
            steric[i].bonds.size() != dry[i].bonds.size() || polar[i].bonds.size() != dry[i].bonds.size())
            throw std::runtime_error("Water changed prepared protein types or connectivity");
    for (sz i = 0; i < sites; ++i) {
        const auto& o = steric[dry.size()+i];
        const sz index = dry.size()+3*i;
        const auto& h2o = polar[index];
        if (o.xs != XS_TYPE_O_P || !o.bonds.empty() || h2o.xs != XS_TYPE_O_DA || h2o.bonds.size() != 2)
            throw std::runtime_error("Wrong XS water oxygen type or extra covalent bond");
        for (sz h : {index+1,index+2})
            if (polar[h].ad != AD_TYPE_HD || polar[h].bonds.size() != 1 ||
                polar[h].bonds[0].connected_atom_index.i != index)
                throw std::runtime_error("Water hydrogen misassigned");
    }
    std::cout << sites << " independently recurrent sites verified\n";
}
