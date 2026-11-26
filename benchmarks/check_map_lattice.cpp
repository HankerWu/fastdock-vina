// Regression for imported multi-channel Vina/Vinardo and AD4 lattices.
// Link against the core library objects (without the command-line main).
#include <iostream>
#include <stdexcept>
#include "cache.h"
#include "ad4cache.h"

namespace fs = boost::filesystem;

static void write_map(const std::string& prefix, const std::string& type,
                      double center = 0, double spacing = 0.5, int nx = 2) {
    ofile out(path(prefix + "." + type + ".map"));
    out << "GRID_PARAMETER_FILE none\nGRID_DATA_FILE none\nMACROMOLECULE none\n"
        << "SPACING " << spacing << "\nNELEMENTS " << nx << " 2 2\nCENTER "
        << center << " 0 0\n";
    for (int i = 0; i < (nx + 1) * 3 * 3; ++i)
        out << i << '\n';
}

template<class Cache>
static void rejects(Cache& cache, const std::string& prefix) {
    try {
        cache.read(prefix);
    } catch (const std::runtime_error& error) {
        if (std::string(error.what()).find("Inconsistent affinity map lattice") != std::string::npos)
            return;
        throw;
    }
    throw std::runtime_error("Incompatible map set accepted: " + prefix);
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    // Each invocation uses a fresh isolated directory to avoid stale fixtures.
    const path root(argv[1]);
    if (!fs::create_directory(root)) {
        std::cerr << "Fixture directory must not exist\n";
        return 2;
    }
    const auto prefix = [&](const std::string& name) { return (root / name).string(); };
    try {
        cache xs;
        write_map(prefix("xs_good"), "C_H");
        write_map(prefix("xs_good"), "C_P");
        xs.read(prefix("xs_good"));
        if (!xs.is_atom_type_grid_initialized(XS_TYPE_C_H) ||
            !xs.is_atom_type_grid_initialized(XS_TYPE_C_P))
            throw std::runtime_error("valid XS maps missing");
        const auto original = xs.gd();

        write_map(prefix("xs_bad"), "C_H", 5);
        write_map(prefix("xs_bad"), "C_P", 5.25);
        rejects(xs, prefix("xs_bad"));
        if (!eq(xs.gd(), original) || !xs.is_atom_type_grid_initialized(XS_TYPE_C_P))
            throw std::runtime_error("failed XS reload changed cache");

        write_map(prefix("xs_new"), "C_H", 3);
        xs.read(prefix("xs_new"));
        if (xs.is_atom_type_grid_initialized(XS_TYPE_C_P) ||
            xs.gd()[0].begin != 2.5)
            throw std::runtime_error("successful XS reload retained old lattice");

        ad4cache ad;
        write_map(prefix("ad_good"), "C");
        write_map(prefix("ad_good"), "OA");
        write_map(prefix("ad_good"), "e");
        write_map(prefix("ad_good"), "d");
        ad.read(prefix("ad_good"));
        if (!ad.are_atom_types_grid_initialized({AD_TYPE_C, AD_TYPE_OA}))
            throw std::runtime_error("valid AD4 maps missing");

        // Three distinct incompatibilities: spatial origin, sample spacing,
        // and voxel count. Exercise every AD4 channel category.
        write_map(prefix("ad_origin"), "C", 0.25);
        write_map(prefix("ad_origin"), "e");
        write_map(prefix("ad_origin"), "d");
        rejects(ad, prefix("ad_origin"));

        write_map(prefix("ad_spacing"), "C");
        write_map(prefix("ad_spacing"), "e", 0, 0.25);
        write_map(prefix("ad_spacing"), "d");
        rejects(ad, prefix("ad_spacing"));

        write_map(prefix("ad_count"), "C");
        write_map(prefix("ad_count"), "e");
        write_map(prefix("ad_count"), "d", 0, 0.25, 4);
        rejects(ad, prefix("ad_count"));
        if (!eq(ad.gd(), original) || !ad.is_atom_type_grid_initialized(AD_TYPE_OA))
            throw std::runtime_error("failed AD4 reload changed cache");

        write_map(prefix("ad_new"), "C", 3);
        write_map(prefix("ad_new"), "e", 3);
        write_map(prefix("ad_new"), "d", 3);
        ad.read(prefix("ad_new"));
        if (ad.is_atom_type_grid_initialized(AD_TYPE_OA) ||
            !ad.are_atom_types_grid_initialized({AD_TYPE_C}) ||
            ad.gd()[0].begin != 2.5)
            throw std::runtime_error("successful AD4 reload retained old lattice");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        fs::remove_all(root);
        return 1;
    }
    fs::remove_all(root);
    std::cout << "Consistent map sets and transactional reload checks passed\n";
}
