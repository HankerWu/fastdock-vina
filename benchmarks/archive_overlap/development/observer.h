// Offline single-worker observation. The production binary never includes this file.
#include "coords.h"
#include "parse_pdbqt.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>

inline std::ostream& archive_trace() {
    static std::ofstream stream(std::getenv("FD_ARCHIVE_TRACE"));
    if (!stream) throw std::runtime_error("FD_ARCHIVE_TRACE must name a writable file");
    return stream << std::setprecision(17);
}

inline void trace_atom_order(const model& m) {
    const auto indices = ligand_heavy_indices_by_serial(m);
    std::map<sz, unsigned> ordered;
    for (const auto& item : indices) ordered.emplace(item.second, item.first);
    auto& stream = archive_trace();
    stream << "atoms";
    for (const auto& atom : ordered) stream << ' ' << atom.second;
    stream << '\n';
}

inline void trace_pool(std::ostream& stream, const output_container& out) {
    for (const auto& pose : out) {
        const auto& atom = pose.coords.front();
        stream << ' ' << pose.e << ',' << atom[0] << ',' << atom[1] << ',' << atom[2];
    }
}

inline void observed_insert(const char* stage, output_container& out, const output_type& pose,
                            fl radius, sz capacity) {
    auto& stream = archive_trace();
    stream << stage << ' ' << radius << ' ' << capacity << ' ' << pose.e;
    for (const auto& atom : pose.coords)
        stream << ' ' << atom[0] << ' ' << atom[1] << ' ' << atom[2];
    stream << '\t';
    trace_pool(stream, out);
    add_to_output_container(out, pose, radius, capacity);
    stream << '\t';
    trace_pool(stream, out);
    stream << '\n';
}
