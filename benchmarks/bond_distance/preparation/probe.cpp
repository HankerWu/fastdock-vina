// Time the same in-memory PDBQT parser used by the CLI and expose the full
// ordered pair tables (including lookup-table indices) as a separate operation.
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include "parse_pdbqt.h"

int main(int argc, char** argv) {
    if (argc != 5) return 2;
    const std::string mode(argv[1]), typing(argv[2]);
    const int repeats = std::stoi(argv[3]);
    if ((mode != "table" && mode != "time") || (typing != "XS" && typing != "AD") || repeats < 1)
        return 2;
    std::ifstream input(argv[4]);
    if (!input) throw std::runtime_error("Cannot open ligand");
    const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    const auto type = typing == "XS" ? atom_type::XS : atom_type::AD;
    const auto parse = [&]() { return parse_ligand_pdbqt_from_string(text, type); };
    if (mode == "table") {
        const model m = parse();
        const auto write = [](const char* name, const interacting_pairs& pairs) {
            for (const auto& p : pairs)
                std::cout << name << '\t' << p.a << '\t' << p.b << '\t' << p.table_index << '\n';
        };
        write("internal", m.get_ligand(0).pairs);
        write("glue", m.get_glue_pairs());
        write("other", m.get_other_pairs());
        write("inter", m.get_inter_pairs());
        return 0;
    }
    for (int i = 0; i < 8; ++i) parse();
    size_t total = 0;
    const auto before = std::chrono::steady_clock::now();
    const auto cpu = std::clock();
    for (int i = 0; i < repeats; ++i) total += parse().num_internal_pairs();
    const double cpu_per_parse = double(std::clock() - cpu) / CLOCKS_PER_SEC / repeats;
    const double wall_per_parse = std::chrono::duration<double>(std::chrono::steady_clock::now() - before).count() / repeats;
    std::cout << std::setprecision(12) << total << ' ' << cpu_per_parse << ' ' << wall_per_parse << '\n';
}
