// Measure the final rank-preserving geometry pass on real emitted pose pools.
#include "coords.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

static output_container read_pool(const char* filename) {
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("missing input");
    output_container pool;
    std::string line;
    while (std::getline(input, line)) {
        if (line.compare(0, 5, "MODEL") == 0) pool.push_back(new output_type(conf{}, 0));
        if (line.compare(0, 19, "REMARK VINA RESULT:") == 0)
            pool.back().e = std::stod(line.substr(19));
        if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0) continue;
        std::istringstream fields(line.substr(77));
        std::string type;
        fields >> type;
        if (type == "HD" || type == "H") continue;
        pool.back().coords.emplace_back(std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                                        std::stod(line.substr(46, 8)));
    }
    if (pool.empty() || pool.front().coords.empty()) throw std::runtime_error("empty pool");
    return pool;
}

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::vector<output_container> pools;
    for (int arg = 1; arg < argc; ++arg) pools.push_back(read_pool(argv[arg]));
    std::cout << "repeat,radius,pools,iterations,selected,us_per_pool\n";
    for (int round=0; round<7; ++round) {
        for (int phase=0; phase<2; ++phase) {
            const double radius = (phase+round)%2;
            std::size_t selected=0;
            const int repeat=200;
            const auto start=std::chrono::steady_clock::now();
            for (int i=0; i<repeat; ++i) {
                for (const auto& pool : pools)
                    selected += ranked_output_indices(pool,pool.size(),radius,boost::none).size();
            }
            const auto end=std::chrono::steady_clock::now();
            const double us=std::chrono::duration<double,std::micro>(end-start).count()/(pools.size()*repeat);
            std::cout << round+1 << ',' << radius << ',' << pools.size() << ',' << repeat << ','
                      << selected << ',' << std::setprecision(9) << us << '\n';
        }
    }
}
