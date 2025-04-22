#include <iostream>
#include <string>
#include <random>
#include "utils/essentials.hpp"

std::string generate_random_genome(size_t N, unsigned int seed) {
    const char nucleotides[] = {'A', 'T', 'C', 'G'};
    std::string genome;
    genome.reserve(N + 20);

    std::mt19937 gen(seed);
    std::uniform_int_distribution<> dis(0, 3);

    genome = ">Random genome: "+std::to_string(N)+"\n";

    for (size_t i = 0; i < N; ++i) {
        genome += nucleotides[dis(gen)];
    }

    return genome;
}


int main(int argc, char **argv) {
    size_t N = std::stoi(argv[1]);
    std::string out_fn = argv[2];

    unsigned int seed = 12345;
    std::string genome = generate_random_genome(N, seed);
    
    std::ofstream os(out_fn.c_str());
    os<<genome;
    os.close();
    
    return 0;
}