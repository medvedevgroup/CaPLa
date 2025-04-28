/**
 * inspired from https://github.com/ksahlin/strobealign/blob/main/src/cmdline.hpp
*/

#pragma once

#include <vector>
#include <string>
#include <utility>

struct CommandLineOptions {
    std::string gn_fn;
    std::string sa_fn;

    std::string rn_gn_fn;
    std::string rn_sa_fn;
    std::string fn_prefix;

    int64_t kmer_size {21};
    std::string count_fn {"-1"};
    
    std::string in_fn {"-1"};

    double theta {0.1};
    bool is_dict_point_seq {false};
    bool is_sketch_data {false};
    bool is_genome_data {true};
    bool is_int_vec {false};
    bool use_all_epsilons {false};

    std::string out_fn {"-1"};

    std::string segment_file;
    double increment {1E-6};
    
};

CommandLineOptions parse_cla_segment_count(int argc, char **argv);
