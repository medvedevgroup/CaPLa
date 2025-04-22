#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <iomanip>
#include "piecewise_linear_model.hpp"
#include "cmdline.hpp"
#include "suffix_array.h"
#include "utils/essentials.hpp"
#include "utils/hasher/MurmurHash2.h"

void check_if_sorted(std::vector<uint64_t> &vec, std::string tag){
    std::cout<<tag<<std::endl;
    for(size_t i=1; i<vec.size(); i++){
        if(vec[i] < vec[i]-1){
            std::cerr<<"Not Sorted!"<<std::endl;
            exit(-1);
        }
    }
    std::cout<<"Sorted!"<<std::endl;
}

void get_hashed_kmer_vec(std::vector<uint64_t>& kmer_vec, std::vector<uint64_t>& hashed_kmer_vec){
    const uint64_t seed = 12345;
    hashed_kmer_vec.resize(kmer_vec.size());
    for(uint64_t i=0; i<kmer_vec.size(); i++){
        hashed_kmer_vec[i] = pla::default_hash64(kmer_vec[i], seed);
    }
    std::sort(hashed_kmer_vec.begin(), hashed_kmer_vec.end());
    check_if_sorted(hashed_kmer_vec, "inside hashed func");
}

void get_random_kmer_vec(uint64_t n, int kmer_size, std::vector<uint64_t>& random_kmer_vec){
    uint64_t _min = 0;
    uint64_t _max = 1;
    for(size_t i=0; i<2*kmer_size; i++) _max *= 2;
    std::random_device rd; // Seed for the random number engine
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_int_distribution<uint64_t> dis(_min, _max);
    
    random_kmer_vec.resize(n);
    for(size_t i=0; i<n; i++){
        random_kmer_vec[i] = dis(gen);
    }
    std::sort(random_kmer_vec.begin(), random_kmer_vec.end());
}

void get_wo_repeats_vec(std::vector<uint64_t>& kmer_vec, std::vector<uint64_t>& wo_rep_kmer_vec){

    wo_rep_kmer_vec.emplace_back(kmer_vec[0]);
    for(uint64_t i=1; i<kmer_vec.size(); i++){
        if(kmer_vec[i] != kmer_vec[i-1]){
            wo_rep_kmer_vec.emplace_back(kmer_vec[i]);
        }
    }
}

void get_kmer_vec(suffix_array<int64_t>& sa, std::vector<uint64_t>& kmer_vec){
    for(size_t i=0; i<sa.size(); i++){
        if(sa[i] != -1){
            kmer_vec.emplace_back(sa[i]);
        }
    }
}

void write_vec_to_file(std::string fn, std::vector<uint64_t>& vec){
    std::ofstream out(fn.c_str(), std::ios::binary);
    essentials::save_vec(out, vec);
    out.close();
}

int main(int argc, char **argv) {
    auto opt = parse_cla_kmer_type(argc, argv);

    std::string gn_fn, sa_fn, random_gn_fn, random_sa_fn, fn_prefix;
    int64_t kmer_size;

    gn_fn = opt.gn_fn;
    sa_fn = opt.sa_fn;
    random_gn_fn = opt.rn_gn_fn;
    random_sa_fn = opt.rn_sa_fn;
    kmer_size = opt.kmer_size;
    fn_prefix = opt.fn_prefix;

    std::string file_path = fn_prefix;

    suffix_array<int64_t> sa(gn_fn, sa_fn, kmer_size);
    suffix_array<int64_t> random_sa(random_gn_fn, random_sa_fn, kmer_size);
    // wr: w/o repeat, h: hashed, rg: random genome
    std::vector<uint64_t> kmer_vec, wr_kmer_vec, h_kmer_vec, wrh_kmer_vec, r_kmer_vec, rg_kmer_vec;

    auto start = std::chrono::high_resolution_clock::now();
    
    get_kmer_vec(sa, kmer_vec);
    get_hashed_kmer_vec(kmer_vec, h_kmer_vec);
    check_if_sorted(h_kmer_vec, "after function");
    get_wo_repeats_vec(kmer_vec, wr_kmer_vec);
    get_wo_repeats_vec(h_kmer_vec, wrh_kmer_vec);

    get_kmer_vec(random_sa, rg_kmer_vec);
    get_random_kmer_vec(sa.size(), kmer_size, r_kmer_vec);

    write_vec_to_file(file_path+".wr_kmer.bin", wr_kmer_vec);
    write_vec_to_file(file_path+".h_kmer.bin", h_kmer_vec);
    write_vec_to_file(file_path+".wrh_kmer.bin", wrh_kmer_vec);
    write_vec_to_file(file_path+".r_kmer.bin", r_kmer_vec);
    write_vec_to_file(file_path+".rg_kmer.bin", rg_kmer_vec);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Completed in " << duration_s << " seconds." << std::endl;

    return 0;
}