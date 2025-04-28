#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <iomanip>
#include <type_traits>
#include "piecewise_linear_model.hpp"
#include "cmdline.hpp"
#include "suffix_array.h"
#include "utils/essentials.hpp"

template <typename T>
uint64_t get_segments_count(size_t eps, const T &sa, bool connected = false) {
    int64_t prev_x = -1;

    auto in_fun = [&](auto i) {
        auto x = sa[i];
        if (x == -1) // x is suffix smaller than the k-mer size
            return std::pair<int64_t, int64_t>(prev_x, i); // segmentation will skip this point

        // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
        // to x=key[i] such that x+1!=key[i+1], we map the values x+1,...,key[i+1]-1 to their correct rank i
        auto flag = i > 0 && i + 1u < sa.size() && x == prev_x && x + 1 < sa[i + 1];
        prev_x = x + flag;

        return std::pair<int64_t, int64_t>(x + flag, i);
    };
    auto end_cond_indx = [&](){return false;};
    auto out_fun = [](auto) {};
    return make_segmentation_mod(sa.size(), eps, in_fun, out_fun, end_cond_indx, connected);
}

template <typename T>
uint64_t get_unique_elements(const T &vec){
    int64_t prev_elem = -1;
    uint64_t _count = 0;
    for(size_t i=0; i<vec.size(); i++){
        if(vec[i] == prev_elem || vec[i] == -1) continue;
        _count++;
        prev_elem = vec[i];
    }
    return _count;
}

template <typename T>
std::vector<size_t> eps_vector_up_to_single_segment(const T &sa, bool is_dict_point_seq = false) {
    std::vector<size_t> result;

    // Binary search for the smallest epsilon value giving one segment
    size_t left = 65536;
    size_t right = sa.size() / 4;
    while (left < right) {
        auto mid = (left + right) / 2;
        auto val = get_segments_count(mid, sa, is_dict_point_seq, false);
        if (val == 1) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    auto max_eps = left;

    // Add powers-of-two and equidistant epsilon values up to max_eps
    auto log2_max_eps = (int) std::log2(max_eps);
    for (size_t log2_eps = 0; log2_eps <= log2_max_eps; ++log2_eps)
        result.push_back(1ULL << log2_eps);
    auto step = std::max<size_t>(1, max_eps / 1024);
    for (size_t eps = step; eps < max_eps; eps += step)
        result.push_back(eps);
    result.push_back(max_eps);
    std::sort(result.begin(), result.end(), std::greater<>());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

template <typename T>
void write_segments_count(string out_fn, const T &sa,  bool use_all_epsilons = false) {
    std::ofstream out_file(out_fn.c_str());
    std::vector<size_t> eps_vec;
    if (use_all_epsilons) {
        for (int64_t i = 1024; i >= 1; i -= 1) {
            eps_vec.push_back(i);
        }
    } else {
        for (int64_t i = 10; i >= 0; i -= 1) {
            eps_vec.push_back(1ULL << i);
        }
    }

    // std::cout << std::setw(10) << "epsilon"
    //           << std::setw(20) << "#segments"
    //          << std::setw(20) << "#segments_connected"
    //           << std::endl;

    #pragma omp parallel for ordered schedule(dynamic)
    for (size_t i = 0; i < eps_vec.size(); i++) {
        auto count_disconnected = get_segments_count(eps_vec[i], sa, false);
        // auto count_connected = get_segments_count(eps_vec[i], sa, true);
        
        #pragma omp ordered
        {
            // std::cout << std::setw(10) << eps_vec[i]
            //           << std::setw(20) << count_disconnected
            //          << std::setw(20) << count_connected
            //           << std::endl;
            if (i == 0)
                out_file << get_unique_elements(sa) << std::endl;
            out_file << eps_vec[i] << " " << count_disconnected << std::endl;
            // out_file << eps_vec[i] << " " << count_connected << std::endl;
        }
    }
    out_file.close();
}

inline void process_suffix_array(std::string gn_fn, std::string sa_fn, int64_t kmer_size,
        std::string out_fn, bool use_all_epsilons) {
    suffix_array<int64_t> sa(gn_fn, sa_fn, kmer_size);
    write_segments_count(out_fn, sa, use_all_epsilons);
}

int main(int argc, char **argv) {
    auto opt = parse_cla_segment_count(argc, argv);

    
    bool is_genome_data = opt.is_genome_data;
    bool use_all_epsilons = opt.use_all_epsilons;
    
    string gn_fn, sa_fn, in_fn, out_fn;
    int64_t kmer_size;
    
    gn_fn = opt.gn_fn;
    sa_fn = opt.sa_fn;
    kmer_size = opt.kmer_size;
    out_fn = opt.count_fn;

    
    

    auto start = std::chrono::high_resolution_clock::now();
    process_suffix_array(gn_fn, sa_fn, kmer_size, out_fn, use_all_epsilons);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Completed in " << duration_s << " seconds." << std::endl;

    return 0;
}
