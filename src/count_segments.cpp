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

// blue definition for dictionary (n points)
template <typename T>
uint64_t get_segments_count(size_t eps, const T &sa, 
    bool is_dict_point_seq = false, bool connected = false) {
    int64_t prev_x = -1;

    // needed for 64 bit values
    // bool is_signed_value = false;
    // using TT = typename decltype(sa)::value_type;
    // if (std::is_signed<TT>::value) {
    //     is_signed_value = true;
    // } 

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
    
    int64_t indx = 0, rank = 0;
    // this will be called first
    auto end_cond = [&](){
        if (indx >= sa.size()) return true;
        auto x = sa[indx++];
        while(x == -1 || x == prev_x) {
            if(indx == sa.size()) return true;
            x = sa[indx++];
        }
        return false;
    };

    // only for dict for now
    // this will be called second
    auto next_fun = [&](auto _){
        auto x = sa[indx-1];
        int64_t diff = x - prev_x;
        prev_x = x;
        return std::pair<int64_t, int64_t>(rank++, x); 
    };

    auto out_fun = [](auto) {};

    return !is_dict_point_seq
        ? make_segmentation_mod(sa.size(), eps, in_fun, out_fun, end_cond_indx, connected)
        : make_segmentation_mod(sa.size(), eps, next_fun, out_fun, end_cond, connected);
}

template <typename T>
uint64_t get_segments_count_int(size_t eps, const T &sa, 
    bool is_dict_point_seq = false, bool connected = false) {
    int64_t prev_x = -1;

    auto in_fun = [&](auto i) {
        auto x = sa[i];
        // if (is_signed_value && x == -1) // x is suffix smaller than the k-mer size
        //     return std::pair<int64_t, int64_t>(prev_x, i); // segmentation will skip this point

        // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
        // to x=key[i] such that x+1!=key[i+1], we map the values x+1,...,key[i+1]-1 to their correct rank i
        auto flag = i > 0 && i + 1u < sa.size() && x == prev_x && x + 1 < sa[i + 1];
        prev_x = x + flag;

        return std::pair<uint64_t, uint64_t>(x + flag, i);
    };
    auto end_cond_indx = [&](){return false;};
    
    int64_t indx = 0, rank = 0;
    // this will be called first
    auto end_cond = [&](){
        if (indx >= sa.size()) return true;
        auto x = sa[indx++];
        while(x == -1 || x == prev_x) {
            if(indx == sa.size()) return true;
            x = sa[indx++];
        }
        return false;
    };

    // only for dict for now
    // this will be called second
    auto next_fun = [&](auto _){
        auto x = sa[indx-1];
        int64_t diff = x - prev_x;
        prev_x = x;
        return std::pair<int64_t, int64_t>(rank++, x); 
    };

    auto out_fun = [](auto) {};

    return !is_dict_point_seq
        ? make_segmentation_mod(sa.size(), eps, in_fun, out_fun, end_cond_indx, connected)
        : make_segmentation_mod(sa.size(), eps, next_fun, out_fun, end_cond, connected);
}

// violet definition for dictionary (N points)
template <typename T>
uint64_t get_segments_count_2(size_t eps, const T &sa, 
    bool is_dict_point_seq = false, bool connected = false) {
    auto out_fun = [](auto) {};
    auto end_cond = [&](){return false;};
    int64_t prev_x = -1;
    if (is_dict_point_seq) {
        int64_t prev_i = -1;
        auto in_fun_dict = [&](auto i) {
            auto x = sa[i];
            auto is_not_full_kmer = x == -1;
            auto is_repetition = i > 0 && x == prev_x;
            if (is_not_full_kmer || is_repetition)
                return std::pair<int64_t, int64_t>(prev_i, -1); // segmentation will skip this point
            prev_x = x;
            prev_i = i;
            return std::pair<int64_t, int64_t>(i, x);
        };
        return make_segmentation_mod(sa.size(), eps, in_fun_dict, out_fun, end_cond, connected);
    } else {
        auto in_fun_index = [&](auto i) {
            auto x = sa[i];
            if (x == -1) // x is suffix smaller than the k-mer size
                return std::pair<int64_t, int64_t>(prev_x, -1); // segmentation will skip this point

            // Here there is an adjustment for inputs with duplicate keys: at the end of a run of duplicate keys equal
            // to x=key[i] such that x+1!=key[i+1], we map the values x+1,...,key[i+1]-1 to their correct rank i
            auto flag = i > 0 && i + 1u < sa.size() && x == prev_x && x + 1 < sa[i + 1];
            prev_x = x + flag;
            return std::pair<int64_t, int64_t>(x + flag, i);
        };
        return make_segmentation_mod(sa.size(), eps, in_fun_index, out_fun, end_cond, connected);
    }
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
void write_segments_count(string out_fn, const T &sa, bool is_dict_point_seq = false, bool use_all_epsilons = false, 
                bool is_int = false) {
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
//              << std::setw(20) << "#segments_connected"
            //   << std::endl;

    #pragma omp parallel for ordered schedule(dynamic)
    for (size_t i = 0; i < eps_vec.size(); i++) {
        uint64_t count_disconnected;
        if(!is_int)
            count_disconnected = get_segments_count(eps_vec[i], sa, is_dict_point_seq, false);
            //    auto count_connected = get_segments_count(eps_vec[i], sa, false, true);
        else
            count_disconnected = get_segments_count_int(eps_vec[i], sa, is_dict_point_seq, false);
    

        #pragma omp ordered
        {
            // std::cout << std::setw(10) << eps_vec[i]
            //           << std::setw(20) << count_disconnected
                    //  << std::setw(20) << count_connected
                    //   << std::endl;
            if (i == 0)
                out_file << get_unique_elements(sa) << std::endl;
            out_file << eps_vec[i] << " " << count_disconnected << std::endl;
            // out_file << eps_vec[i] << " " << count_connected << std::endl;
        }
    }
    out_file.close();
}

inline void process_suffix_array(std::string gn_fn, std::string sa_fn, int64_t kmer_size,
        std::string out_fn, bool is_dict_point_seq, bool use_all_epsilons) {
    suffix_array<int64_t> sa(gn_fn, sa_fn, kmer_size);
    write_segments_count(out_fn, sa, is_dict_point_seq, use_all_epsilons, false);
}

inline void process_fracminhash_vec(std::string sketch_fn, 
        std::string out_fn, bool is_dict_point_seq, bool use_all_epsilons){
    std::vector<int64_t> kmer_sketch_vec;
    double theta;
    ifstream is(sketch_fn, std::ios::binary);
    essentials::load_pod(is, theta);
    essentials::load_vec(is, kmer_sketch_vec);
    write_segments_count(out_fn, kmer_sketch_vec, is_dict_point_seq, use_all_epsilons, false);
}

inline void process_int_vec(std::string int_vec_fn, 
        std::string out_fn, bool use_all_epsilons){
    std::vector<uint64_t> kmer_vec;
    ifstream is(int_vec_fn, std::ios::binary);
    essentials::load_vec(is, kmer_vec);
    write_segments_count(out_fn, kmer_vec, false, use_all_epsilons, true);
}

int main(int argc, char **argv) {
    auto opt = parse_cla_count_segments(argc, argv);

    bool is_dict_point_seq = opt.is_dict_point_seq;
    bool is_sketch_data = opt.is_sketch_data;
    bool is_genome_data = opt.is_genome_data;
    bool use_all_epsilons = opt.use_all_epsilons;
    bool is_int_vec = opt.is_int_vec;

    string gn_fn, sa_fn, in_fn, out_fn;
    int64_t kmer_size;
    if(is_genome_data){
        gn_fn = opt.gn_fn;
        sa_fn = opt.sa_fn;
        kmer_size = opt.kmer_size;
        out_fn = opt.count_fn;
    }
    else if(is_sketch_data){
        in_fn = opt.in_fn;
        out_fn = opt.count_fn;
    }
    else if(is_int_vec){
        in_fn = opt.in_fn;
        out_fn = opt.count_fn;
    }
    

    auto start = std::chrono::high_resolution_clock::now();
    if(is_genome_data)
        process_suffix_array(gn_fn, sa_fn, kmer_size, out_fn, is_dict_point_seq, use_all_epsilons);
    else if(is_sketch_data)
        process_fracminhash_vec(in_fn, out_fn, true, use_all_epsilons);
    else if(is_int_vec)
        process_int_vec(in_fn, out_fn, use_all_epsilons);
    

    auto end = std::chrono::high_resolution_clock::now();
    auto duration_s = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    std::cout << "Completed in " << duration_s << " seconds." << std::endl;

    return 0;
}
