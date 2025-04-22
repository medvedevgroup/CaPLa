#pragma once
#include <vector>
#include <string>
#include "utils/util.hpp"

using namespace std;

template <typename T>
class suffix_array {
private:
    std::vector<T> sa;
    std::vector<char> str;
    size_t fSize;
    int64_t kmer_size;
    uint64_t num_uniq_kmers;

public:
    class Iterator {
    private:
        typename std::vector<T>::iterator it;
        suffix_array& sa;

    public:
        Iterator(typename std::vector<T>::iterator iter, suffix_array& sa) 
                : it(iter), sa(sa) {}

        // Overload prefix increment operator (++it)
        Iterator& operator++() {
            ++it;
            return *this;
        }

        // Overload postfix increment operator (it++)
        Iterator operator++(int) {
            Iterator temp = *this;
            ++it;
            return temp;
        }

        // Overload prefix decrement operator (++it)
        Iterator& operator--() {
            --it;
            return *this;
        }

        // Overload postfix increment operator (it++)
        Iterator operator--(int) {
            Iterator temp = *this;
            --it;
            return temp;
        }

        // Overload -= operator (it -= n)
        Iterator& operator-=(std::ptrdiff_t n) {
            it -= n;
            return *this;
        }

        // Overload equality operator (it1 == it2)
        bool operator==(const Iterator& other) const {
            return it == other.it;
        }

        // Overload inequality operator (it1 != it2)
        bool operator!=(const Iterator& other) const {
            return it != other.it;
        }

        // Overload - operator (it1 - it2)
        std::ptrdiff_t operator-(const Iterator& other) const {
            return it - other.it;
        }

        // Optionally, overload dereference operator (*it)
        T operator*() {
            // std::cout<<"iterator val: "<<*it
            // <<" "<<sa.size()<<std::endl;

            return sa.get_kmer_val_at_str_index(*it);
            // return *it;
        }

        T operator[](std::ptrdiff_t index) const{
            return sa.get_kmer_val_at_str_index(*(it + index));
            // return *(it + index);
        }
    };


    suffix_array(){}

    suffix_array(const std::string &gn_fn, const std::string &sa_fn, int kmer_size):
        kmer_size(kmer_size)
    {
        load(gn_fn, sa_fn);
        num_uniq_kmers = 0;
    }

    // Overloading the [] operator
    T operator[](std::size_t index) const{
        if (index < 0 ||  index >= sa.size()) {
            throw std::out_of_range("Index out of range for SA");
        }
        return get_kmer_val_at_sa_index(index);
    }

    Iterator begin() {
        return Iterator(sa.begin(), *this);
    }

    Iterator end() {
        return Iterator(sa.end(), *this);
    }

    // Additional member functions for the vector class
    void push_back(const T& value) {
        sa.push_back(value);
    }

    std::size_t size() const {
        return sa.size();
    }

    uint64_t get_num_uniq_kmers() {
        if (num_uniq_kmers) return num_uniq_kmers;
        int64_t curr_kmer, prev_kmer = -1;
        for (int i = 0; i < sa.size(); i++){
            curr_kmer = get_kmer_val_at_sa_index(i);
            if (curr_kmer == prev_kmer || curr_kmer == -1) continue;
            num_uniq_kmers++;
            prev_kmer = curr_kmer;
        }
        return num_uniq_kmers;
    }

    int64_t get_kmer_val_at_sa_index(const int64_t i) const {
        return get_kmer_val_at_str_index(sa[i]);
    }

    int64_t get_kmer_val_at_str_index(int64_t i) const {
        auto kmer_start = str.begin() + i;
        auto kmer_end = kmer_start + kmer_size;
        if (kmer_end > str.end()) {
            return -1;
        }

        int64_t result = 0;
        for (auto it = kmer_start; it != kmer_end; ++it) {
            result <<= 2;
            switch (*it) {
                case 'A':
                    break;
                case 'C':
                    result |= 1;
                    break;
                case 'G':
                    result |= 2;
                    break;
                case 'T':
                    result |= 3;
                    break;
                default:
                    throw std::invalid_argument("Invalid character in sequence");
            }
        }
        return result;
    }

    const int64_t get_largest() const {
        int64_t last_idx = fSize-1;
        while (get_kmer_val_at_sa_index(last_idx) == -1) {
            last_idx--;
        }
        return get_kmer_val_at_sa_index(last_idx);
    }

    inline const T get_sa_val(uint64_t idx) const {
        return sa[idx];
    }

    inline const int64_t get_kmer_val(const string &s) const {
        // if we have a smaller suffix than k-mer size,
        // we want to ignore it so that it never becomes a breakpoint
        if(s.compare("-1") == 0) return -1; 
        int64_t  val=0;
        int i=0;
        while(s[i]){
            val=val<<2;
            if(s[i]=='A'){
                val=val|0;
            }
            else if(s[i]=='C'){
                val=val|1;
            }
            else if(s[i]=='G'){
                val=val|2;
            }
            else{
                val=val|3;
            }
            i++;
        }
        return val;
    }

    const std::string get_kmer_at_index(const int64_t i) const{
        int64_t suff_i_start = sa[i];
        if(str.begin()+suff_i_start+kmer_size <= str.end()){
            std::string suffix_i(str.begin()+suff_i_start, str.begin()+suff_i_start+kmer_size);
            return suffix_i;
        }
        return "-1"; // substring size < k: ignore these
    }
    
    const int64_t binary_search(const std::string &s, int64_t lo, int64_t hi, bool isFirstIndexReturned) const
    {
        int64_t mid, idx, nLcp, gLcp, loLcp = 0, hiLcp = 0;
        char ch = ' ';
        while(1){
            // nLcp = min(loLcp, hiLcp);  
            // cout<<lo<<" "<<(lo+hi)/2<<" "<<hi<<endl;
            // std::cout<<get_kmer_val_at_sa_indx(lo)<<" "
            // <<get_kmer_val(s)<<" "
            // <<get_kmer_val_at_sa_indx(hi)<<endl;
            nLcp = loLcp < hiLcp ? loLcp : hiLcp;      
            if(hi - lo <= 1){
                for(int64_t i=lo; i<=hi; i++){
                    // sa_bs_cnt++;
                    idx = sa[i];
                    gLcp = nLcp;
                    for(; gLcp<kmer_size && idx+gLcp<fSize; gLcp++ ) {
                        ch = str[idx+gLcp];
                        if(s[gLcp] != ch) break;
                    }
                    if(gLcp == kmer_size) {
                        if(isFirstIndexReturned) return i;
                        return idx;
                    }
                }
                return -1; // not found
            }
            mid = (lo + hi) >> 1;
            idx = sa[mid]; 
            for(; nLcp<kmer_size && idx+nLcp<fSize; nLcp++ ) {
                ch = str[idx+nLcp];
                if(s[nLcp] != ch) break;
            }
            if(nLcp == kmer_size) {
                if(isFirstIndexReturned) return mid;
                return idx;
            }
            else if(nLcp + idx == fSize || s[nLcp] > ch)
            {
                lo = mid;
                loLcp = nLcp;
            }
            else
            {
                hi = mid;
                hiLcp = nLcp;
            }   
        }
    }

    std::string get_kmer_at_sa_index_it(int64_t suff_i_start){
        // int64_t suff_i_start = sa[i];
        if(str.begin()+suff_i_start+kmer_size <= str.end()){
            std::string suffix_i(str.begin()+suff_i_start, str.begin()+suff_i_start+kmer_size);
            return suffix_i;
        }
        return "-1"; // substring size < k: ignore these
    }

    std::string get_kmer_at_str_pos(size_t strPos){
        if(strPos + kmer_size >= str.size()){
            // cout<<strPos<<" " <<saq.str.size()<<endl;
            return "-1";
        }
        if(strPos + kmer_size < str.size()){
            std::string kmer(str.begin()+strPos, str.begin()+strPos+kmer_size);
            return kmer;
        }
        std::string kmer(str.begin()+strPos, str.end());
        return kmer;
    }   

    void load_genome(std::string genomeFile){
        str.reserve(fSize);
        string line;
        ifstream in_file(genomeFile.c_str());
        // skipping header
        getline(in_file, line);
        while(getline(in_file, line)){
            line = pla::util::trim_string(line);
            std::copy(line.begin(), line.end(), std::back_inserter(str));
        }
        fSize = str.size();
        printf("#Bases: %ld\n",fSize);
        // str.resize(fSize);
        // FILE *fp = fopen(genomeFile.c_str(), "rb");        
        // int err = fread(&str[0], 1, fSize, fp);
        // fclose(fp);
    }

    void load_sa(std::string sa_file){
        sa.resize(fSize);
        FILE *fp = fopen(sa_file.c_str(), "rb");
        int err = fread(&sa[0], sizeof(int64_t), fSize, fp);
        fclose(fp);
    }

    inline bool is_kmer_size_needed(){return true;}

    void set_kmer_size(int64_t km_size){
        kmer_size = km_size;
    }

    inline int64_t get_kmer_size(){return kmer_size;}

    void load(std::string genomeF, std::string saF){
        std::ifstream in(genomeF.c_str(), std::ifstream::ate | std::ifstream::binary);
        fSize = in.tellg();
        this->load_genome(genomeF);        
        this->load_sa(saF);
    }
};