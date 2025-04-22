// This file is part of PGM-index <https://github.com/gvinciguerra/PGM-index>.
// Copyright (c) 2018 Giorgio Vinciguerra.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <stdexcept>
#include <type_traits>

template<typename T>
using LargeSigned = typename std::conditional_t<std::is_floating_point_v<T>,
                                                long double,
                                                std::conditional_t<(sizeof(T) < 8), int64_t, __int128>>;

template<typename X, typename Y>
class OptimalPiecewiseLinearModel {
private:
    using SX = LargeSigned<X>;
    using SY = LargeSigned<Y>;

    struct Slope {
        SX dx{};
        SY dy{};

        bool operator<(const Slope &p) const { return dy * p.dx < dx * p.dy; }
        bool operator>(const Slope &p) const { return dy * p.dx > dx * p.dy; }
        bool operator==(const Slope &p) const { return dy * p.dx == dx * p.dy; }
        bool operator!=(const Slope &p) const { return dy * p.dx != dx * p.dy; }
        explicit operator long double() const { return dy / (long double) dx; }
    };

    struct Point {
        X x{};
        Y y{};

        Slope operator-(const Point &p) const { return {SX(x) - p.x, SY(y) - p.y}; }
    };

    const Y epsilon;
    std::vector<Point> lower;
    std::vector<Point> upper;
    X first_x = 0;
    X last_x = 0;    
    size_t lower_start = 0;
    size_t upper_start = 0;
    size_t points_in_hull = 0;
    Point rectangle[4];

    Y act_sa_indx;
    Y repeat_count;

    auto cross(const Point &O, const Point &A, const Point &B) const {
        auto OA = A - O;
        auto OB = B - O;
        // std::cout<<OA.dx <<" "<< OB.dy<<" "<<OA.dy <<" "<< OB.dx<<std::endl;
        // std::cout<<OA.dx * OB.dy<<" "<<OA.dy * OB.dx<<std::endl;
        return (OA.dx * OB.dy) - (OA.dy * OB.dx);
    }

public:
    uint64_t seg_start_indx;
    class CanonicalSegment;
    /**
     * The explicit keyword in C++ is primarily used in constructors to indicate that the constructor 
     * should not be implicitly used for implicit type conversions. When you mark a constructor with 
     * explicit, you're telling the compiler not to perform automatic type conversions using that constructor.
    */
    explicit OptimalPiecewiseLinearModel(Y epsilon) : epsilon(epsilon), lower(), upper() {
        if (epsilon < 0)
            throw std::invalid_argument("epsilon cannot be negative");

        upper.reserve(1u << 8);
        lower.reserve(1u << 8);
    }

    bool add_point(const X &x, const Y &y, Y occurenceCount = 1) {
        if (points_in_hull > 0 && x <= last_x){
            // std::std::cout<<x<< " "<<last_x<<std::endl;
            throw std::logic_error("Points must be increasing by x.: "
                +std::to_string(uint64_t(x))+" "+std::to_string(uint64_t(last_x)));
        }
            

        occurenceCount = (int64_t)1;
        // std::cout<<"curr point: "<<x<<std::endl;

        last_x = x;
        auto max_y = std::numeric_limits<Y>::max();
        auto min_y = std::numeric_limits<Y>::lowest();
        Point p1{x, y >= max_y - epsilon - occurenceCount + 1 ? max_y : y + epsilon + occurenceCount - 1};
        Point p2{x, y <= min_y + epsilon ? min_y : y - epsilon};

        // cout<<"("<<x<<", "<<int64_t(y) + epsilon + occurenceCount -1<<") , ("<<x<<", "<<int64_t(y) - epsilon<<")"<<endl;
        // if(x >= 2188 && x<= 2334){
        //     std::cout<<"["<<(int64_t)SY(y) + epsilon + occurenceCount -1<<", ";
        //     std::cout<<(int64_t)SY(y) - epsilon<<"]"<<std::endl;
        // }
        if (points_in_hull == 0) {
            first_x = x;
            act_sa_indx = y;
            repeat_count = occurenceCount;
            rectangle[0] = p1;
            rectangle[1] = p2;
            upper.clear();
            lower.clear();
            upper.push_back(p1);
            lower.push_back(p2);
            upper_start = lower_start = 0;
            ++points_in_hull;
            return true;
        }

        if (points_in_hull == 1) {
            rectangle[2] = p2;
            rectangle[3] = p1;
            upper.push_back(p1);
            lower.push_back(p2);
            ++points_in_hull;
            return true;
        }

//        if (epsilon == 0) {
//            auto p1_on_line1 = p1 - rectangle[0] == rectangle[2] - rectangle[0];
//            points_in_hull = p1_on_line1 ? points_in_hull + 1 : 0;
//            last_x = x;
//            return p1_on_line1;
//        }

        auto slope1 = rectangle[2] - rectangle[0];
        auto slope2 = rectangle[3] - rectangle[1];
        bool outside_line1 = p1 - rectangle[2] < slope1; // minus has higher precedence
        bool outside_line2 = p2 - rectangle[3] > slope2;

        if (outside_line1 || outside_line2) {
            points_in_hull = 0;
            return false;
        }

        if (p1 - rectangle[1] < slope2) {
            // Find extreme slope
            auto min = lower[lower_start] - p1;
            auto min_i = lower_start;
            for (auto i = lower_start + 1; i < lower.size(); i++) {
                auto val = lower[i] - p1;
                if (val > min)
                    break;
                min = val;
                min_i = i;
            }

            rectangle[1] = lower[min_i];
            rectangle[3] = p1;
            lower_start = min_i;

            // Hull update
            auto end = upper.size();
            for (; end >= upper_start + 2 && cross(upper[end - 2], upper[end - 1], p1) <= 0; --end)
                continue;
            upper.resize(end);
            upper.push_back(p1);
        }

        if (p2 - rectangle[0] > slope1) {
            // Find extreme slope
            auto max = upper[upper_start] - p2;
            auto max_i = upper_start;
            for (auto i = upper_start + 1; i < upper.size(); i++) {
                auto val = upper[i] - p2;
                if (val < max)
                    break;
                max = val;
                max_i = i;
            }

            rectangle[0] = upper[max_i];
            rectangle[2] = p2;
            upper_start = max_i;

            // Hull update
            auto end = lower.size();
            for (; end >= lower_start + 2 && cross(lower[end - 2], lower[end - 1], p2) >= 0; --end)
                continue;
            lower.resize(end);
            lower.push_back(p2);
        }

        ++points_in_hull;
        // for(int i=0; i<4; i++){
        //     std::cout<<rectangle[i].x<<" "<<(int64_t)rectangle[i].y<<std::endl;
        // }
        // for(int i=0; i<lower.; i++){
        //     std::cout<<lower[i].x<<" "<<(int64_t)lower[i].y<<std::endl;
        // }
        // for(int i=0; i<4; i++){
        //     std::cout<<lower[i].x<<" "<<(int64_t)lower[i].y<<std::endl;
        // }
        // std::cout<<std::endl;
        // if(x >= 2188 && x<=2334){
        //     std::cout<<x<<" "<<(int64_t) y<<std::endl;
        //     // sanityCheck(x);
        // }
        // for(int i=0; i<4; i++){
        //     if(rectangle[i].x == 65 || rectangle[i].x == 155375131957){
        //         std::cout<<rectangle[i].x<<" "<<(int64_t)rectangle[i].y<<std::endl;
        //         break;
        //     }
        // }

        return true;
    }

    CanonicalSegment get_segment() const {
        if (points_in_hull == 1)
            return CanonicalSegment(rectangle[0], rectangle[1], first_x, last_x, act_sa_indx, repeat_count);
        return CanonicalSegment(rectangle, first_x, last_x, act_sa_indx, repeat_count);
    }

    void reset() {
        points_in_hull = 0;
        lower.clear();
        upper.clear();
    }
};

template<typename X, typename Y>
class OptimalPiecewiseLinearModel<X, Y>::CanonicalSegment {
    friend class OptimalPiecewiseLinearModel;

    Point rectangle[4];
    X first, last;

    Y act_sa_indx, rep_count;

    CanonicalSegment(const Point &p0, const Point &p1, X first, X last, Y act_sa_indx, Y rep_count) 
        : rectangle{p0, p1, p0, p1}, first(first), last(last),
        act_sa_indx(act_sa_indx), rep_count(rep_count)
        {};

    CanonicalSegment(const Point (&rectangle)[4], X first, X last, Y act_sa_indx, Y rep_count)
        : rectangle{rectangle[0], rectangle[1], rectangle[2], rectangle[3]}, 
        first(first), last(last),
        act_sa_indx(act_sa_indx), rep_count(rep_count)
        {};

    bool one_point() const {
        return rectangle[0].x == rectangle[2].x && rectangle[0].y == rectangle[2].y
            && rectangle[1].x == rectangle[3].x && rectangle[1].y == rectangle[3].y;
    }

public:

    CanonicalSegment() = default;

    explicit CanonicalSegment(X first) : CanonicalSegment({first, 0}, {first, 0}, first) {};

    X get_first_x() const { return first; }
    X get_last_x() const { return last; }

    std::tuple<int64_t, int64_t> get_sa_rep(){
        return {act_sa_indx, rep_count};
    }

    CanonicalSegment copy(X x) const {
        auto c(*this);
        c.first = x;
        return c;
    }

    std::tuple<int64_t, uint8_t, SY> get_fixed_point_segment(X origin, X max_input) const {
        if (one_point())
            return {0, 0, (rectangle[0].y + rectangle[1].y) / 2};
        // std::cout<<"in get_fixed_ps: "<<max_input<<std::endl;
        auto &p1 = rectangle[1];
        auto max_slope = rectangle[3] - rectangle[1];

        auto is_slope_integral = max_slope.dy % max_slope.dx == 0;
        auto slope_exponent = is_slope_integral ? 0 : (uint8_t) std::ceil(std::log2(max_input)) + 1;
        auto slope_significand = (max_slope.dy << slope_exponent) / max_slope.dx;

        auto intercept_n = max_slope.dy * (SX(origin) - p1.x);
        auto intercept_d = max_slope.dx;
        auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
        auto intercept = (intercept_n + rounding_term) / intercept_d + p1.y;

        return {slope_significand, slope_exponent, intercept};
    }
    
    SY get_prediction(X kmer) const{
        auto &p1 = rectangle[1];
        auto max_slope = rectangle[3] - rectangle[1];
        auto origin = first;
        auto max_input = last-first+1;

        auto is_slope_integral = max_slope.dy % max_slope.dx == 0;
        auto slope_exponent = is_slope_integral ? 0 : (uint8_t) std::ceil(std::log2(max_input)) + 1;
        auto slope_significand = (max_slope.dy << slope_exponent) / max_slope.dx;

        auto intercept_n = max_slope.dy * (SX(origin) - p1.x);
        auto intercept_d = max_slope.dx;
        auto rounding_term = ((intercept_n < 0) ^ (intercept_d < 0) ? -1 : +1) * intercept_d / 2;
        auto intercept = (intercept_n + rounding_term) / intercept_d + p1.y;

        auto prediction = (((slope_significand) * (kmer - first)) >> slope_exponent) + intercept;
        return prediction;
    }

    std::string ConvertToString(__int128 num) {
        std::string str;
        do {
            int digit = num % 10;
            str = std::to_string(digit) + str;
            num = (num - digit) / 10;
        } while (num != 0);
        return str;
    }

    
    std::tuple<int64_t, int64_t> get_knot_intersection(){
        if(one_point()){
            int64_t knot_si = rectangle[1].y;
            int64_t knot_ei = rectangle[2].y;
        }
        X knot_start = first;
        X knot_end = last;
        auto max_slope = rectangle[3] - rectangle[1];
        int64_t knot_si = round(int64_t(rectangle[1].y) +
                    (max_slope.dy*(__int128_t(knot_start)-__int128_t(rectangle[1].x)))
                    /(double)max_slope.dx);
        int64_t knot_ei = round(int64_t(rectangle[1].y)+
                    (uint64_t(max_slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)max_slope.dx));
        return {knot_si, knot_ei};
    }

    bool isSlopeAtHalfPoint(){
        X knot_start = first;
        X knot_end = last;
        auto max_slope = rectangle[3] - rectangle[1];
        double knot_si = int64_t(rectangle[1].y) +
                    (max_slope.dy*(__int128_t(knot_start)-__int128_t(rectangle[1].x)))
                    /(double)max_slope.dx;
        double knot_ei = int64_t(rectangle[1].y)+
                    (uint64_t(max_slope.dy))*((__int128_t(knot_end)-__int128_t(rectangle[1].x))
                    /(double)max_slope.dx);
        
        // std::cout<<"si: "<<knot_si<<" ei: "<<knot_ei
        //     <<" bool: "
        //     <<(floor(knot_si-0.5)== floor(knot_si)) <<" "<<(floor(knot_ei - 0.5)==floor(knot_ei))
        //     <<" fsi: "<<floor(knot_si-0.5) <<" fei: "<<floor(knot_ei - 0.5)
        //     <<" full: "
        //     <<((floor(knot_si-0.5)== floor(knot_si))&& (floor(knot_ei - 0.5)==floor(knot_ei)))
        //     <<std::endl;
    
        return ((floor(knot_si-0.5)== floor(knot_si))&& (floor(knot_ei - 0.5)==floor(knot_ei)));
    }

    int64_t Pred_idx(int64_t brk_beg_sa_indx, int64_t brk_end_sa_indx,
                uint64_t query_val, uint64_t brk_beg_kval, 
                uint64_t brk_end_kval){
        return round(brk_beg_sa_indx + 
                    (brk_end_sa_indx - brk_beg_sa_indx)*
                    ((uint64_t(query_val) - uint64_t(brk_beg_kval)) /
                    (double)(brk_end_kval - brk_beg_kval)) );
    }
    
};


template<typename Fin, typename Fout, typename Fec>
size_t make_segmentation_mod(size_t n, size_t epsilon, Fin in, Fout out,
    Fec end_condition, bool connect) {
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t c = 0;
    if(end_condition()){
         throw std::invalid_argument("Invalid input points");
    }
    auto p = in(0);

    OptimalPiecewiseLinearModel<X, Y> opt(epsilon);
    opt.add_point(p.first, p.second);

    for (size_t i = 1; i < n; ++i) {
        if(end_condition()) break;
        auto next_p = in(i);
        if (next_p.first == p.first)
            continue;
        if (!opt.add_point(next_p.first, next_p.second)) {
            out(opt.get_segment());
            if (connect) {
                opt.add_point(p.first, p.second);
                opt.add_point(next_p.first, next_p.second);
            }
            else {
                opt.add_point(next_p.first, next_p.second);
            }
            ++c;
        }
        p = next_p;
    }

    out(opt.get_segment());
    return ++c;
}

template<typename Fin, typename Fout>
size_t make_segmentation(int64_t n, int64_t epsilon, Fin in, Fout out) {
    
    if (n == 0)
        return 0;

    using X = typename std::invoke_result_t<Fin, size_t>::first_type;
    using Y = typename std::invoke_result_t<Fin, size_t>::second_type;
    size_t c = 0;// count variable
    size_t start = 0;
    
    bool isIncluded = false;
    OptimalPiecewiseLinearModel<X, Y> opt(epsilon);
    
    auto p = in(start);
    
    while(1){
        p = in(start); // pair of (x,y)
        if(p.first != -1){
            break;
        }
        start++;
    }
    opt.seg_start_indx = start;
    
    int64_t curr_sa_indx = start, occurenceCount = 1;
    while(1){
        curr_sa_indx++;
        auto temp = in(curr_sa_indx);
        if(temp.first == p.first) occurenceCount++;
        else break;
    }
    
    start = curr_sa_indx; 
    opt.add_point(p.first, p.second+epsilon, 1);
    
    for (int64_t i = start; i < n; ++i) {  
        auto next_p = in(i);
        if ((i != start && next_p.first == p.first) || next_p.first == -1)
            continue;
    
        p = next_p;
        curr_sa_indx = i;
        occurenceCount = 1;
        while(1){
            curr_sa_indx++;
            if(curr_sa_indx == n) break;
            auto temp = in(curr_sa_indx);
            if(temp.first == p.first) occurenceCount++;
            else break;
        }
        
        i+=(occurenceCount-1);
        
        isIncluded = opt.add_point(p.first, p.second+epsilon, occurenceCount);
        if (!isIncluded) {
            OptimalPiecewiseLinearModel<int64_t, uint64_t>::CanonicalSegment cs = opt.get_segment();
            out(opt.get_segment());
            start = i - occurenceCount;
            opt.seg_start_indx = start+1; 
            i = start++;                   
            ++c;
        }
    }
    out(opt.get_segment());
    ++c;
    return c;
}
