//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

namespace rsa_helper {

template<int DIM>
template<typename... Args>
void helper_setter<DIM>::operator()(const int a_idx, arr_vec_double<DIM>& a_to, const double a_arg, const Args... a_from) {
    auxi_operator<DIM>(a_idx, a_to, a_arg, a_from...);
}

template<int DIM>
template<int D, typename... Args>
void helper_setter<DIM>::auxi_operator(const int a_idx, arr_vec_double<DIM>& a_to, const double a_arg, const Args... a_from) {
    constexpr int d = DIM - (sizeof...(Args) + 1);
    static_assert(D == DIM - d);
    a_to[a_idx][d] = a_arg;
    if constexpr (D > 1) {
        this->auxi_operator<D - 1>(a_idx, a_to, a_from...);
    }
}

template<int DIM>
void helper_setter<DIM>::operator()(const int a_idx, arr_vec_double<DIM>& a_to, const vec_d<DIM>& a_from) {
    a_to[a_idx] = a_from;
}

template<int DIM>
bool helper_check_interval<DIM>::operator()(const std::array<double, DIM>& a_lower,
    const std::array<double, DIM>& a_upper,
    const vec_d<DIM>& a_position) const {
    for (int d = 0; d < DIM; d++) {
        if (a_position[d] < a_lower[d]) return false;
        if (a_position[d] > a_upper[d]) return false;
    }
    return true;
}

template<int DIM>
void helper_extractor<DIM>::operator()(const arr_vec_double<DIM>& a_from, const int a_idx, std::array<double, DIM>& a_to) {
    a_to = a_from[a_idx];
}

} // namespace rsa_helper 
