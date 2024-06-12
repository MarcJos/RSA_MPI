//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief Types for storage of vectors etc.

#pragma once


#include <cmath>
#include<vector>
#include<array>

#include "link_to_auxi.hxx"

//! long vector storages

using vec_double = std::vector<double>;

using  vec_int = std::vector<uint64_t>;


//! D-DIM array storages
template<int DIM>
using vec_d = sac_de_billes::Point<DIM>;

template<int DIM>
using  arr_vec_double = std::vector<vec_d<DIM>>;

template<int DIM>
using vec_i = std::array<int, DIM>;

//! sphere

template<int DIM>
using rsa_sphere = sac_de_billes::Sphere<DIM>;


class BooleanVector : public std::vector<char> {
public:
    BooleanVector() {}
    BooleanVector(uint64_t size_) : std::vector<char>(size_) {}
    BooleanVector(uint64_t size_, bool val_) : std::vector<char>(size_, val_) {}
};
