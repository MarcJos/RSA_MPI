//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
/**
 * \file helper.hxx
 * \brief This file contains helper structures. These structures are designed to deal with the DIM template.
 * \author Raphael Prat
 * \version 0.1
 * \date 22 avril 2023
 */

#include <array>
#include <vector>

#include <basic_types.hxx>

namespace rsa_helper {
//! brief : This functor sets the values of args in a_to for the index a_idx
template<int DIM>
struct helper_setter {
	template<typename... Args>
	void operator()(const int a_idx, arr_vec_double<DIM>& a_to, const double a_arg, const Args... a_from);
	void operator()(const int a_idx, arr_vec_double<DIM>& a_to, const vec_d<DIM>& a_from);
private:
	template<int D, typename... Args>
	void auxi_operator(const int a_idx, arr_vec_double<DIM>& a_to, const double a_arg, const Args... a_from);
};

//! @brief  : This functor checks if args are between a_lower and a_upper for every DIMensions
template<int DIM>
struct helper_check_interval {
	bool operator()(const std::array<double, DIM>& a_lower,
		const std::array<double, DIM>& a_upper,
		const vec_d<DIM>& position) const;
};

template<int DIM>
class rsa_data_storage;

//! @brief  : This functor gets the values of a_from for the index a_idx and it stores thme in a_to
template<int DIM>
struct helper_extractor {
	void operator()(const arr_vec_double<DIM>& a_from, const int a_idx, std::array<double, DIM>& a_to);
};

} // namespace rsa_helper


#include<helper.ixx>
