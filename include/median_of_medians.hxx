//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
//! @brief : find the indices of the n first items of a vector distributed accross MPI processes
//!

#include "link_to_auxi.hxx"

namespace median_of_medians {

//! @return  find a value "pivot" such that the vector v=a_vector_of_values, supposedly distributed accross MPI processes
//! is split in a part v[i] < pivot of size a_nb_first_elements, and a part v[i] >= pivot
//! @param a_vector_of_values : a vector containing values that can be converted to double
//! @param a_nb_first_elements : desired number of elements
//! @param order : order < or >
template<class TypeVector, class Order>
double find_pivot_for_first_elements(const TypeVector& a_vector_of_values,
    int64_t a_nb_first_elements, const Order& order);

namespace auxi {
//! @brief find a value "pivot" such that :
//! the percentage of (values < pivot) is < a_percentage
//! the percentage of (values >= pivot) is > a_percentage
//! @tparam TypeVector
//! @param a_vector_of_values_sorted assumed to be sorted
//! @param a_percentage
//! @return the value
template<class TypeVector>
std::tuple<double, double>  find_pivot(const TypeVector& a_vector_of_values_sorted, double a_percentage);
//! @return true if the proposed pivot is too small to obtain a_nb_first_elements smaller than it
//! @tparam TypeVector
//! @param a_vector_of_values_sorted
//! @param a_nb_first_elements
//! @param a_global_pivot
template<class TypeVector, class Order>
int64_t local_nb_smaller_than_pivot(const TypeVector& a_vector_of_values_sorted,
    double a_candidate_pivot, const Order& order);
template<class TypeVector, class Order>
int64_t global_nb_smaller_than_pivot(const TypeVector& a_vector_of_values_sorted,
    double a_candidate_pivot, const Order& order);
} // namespace  auxi

} // namespace  median_of_medians

#include <median_of_medians.ixx>
