//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <rsa_decoration.hxx>
#include <limits>

template<class TypeVector, class Order>
double median_of_medians::find_pivot_for_first_elements(const TypeVector& a_vector_of_values,
    int64_t a_nb_first_elements, const Order& order) {
    int64_t local_size = a_vector_of_values.size();
    int64_t global_size = rsa_mpi::compute_mpi_sum(local_size);

    // Case errors
    if (a_nb_first_elements > global_size) {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        throw runtime_error("Impossible to find the pivot for N first elements if N is smaller than the size of the array");
    }
    if (a_nb_first_elements <= 0) {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        throw runtime_error("Impossible to find the pivot for N first elements if N  <= 0");
    }
    if (global_size == 0) {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        throw runtime_error("Impossible to find the pivot for N first elements if no element");
    }
    // Case errors

    double average_percentage = a_nb_first_elements / (1. * global_size);
    //
    auto vector_copy = a_vector_of_values;
    std::sort(vector_copy.begin(), vector_copy.end(), order);
    //
    auto [local_candidate_0, local_candidate_1] = auxi::find_pivot(vector_copy, average_percentage);
    //
    int mpi_size = rsa_mpi::get_number_of_mpi_processes();
    vector<double> all_candidates(2 * mpi_size);
    double* pointer_to_data = all_candidates.data();
    MPI_Allgather(&local_candidate_0, 1, MPI_DOUBLE, pointer_to_data,
        1, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&local_candidate_1, 1, MPI_DOUBLE, pointer_to_data + mpi_size,
        1, MPI_DOUBLE, MPI_COMM_WORLD);
    std::sort(all_candidates.begin(), all_candidates.end(), order);
    all_candidates.erase(std::remove_if(all_candidates.begin(), all_candidates.end(), [](auto elem) {return std::isnan(elem);}), all_candidates.end());
    //
    auto compare_with_a_nb_first_elements = [&vector_copy, &order](int64_t nb_first_elements, auto a_candidate_pivot) {
        return nb_first_elements <= median_of_medians::auxi::global_nb_smaller_than_pivot(vector_copy,
            a_candidate_pivot, order);
        };
    auto pointer_to_pivot = std::upper_bound(all_candidates.begin(), all_candidates.end(),
        a_nb_first_elements, compare_with_a_nb_first_elements);
    // ERROR
    if (pointer_to_pivot == all_candidates.end()) {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        throw runtime_error("Unexpected - 0");
    }
    //
    //
    double global_pivot_max = *pointer_to_pivot;
    double global_pivot = global_pivot_max;
    int64_t nb_elements_below = median_of_medians::auxi::global_nb_smaller_than_pivot(vector_copy,
        global_pivot, order);
    if (nb_elements_below != a_nb_first_elements) {
        // dicotomy step
        // ERROR
        if (pointer_to_pivot == all_candidates.begin()) {
            // error message
            int64_t my_barrier = 0, end_barrier = 0;
            for (int j = 0; j < rsa_mpi::get_number_of_mpi_processes(); j++) {
                MPI_Barrier(MPI_COMM_WORLD);
                if (j == rsa_mpi::get_my_rank()) {
                    std::cerr << "####### BEGIN " << j << std::endl;
                    for (size_t i = 0; i < vector_copy.size(); i++) {
                        std::cerr << to_string(vector_copy[i]) << endl;;
                    }
                    std::cerr << "####### END " << j << std::endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            //
            MPI_Barrier(MPI_COMM_WORLD);
            rsa_mpi::message("####### BEGIN \n");
            for (size_t i = 0; i < all_candidates.size(); i++) {
                double a_candidate_pivot = all_candidates[i];
                rsa_mpi::message(to_string(all_candidates[i]) + " : " + to_string(median_of_medians::auxi::global_nb_smaller_than_pivot(vector_copy,
                    a_candidate_pivot, order)) + "\n");
            }
            rsa_mpi::message("####### END \n");
            MPI_Barrier(MPI_COMM_WORLD);
            //

            MPI_Barrier(MPI_COMM_WORLD);
            rsa_mpi::message(to_string(average_percentage));
            rsa_mpi::message(to_string(nb_elements_below) + " " + to_string(a_nb_first_elements));
            end_barrier = rsa_mpi::compute_mpi_sum(my_barrier);
            std::cerr << end_barrier;
            MPI_Barrier(MPI_COMM_WORLD);
            // error message
            rsa_mpi::message(__PRETTY_FUNCTION__);
            throw runtime_error("Unexpected - 1");
        }
        // ERROR
        else {
            double global_pivot_min = *(pointer_to_pivot - 1);
            double step = global_pivot_max - global_pivot_min;
            double local_add = 0.5 * step;
            int j_max = 20;
            for (int j = 0; j < j_max; j++) {
                global_pivot = global_pivot_min + local_add;
                nb_elements_below = median_of_medians::auxi::global_nb_smaller_than_pivot(vector_copy,
                    global_pivot, order);
                if (nb_elements_below == a_nb_first_elements) { break; }
                if (nb_elements_below > a_nb_first_elements) { local_add = 0.5 * std::abs(local_add); }
                if (nb_elements_below < a_nb_first_elements) { local_add = -0.5 * std::abs(local_add); }
                if (j == j_max) {
                    std::cerr << __PRETTY_FUNCTION__ << std::endl;
                    std::cerr << "Warning : I could not find the correct pivot" << std::endl;
                }
            }
        }
    }
    return global_pivot;
}




template<class TypeVector>
std::tuple<double, double> median_of_medians::auxi::find_pivot(const TypeVector& a_vector_of_values_sorted,
    double a_percentage) {
    double first = 0., secnd = 0.;
    if (a_vector_of_values_sorted.size() == 0) {
        first = std::nan("0");
        secnd = std::nan("0");
    } else {
        int64_t piv_index = static_cast<int64_t>(a_vector_of_values_sorted.size() * a_percentage) - 1;
        int64_t min_index = piv_index;
        int64_t max_index = piv_index + 1;
        if (abs(a_vector_of_values_sorted.size() * a_percentage - 1 - piv_index) < 1e-8) {
            min_index -= 1;
        }
        if (min_index < 0) {
            min_index = 0;
        }
        if (abs(a_vector_of_values_sorted.size() * a_percentage - 1 - piv_index - 1) < 1e-8) {
            max_index += 1;
        }
        if (max_index > int64_t(a_vector_of_values_sorted.size()) - 1) {
            max_index = a_vector_of_values_sorted.size() - 1;
        }
        //////
        first = a_vector_of_values_sorted[min_index] - std::numeric_limits<double>::epsilon();
        secnd = a_vector_of_values_sorted[max_index] + std::numeric_limits<double>::epsilon();
    }
    //
    return std::tuple<double, double>(first, secnd);
}

template<class TypeVector, class Order>
int64_t median_of_medians::auxi::global_nb_smaller_than_pivot(const TypeVector& a_vector_of_values_sorted,
    double a_candidate_pivot, const Order& order) {
    int64_t local_small = median_of_medians::auxi::local_nb_smaller_than_pivot(a_vector_of_values_sorted,
        a_candidate_pivot, order);
    int64_t global_small = rsa_mpi::compute_mpi_sum(local_small);
    return global_small;
}

template<class TypeVector, class Order>
int64_t median_of_medians::auxi::local_nb_smaller_than_pivot(const TypeVector& a_vector_of_values_sorted,
    double a_candidate_pivot, const Order& order) {
    if (a_vector_of_values_sorted.size() == 0) {
        return 0;
    }
    auto ptr_to_element = std::upper_bound(a_vector_of_values_sorted.begin(), a_vector_of_values_sorted.end(),
        a_candidate_pivot, order);
    int local_number = ptr_to_element - a_vector_of_values_sorted.begin();
    return local_number;
}
