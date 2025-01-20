//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <mpi.h>
#include <cassert>
#include <iostream>
#include <cstdint>

namespace rsa_mpi {


inline uint64_t compute_mpi_sum(uint64_t local_quantity) {
	uint64_t ret = 0;
	MPI_Allreduce(&local_quantity, &ret, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
	return ret;
}

inline int64_t compute_mpi_sum(int64_t local_quantity) {
	int64_t ret = 0;
	MPI_Allreduce(&local_quantity, &ret, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	return ret;
}

inline double compute_mpi_sum(double local_quantity) {
	double ret = 0;
	MPI_Allreduce(&local_quantity, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return ret;
}

template<class T>
inline int compute_mpi_sum(T) = delete;

inline double compute_mpi_max(double local_quantity) {
	double ret = 0;
	MPI_Allreduce(&local_quantity, &ret, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	return ret;
}

inline double compute_mpi_average(double local_quantity) {
	double average_quantity = 0.;
	auto comm = MPI_COMM_WORLD;
	int mpi_size = 0;
	MPI_Comm_size(comm, &mpi_size);
	MPI_Allreduce(&local_quantity, &average_quantity, 1, MPI_DOUBLE, MPI_SUM, comm);
	average_quantity *= 1. / mpi_size;
	return average_quantity;
}

inline int get_number_of_mpi_processes() {
	auto comm = MPI_COMM_WORLD;
	int mpi_size = 0;
	MPI_Comm_size(comm, &mpi_size);
	return mpi_size;
}

/** @brief This function intializes the variable to -1 if it's the first call, otherwise it return my_rank. The default value is -1 to force the user to set it to the mpi rank with the function set_my_rank. */
inline int get_my_rank() {
	int my_rank = -1;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	return my_rank;
}

/** @brief This functions return true if the current mpi process is the master rank, ie the mpi_rank == 0. */
inline bool is_master_rank() {
	constexpr int master_rank = 0;
	bool ret = (master_rank == get_my_rank());
	return ret;
}


/** @brief This function returns the result of a boolean test only made on the master processor **/
template<class TEST>
inline bool test_only_on_master(TEST test) {
	bool result = false;
	if (is_master_rank()) {
		result = test();
	}
	MPI_Request request;
	MPI_Ibcast(&result, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD, &request);
	return result;
}

/** @brief displays a parameter according to its type.*/
template<typename Arg>
void message(Arg a_arg) {
	if (is_master_rank()) {
		std::cout << a_arg << std::endl;
	}
}

/** @brief displays a range of parameters according to their types. A space is added between two expressions.*/
template<typename Arg, typename... Args>
void message(Arg a_arg, Args... a_args) {
	if (is_master_rank()) {
		std::cout << a_arg << " ";
		message(a_args...);
	}
}

inline void banner() {
	message("================= ===== ================= ");
	message("==== rsa_algo with mpi parallelization == ");
	message("================= start ================= ");
}

inline void bye_bye_message() {
	message("=========== end ===========");
}

/** @brief init function initializes MPI, sets the mpi rank to the static variable my_rank and displays the banner. */
template<typename... Args>
void init(Args... a_args) {
	assert(sizeof...(Args) == 2);
	MPI_Init(a_args...);
	banner();
}

inline void finalize() {
	bye_bye_message();
	[[maybe_unused]] auto success = MPI_Finalize();
	assert(success == MPI_SUCCESS);
}
} //  namespace rsa_mpi
