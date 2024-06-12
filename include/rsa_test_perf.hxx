//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include<chrono>
#include <omp.h>

#include<test.hxx>
#include<rsa_parameters.hxx>

template<class TypeTest, class ...Args>
void make_XP(int nb_threads, TypeTest type_test, Args... args) {
    omp_set_num_threads(nb_threads);
#pragma omp parallel
    { // begin parallel section
#pragma omp master
        {
            rsa_mpi::message("Number of openmp threads = " + std::to_string(omp_get_num_threads()) + "\n");
        }
    } // end parallel section
    auto start = std::chrono::system_clock::now();
    type_test(args...);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    rsa_mpi::message("Time for " + std::to_string(nb_threads) + " threads :  " + std::to_string(elapsed_seconds.count()) + "s");
}

template<class TYPE>
void write_to_file(const std::string& file_path, TYPE value) {
    std::ofstream file(file_path, std::ios::out | std::ios::app);
    if (file.is_open()) {
        file << value;
        file.close();
    }
}

template<int DIM, int method>
void test_perf(double rad, vector<double> l_min_, vector<double> l_max_, size_t seed, size_t size, size_t n_ndraw, bool check_vox_time = false, bool see_paraview = false) {
    std::string outputFile = "outputFile.txt"; //  magical
    const int ghost_layer = 1;

    array<double, DIM> l_min, l_max;
    if (l_min_.size() != DIM or l_min_.size() != DIM) {
        throw runtime_error("Incoherent dimensions!");
    }
    for (size_t d = 0; d < DIM; d++) {
        l_min[d] = l_min_[d];
        l_max[d] = l_max_[d];
    }

    rsa_domain<DIM> domain(l_min, l_max, ghost_layer, rad);

    domain.domain_log();

    auto start = std::chrono::system_clock::now();
    //
    auto algo = algorithm::uniform_generate<DIM, method>(domain, rad, size, n_ndraw, seed);
    //
    auto end_algo_time = std::chrono::system_clock::now();
    std::chrono::duration<double> seconds_for_algo = end_algo_time - start;
    if (check_vox_time) {
        algo.check_vox_time();
    }
    std::chrono::duration<double> seconds_for_vox = std::chrono::system_clock::now() - end_algo_time;

    // total nb spheres
    int64_t nb_spheres = algo.get_grid().get_number_of_spheres_fast();
    int64_t total_nb_spheres = rsa_mpi::compute_mpi_sum(nb_spheres);
    //
    if (rsa_mpi::is_master_rank()) {
        write_to_file(outputFile, rsa_mpi::get_number_of_mpi_processes());
        write_to_file(outputFile, " ");
        write_to_file(outputFile, DIM);
        write_to_file(outputFile, " ");
        write_to_file(outputFile, rad);
        write_to_file(outputFile, " ");
        for (size_t i = 0; i < DIM; i++) {
            write_to_file(outputFile, l_min[i]);
            write_to_file(outputFile, " ");
        }
        for (size_t i = 0; i < DIM; i++) {
            write_to_file(outputFile, l_max[i]);
            write_to_file(outputFile, " ");
        }
        write_to_file(outputFile, seed);
        write_to_file(outputFile, " ");
        write_to_file(outputFile, total_nb_spheres);
        //
        write_to_file(outputFile, " ");
        write_to_file(outputFile, to_string(seconds_for_algo.count()));
        write_to_file(outputFile, " ");
        write_to_file(outputFile, to_string(seconds_for_vox.count()));
        write_to_file(outputFile, "\n");
    }

    if (see_paraview == 1) {
        rsa_paraview::paraview("ParaviewOutput", domain.get_grid());
    }
}

template<int DIM, class ...Args>
void test_perf_(int dim, Args... args) {
    if constexpr (DIM > 9) {
        throw std::runtime_error("Impossible!");
    } else {
        if (dim == DIM) {
            test_perf<DIM, 1>(args...);
        } else {
            test_perf_<DIM + 1>(dim, args...);
        }
    }
}

inline void test_perf(rsa_mpi::rsa_parameters& p) {
    test_perf_<1>(p.DIM, p.radius, p.l_min, p.l_max, p.seed, p.size, p.n_draw, false, static_cast<bool>(p.paraview));
}

template<class ...Args>
void test_perf(int dim, Args... args) {
    test_perf_<1>(dim, args...);
}
