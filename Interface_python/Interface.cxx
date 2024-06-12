//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief 
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include<iostream>
#include<functional>
#include<omp.h>

#include "../include/test.hxx"
#include "../include/user_interface.hxx"

namespace py = pybind11;

void set_nb_threads(int nb_threads) {
    omp_set_num_threads(nb_threads);
#pragma omp parallel
    { // begin parallel section
#pragma omp master
        {
            std::cerr << ("Number of openmp threads = " + std::to_string(omp_get_num_threads()) + "\n");
        }
    } // end parallel section
}


void createModule_rsa_mpi_py(py::module_& rsa_mpi_py) {
    rsa_mpi_py.def("set_nb_threads", &set_nb_threads);
}



template<unsigned short DIM>
inline void create_functions_depending_on_dimension(py::module_& my_module, std::string DIM_S) {

    auto my_string_rsa_sphere = "sphere" + DIM_S;
    py::class_<rsa_sphere<DIM>>(my_module, my_string_rsa_sphere.c_str(),
        "Representation of sphere.")
        .def(py::init<>())
        .def(py::init<const Point<DIM> &, double, int>())
        .def_readwrite("center", &Sphere<DIM>::center)
        .def_readwrite("radius", &Sphere<DIM>::radius)
        .def_readwrite("phase", &Sphere<DIM>::phase)
        ;

    auto my_string_thowSpheres = "throw_spheres" + DIM_S;
    my_module.def(my_string_thowSpheres.c_str(), &user_interface::throw_spheres<DIM>,
        "([double, ...] L, [[double, double, int]] desired_radius_volumeFraction_phase, seed, write_paraview)\n"
        "L : size of the box [0, L]^3\n"
        "desired_radius_volumeFraction_phase : vector of triples { (double) radius , (double) volumeFraction, (int) phase}\n"
        "seed : indicates the random seed for the pseudo-random generator. \n"
        "write_paraview : write a paraview output \n"
        "return The local list of spheres inside the domain (NOT gathered over the MPI processes)");

    auto my_string_rsa_domain = "rsa_domain" + DIM_S;
    py::class_<user_interface::rsa_domain<DIM>>(my_module, my_string_rsa_domain.c_str(),
        "Data storage for the spheres.\n"
        "Pave a cuboid $`[l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]]`$.\n")
        .def(py::init<const vec_d<DIM>&, const vec_d<DIM>&, const double>(),
            "Constructor\n."
            "parameters : ([double, ...] global_inf, [double, ...] global_sup, maximal_radius)")
        .def("domain_log", &user_interface::rsa_domain<DIM>::domain_log)
        .def("get_total_volume", &user_interface::rsa_domain<DIM>::get_total_volume,
            "Return the total volume of the (global) cuboid.")
        .def("paraview", &user_interface::rsa_domain<DIM>::paraview,
            "Print paraview output.")
        .def("compute_total_volume_of_spheres", &user_interface::rsa_domain<DIM>::compute_total_volume_of_spheres,
            "Compute the cumulated volume of all spheres inside the whole domains")
        .def("extract_spheres", &user_interface::rsa_domain<DIM>::extract_spheres,
            "return list of spheres inside the domain"
            "warning: only on a single MPI processe")
        ;

    auto my_string_radius_generator = "radius_generator" + DIM_S;
    py::class_<user_interface::RadiusGenerator<DIM>>(my_module, my_string_radius_generator.c_str(),
        "Store the prescribed radius distribution fixed by the user.\n"
        "It evolves dynamically when algorithm is performed.")
        .def(py::init<vector<tuple<double, uint64_t, int>>, double>(),
            "`radius_generator_3D( [[rad[0], N[0], p[0]], ...], exclusion_distance)`:\n"
            "First place `N[0]` spheres of radius `rad[0]`, denoted by phase `p[0]`, then `N[1]` spheres...\n"
            "All spheres are separated by a distance `exclusion_distance`.")
        .def(py::init<vector<tuple<double, double, int>>, double, double>(),
            "`radius_generator_3D([[rad[0], phi[0], p[0]], ...], volume_total, exclusion_distance)`:\n"
            "First place `N` spheres of radius `rad[0]` denoted by phase `p[0]`  such that the resulting volume fraction is `phi[0]`, then `N` spheres...\n"
            "All spheres are separated by a distance `exclusion_distance`. \n"
            "The total volume of the domain is `volume_total`.")
        ;

    auto my_string_rsa_algo = "rsa_algo" + DIM_S;
    py::class_<user_interface::rsa_algo<DIM>>(my_module, my_string_rsa_algo.c_str(),
        "Class for execution of the algorithm.")
        .def(py::init<user_interface::rsa_domain<DIM>&, user_interface::RadiusGenerator<DIM>&, int>(),
            "Constructor")
        .def("proceed", &user_interface::rsa_algo<DIM>::proceed,
            "proceed(seed). Execute the RSA algorithm, with the seed for the pseudo-random generator.")
        .def("reset_radius_generator", &user_interface::rsa_algo<DIM>::reset_radius_generator,
            "reset_radius_generator(a_radius_generator)\n"
            "Reset the radius generator, but retains the already placed spheres.\n"
            "Proceed can then be appealed once more.")
        ;

}


PYBIND11_MODULE(rsa_mpi_py, rsa_mpi_py_) {
    createModule_rsa_mpi_py(rsa_mpi_py_);
    create_functions_depending_on_dimension<2>(rsa_mpi_py_, std::string("_2D"));
    create_functions_depending_on_dimension<3>(rsa_mpi_py_, std::string("_3D"));
    create_functions_depending_on_dimension<4>(rsa_mpi_py_, std::string("_4D"));
    create_functions_depending_on_dimension<5>(rsa_mpi_py_, std::string("_5D"));
    create_functions_depending_on_dimension<6>(rsa_mpi_py_, std::string("_6D"));
    create_functions_depending_on_dimension<7>(rsa_mpi_py_, std::string("_7D"));
}
