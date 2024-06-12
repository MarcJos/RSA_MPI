//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <rsa_data_storage.hxx>
#include <rsa_random.hxx>
#include <rsa_domain.hxx>
#include <rsa_decoration.hxx>
#include <operator_algorithm.hxx>
#include <operator_paraview.hxx>
#include <radius_generator.hxx>

template<int DIM, int method>
vector<rsa_sphere<DIM>> throw_spheres(sac_de_billes::Point<DIM> L,
    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
    size_t seed, bool write_paraview = true) {
    std::array<double, DIM> domain_inf = sac_de_billes::create_array<DIM>(0.);
    std::array<double, DIM> domain_sup = domain_inf + L;
    const int ghost_layer = 1;
    double volume = sac_de_billes::auxi_function::productOf<double>(L);
    sac_de_billes::RadiusGenerator<DIM> radius_generator(desired_radius_volumeFraction_phase, 1);
    rsa_domain<DIM> domain(domain_inf, domain_sup, ghost_layer, radius_generator.get_max_radius());

    domain.domain_log();
    algorithm::uniform_generate<DIM, method>(domain, radius_generator, 6000, 10, seed);

    if (write_paraview) {
        rsa_mpi::message("====> write mpi output file <=====");
        rsa_paraview::paraview(domain);
    }
    return domain.extract_spheres();
}

template<unsigned short DIM, int method>
void test_ND_radius(vector<tuple<double, double, int>> desired_radius_volumeFraction_phase) {
    std::array<double, DIM> domain_inf = create_array<DIM>(0.);
    std::array<double, DIM> domain_sup = create_array<DIM>(1.);
    const int ghost_layer = 1;
    double volume = 1.;
    sac_de_billes::RadiusGenerator<DIM> radius_generator(desired_radius_volumeFraction_phase, 1);
    rsa_domain<DIM> domain(domain_inf, domain_sup, ghost_layer, radius_generator.get_max_radius());

    domain.domain_log();
    size_t seed = 0;
    algorithm::uniform_generate<DIM, method>(domain, radius_generator, 6000, 10, seed);

    if constexpr (DIM == 2 or DIM == 3) {
        rsa_mpi::message("====> write mpi output file <=====");
        rsa_paraview::paraview(domain);
    }
}

template<unsigned short DIM, int method>
void test_ND(double rad) {
    std::array<double, DIM> domain_inf = create_array<DIM>(0.);
    std::array<double, DIM> domain_sup = create_array<DIM>(1.);
    const int ghost_layer = 1;
    rsa_domain<DIM> domain(domain_inf, domain_sup, ghost_layer, rad);

    domain.domain_log();
    size_t seed = 0;
    algorithm::uniform_generate<DIM, method>(domain, rad, 6000, 10, seed);
    if constexpr (DIM == 2 or DIM == 3) {
        rsa_mpi::message("====> write mpi output file <=====");
        rsa_paraview::paraview(domain);
    }
}

void test_radius_generator() {
    std::vector<std::tuple<double, uint64_t, int>> desired_radius_nb_phase =
    { std::make_tuple<double, uint64_t, int>(0.1, 10, 1) };
    sac_de_billes::RadiusGenerator<3> RadiusGenerator(desired_radius_nb_phase, 0);
}

//! @brief provoke subdivision and cause a bug for int might be too small to store voxel coordinates
void test_bug_2D() {
    Point<2> L = { 1, 1 };
    int seed = 0;
    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase = { std::make_tuple<double, double, int>(0.0005, 1, 0) };
    throw_spheres<2, 1>(L, desired_radius_volumeFraction_phase, seed);
}
