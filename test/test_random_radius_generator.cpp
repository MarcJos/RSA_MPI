//! Copyright : Apache 2.0, see LICENSE 
//! 
#include<test.hxx>

int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);
    constexpr int DIM = 3;
    std::array<double, DIM> domain_inf = create_array<DIM>(0.);
    std::array<double, DIM> domain_sup = create_array<DIM>(1.);
    const int ghost_layer = 1;


    // uniform law between (0.025 and 0.1)
    uint64_t phase = 0;
    double r_min = 0.025;
    double r_max = 0.1;
    auto nonlinear_transform = [r_min, r_max](double x) {return r_min + (r_max - r_min) * x;};
    uint64_t nb_spheres_max = 1e16; // full configuration
    sac_de_billes::RandomRadiusGenerator random_radius_generator(nonlinear_transform, phase);

    double exclusion_distance = 0.;

    sac_de_billes::RadiusGenerator<DIM> radius_generator(r_min, r_max, &random_radius_generator, nb_spheres_max, exclusion_distance);
    rsa_domain<DIM> domain(domain_inf, domain_sup, ghost_layer, radius_generator.get_max_radius());

    domain.domain_log();
    size_t seed = 0;
    algorithm::uniform_generate<DIM>(domain, radius_generator, 6000, 10, seed);

    if constexpr (DIM == 2 or DIM == 3) {
        rsa_mpi::message("====> write mpi output file <=====");
        rsa_paraview::paraview(domain);
    }

    rsa_mpi::finalize();
    return EXIT_SUCCESS;
}
