//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief

#include <user_interface.hxx>

int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);
    // NECESSARY for MPI

    omp_set_num_threads(1);

    // data
    constexpr int DIM = 3;
    Point<DIM> l_inf = { 0., 0., 0. };
    Point<DIM> l_sup = { 1., 1., 1. };
    // domain is [l_inf[0], l_sup[0]] x ... x [l_inf[DIM-1], l_sup[DIM-1]]

    vector<tuple<double, double, int>>
        desired_radius_volumeFraction_phase = { {0.1, 0.25, 0}, {0.05, 1., 0} };
       // desired_radius_volumeFraction_phase = { {0.05, 0.25, 0}, {0.025, 1., 0} };
    // volume fraction {{radius_0, volume fraction_0, phase_id_0}, ...}

    double rad_max = 0.;
    for (const auto rvp : desired_radius_volumeFraction_phase) {
        rad_max = std::max(rad_max, std::get<0>(rvp));
    }
    // r_max should be larger than the maximal radius in the simulation
    
    int seed = 0;
    // random seed
    
    double distance_min = 0.;
    // minimal distance between 2 different spheres

    int nb_shots_single_draw = 6000; 
    // number of shots before MPI communications

    ///////////////////////
    // launch the algorithm
    ///////////////////////

    user_interface::rsa_domain<DIM> domain{ l_inf, l_sup, rad_max };
    // domain
    user_interface::RadiusGenerator<DIM> radius_generator(
        desired_radius_volumeFraction_phase, domain.get_total_volume(), distance_min);
    // for generating radii along the simulation
    user_interface::rsa_algo<DIM> rsaalgo(domain, radius_generator, nb_shots_single_draw);
    // object for algorithm
    rsaalgo.proceed(seed);
    // launch the algorithm


    auto the_spheres = domain.extract_spheres();
    // get the output
    domain.paraview();
    // print the output

    rsa_mpi::message(std::to_string(the_spheres.size()));
    // for using the_spheres (only for compilation)


    rsa_mpi::finalize();
    // NECESSARY for MPI
    return EXIT_SUCCESS;
}
