//! Copyright : see license.txt
//!
//! \brief

#include<user_interface.hxx>
using namespace sac_de_billes;

int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);
    {
        constexpr int DIM = 3;
        Point<DIM> l = { 0., 0., 0. };
        Point<DIM> L = { 3.5, 3.5, 3.5 };
        vector<tuple<double, double, int>> desired_radius_volumeFraction_phase{};
        //desired_radius_volumeFraction_phase.push_back({ 0.05, 0.25, 0 });
        desired_radius_volumeFraction_phase.push_back({ 0.1, 1, 0 });
        double rad_max = 0.1;//0.05;
        int seed = 0;
        double distance_min = 0.;
        int nb_shots_single_draw = 6000;

        user_interface::rsa_domain<DIM> domain(l, L, rad_max);
        user_interface::RadiusGenerator<DIM> radius_generator(desired_radius_volumeFraction_phase, domain.get_total_volume(), distance_min);
        user_interface::rsa_algo<DIM> rsaalgo(domain, radius_generator, nb_shots_single_draw);
        rsaalgo.proceed(seed);

        rsa_paraview::paraview(domain.get_pointed_to());
    }
    rsa_mpi::finalize();

    return EXIT_SUCCESS;
}
