//! Copyright : Apache 2.0, see LICENSE 
//! 
#include<test.hxx>
//#define VERBOSE_FILL 1

int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);

    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase = { tuple<double, double, int> {0.25, 0.1, 0}, tuple<double, double, int> {0.125, 0.1, 1} };
    test_ND_radius<3, 1>(desired_radius_volumeFraction_phase);

    rsa_mpi::finalize();
    return EXIT_SUCCESS;
}
