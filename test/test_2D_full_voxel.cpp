//! Copyright : Apache 2.0, see LICENSE 
//! 
#include<test.hxx>
//#define VERBOSE_FILL 1

int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);

    test_ND<2, 1>(0.05);

    rsa_mpi::finalize();
    return EXIT_SUCCESS;
}
