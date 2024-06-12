//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <stdlib.h>     /* atof */
#include<chrono>
#include <omp.h>

#include<test.hxx>
#include <write_xyz.hpp>
//#define VERBOSE_FILL 1

template<int method>
void test_plugin_write_xyz(double rad, double x, double y, double z) {
	omp_set_num_threads(1);
	constexpr int DIM = 3;
	std::array<double, DIM> domain_inf = { 0.0,0.0,0.0 };
	std::array<double, DIM> domain_sup = { x, y, z };
	const int ghost_layer = 1;
	rsa_domain<DIM> domain(domain_inf, domain_sup, ghost_layer, rad);
	size_t seed = 0;
	auto algo = algorithm::uniform_generate<DIM, method>(domain, rad, 6000, 10, seed);
	auto spheres = domain.get_grid().template extract_data<TypeCell::Real>();;
	rsa_mpi::message("====> write xyz output file <=====");
	write_xyz<DIM>("input.xyz", spheres, domain_sup);
}

int main(int argc, char** argv) {
	assert(argc == 5);
	// argv 2 = radius
	// argv 3 = sup_x
	// argv 4 = sup_y
	// argv 5 = sup_z
	rsa_mpi::init(&argc, &argv);
	int nb_threads = 1;
	test_plugin_write_xyz<1>(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]));
	rsa_mpi::finalize();
	return EXIT_SUCCESS;
}
