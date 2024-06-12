//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <rsa_sphere.hxx>
#include <rsa_grid.hxx>
#include <rsa_random.hxx>
#include <rsa_domain.hxx>
#include <rsa_decoration.hxx>
//#define VERBOSE_FILL 1

int main(int argc, char** argv) {
	rsa_mpi::init(&argc, &argv);

	std::cout << "================= ===== ================= " << std::endl;
	std::cout << "==== rsa_algo with mpi parallelization == " << std::endl;
	std::cout << "================= start ================= " << std::endl;

	constexpr int n = 10;
	const int random_try = 100;
	constexpr int dim = 3;
	std::array<double, dim> domain_inf = { 0.0, 0.0, 0.0 };
	std::array<double, dim> domain_sup = { 1.0, 1.0, 1.0 };
	std::array<int, dim> sizes = { n, n, n };
	std::array<double, dim>  dx = { 1.0 / (n + 2), 1.0 / (n + 2), 1.0 / (n + 2) };
	const double rad = 0.01;
	rsa_domain<dim> domain(domain_inf, domain_sup);
	rsa_grid* grid = domain.get_grid();

	std::cout << "= option -- domain size: [" << domain_sup[0] << "," << domain_sup[1] << "," << domain_sup[2] << "]" << std::endl;
	std::cout << "= option -- discretization : [" << sizes[0] << "," << sizes[0] << "," << sizes[0] << "]" << std::endl;
	std::cout << "= option -- regular_try: " << n * n * n << std::endl;
	std::cout << "= option -- random_try: " << random_try << std::endl;
	std::cout << "= option -- default radius: " << rad << std::endl;


	domain.uniform_generate(rad, n * n * n);

	grid->print_log();
	//grid.print_spheres();
	std::cout << "====> write output file <=====" << std::endl;
	grid->write_spheres();


	rsa_mpi::finalize();
	std::cout << "=========== end ===========" << std::endl;
	return EXIT_SUCCESS;
}
