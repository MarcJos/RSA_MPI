//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <basic_types.hxx>
#include <rsa_data_storage.hxx>
#include <rsa_random.hxx>
#include <rsa_domain.hxx>

//#define VERBOSE_FILL 1

int main(int argc, char** argv) {
	std::cout << "================= ===== ================= " << std::endl;
	std::cout << "==== rsa_algo with mpi parallelization == " << std::endl;
	std::cout << "================= start ================= " << std::endl;
	rsa_mpi::init(&argc, &argv);

	constexpr int n = 10;
	const int random_try = 100;
	constexpr int DIM = 3;
	std::array<double, DIM> domain = { 1.0, 1.0, 1.0 };
	std::array<int, DIM> sizes = { n, n, n };
	std::array<double, DIM>  dx = { 1.0 / (n + 2), 1.0 / (n + 2), 1.0 / (n + 2) };
	const double rad = 0.01;
	rsa_data_storage<DIM> data_storage;

	std::cout << "= option -- domain size: [" << domain[0] << "," << domain[1] << "," << domain[2] << "]" << std::endl;
	std::cout << "= option -- discretization : [" << sizes[0] << "," << sizes[0] << "," << sizes[0] << "]" << std::endl;
	std::cout << "= option -- regular_try: " << n * n * n << std::endl;
	std::cout << "= option -- random_try: " << random_try << std::endl;
	std::cout << "= option -- default radius: " << rad << std::endl;

	auto fill = [rad](rsa_data_storage<DIM>& a_data_storage, const double a_rx, const double a_ry, const double a_rz, const uint64_t a_id)->bool {
		rsa_sphere<DIM> sphere(std::array<double, DIM>{a_rx, a_ry, a_rz}, rad, a_id);
#ifdef VERBOSE_FILL
		std::cout << "= -> trying to add :" << " r " << rad << " rx " << a_rx << " ry " << a_ry << " rz " << a_rz << " id " << a_id << std::endl;
#endif
		if (a_data_storage.check_intersection(sphere)) {
			a_data_storage.add_sphere(a_id, sphere);
#ifdef VERBOSE_FILL
			std::cout << "=  --> added" << std::endl;
#endif
			return true;
		} else {
#ifdef VERBOSE_FILL
			std::cout << "=  --> not added" << std::endl;
#endif
			return false;
		}
		};


	std::cout << "====> predefined add <====" << std::endl;
	for (int z = 0; z < sizes[2]; z++) {
		for (int y = 0; y < sizes[1]; y++) {
			for (int x = 0; x < sizes[0]; x++) {
				const double rx = (x + 1) * dx[0];
				const double ry = (y + 1) * dx[1];
				const double rz = (z + 1) * dx[2];
				const uint64_t id = x + sizes[0] * (y + z * sizes[1]);
				fill(data_storage, rx, ry, rz, id);
			}
		}
	}

	const int added_spheres = data_storage.size();
	std::cout << "=----> numbero of spheres added: " << added_spheres << "/" << sizes[0] * sizes[1] * sizes[2] << std::endl;

	std::cout << "====> randomized add <====" << std::endl;
	// random add
	std::random_device rd;
	std::vector<law::uniform<double>> ulaw{};

	for (int i = 0; i < DIM; i++) {
		std::mt19937 tmp_gen(rd());
		std::uniform_real_distribution<> tmp_dis(0, domain[i]);
		ulaw.emplace_back(law::uniform<double>(tmp_gen, tmp_dis));
	}

	uint64_t current_number_of_spheres = data_storage.size();
	for (int it = 0; it < random_try; it++) {
		const double rx = ulaw[0]();
		const double ry = ulaw[0]();
		const double rz = ulaw[0]();
		auto ok = fill(data_storage, rx, ry, rz, current_number_of_spheres);
		if (ok) current_number_of_spheres++;
	}

	std::cout << "=----> numbero of spheres added: " << (current_number_of_spheres - added_spheres) << "/" << random_try << std::endl;

	data_storage.print_log();
	//data_storage.print_spheres();
	std::cout << "====> write output file <=====" << std::endl;
	data_storage.write_spheres();


	rsa_mpi::finalize();
	std::cout << "=========== end ===========" << std::endl;
	return EXIT_SUCCESS;
}
