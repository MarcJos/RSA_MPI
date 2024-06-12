//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief

#include <rsa_data_storage.hxx>
#include <rsa_random.hxx>
#include <rsa_domain.hxx>
#include <rsa_decoration.hxx>
#include <operator_algorithm.hxx>
#include <radius_generator.hxx>

template<int DIM, int method>
void test_max_frac(double rad, vector<double> l_min_, vector<double> l_max_, size_t seed,
    size_t size, size_t n_draw, double multiplier, int jmax, bool see_paraview = false) {
    // domain
    const int ghost_layer = 1;
    array<double, DIM> l_min, l_max;
    if (l_min_.size() != DIM or l_min_.size() != DIM) {
        throw runtime_error("Incoherent dimensions!");
    }
    for (size_t d = 0; d < DIM; d++) {
        l_min[d] = l_min_[d];
        l_max[d] = l_max_[d];
    }
    rsa_domain<DIM> domain(l_min, l_max, ghost_layer, rad);
    domain.domain_log();
    //
    sac_de_billes::RadiusGenerator<DIM> radius_generator(
        vector<tuple<double, double, int>>{ {rad, 1., 0} }, domain.get_total_volume()
    );

    auto print_here = [&]() {
        double total_volume = domain.compute_total_volume_of_spheres();
        rsa_mpi::message(" Total volume = " + to_string(total_volume));
        };

    algorithm::rsa_algo<DIM> rsaalgo(domain, radius_generator, size, n_draw);
    rsaalgo.template proceed<method>(seed);
    print_here();

    double current_rad = rad;
    for (int j = 0; j < jmax; j++) {
        auto start = std::chrono::system_clock::now();
        {
            current_rad *= multiplier;
            sac_de_billes::RadiusGenerator<DIM> radius_generator(
                vector<tuple<double, double, int>>{ {current_rad, 1., 0} }, domain.get_total_volume()
            );
            rsaalgo.reset_radius_generator(radius_generator);
            rsaalgo.template proceed<method>(seed);
            print_here();
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        rsa_mpi::message("Time for " + std::to_string(j) + ":  " + std::to_string(elapsed_seconds.count()) + "s");

        if (see_paraview == 1) {
            rsa_paraview::paraview(domain);
        }
    }
}



int main(int argc, char** argv) {
    rsa_mpi::init(&argc, &argv);

    int dim = abs(std::atof(argv[1]));
    double rad = std::atof(argv[2]);
    double L = std::atof(argv[3]);
    int seed = std::atof(argv[4]);
    int jmax = std::atof(argv[5]);
    double multiplier = std::atof(argv[6]);

    omp_set_num_threads(1);

    std::vector<double> l_inf, l_sup;
    l_inf.resize(dim);
    l_sup.resize(dim);
    for (int d = 0; d < dim; d++) {
        l_inf[d] = 0;
        l_sup[d] = L;
    }

    if (dim == 2) {
        test_max_frac<2, 1>(rad, l_inf, l_sup, seed,
            500, 10, multiplier, jmax, true);
    } else if (dim == 3) {
        test_max_frac<3, 1>(rad, l_inf, l_sup, seed,
            500, 10, multiplier, jmax, true);
    }

    rsa_mpi::finalize();
    return EXIT_SUCCESS;
}

