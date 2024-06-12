//! Copyright : Apache 2.0, see LICENSE 

#pragma once

#include<rsa_domain.hxx>
#include<radius_generator.hxx>
#include<operator_algorithm.hxx>
#include<test.hxx>


namespace user_interface {

template<int DIM>
class rsa_domain {
    //! @brief: interface class for ::rsa_domain<DIM>. Acts as a reference.
    //! Pave a cuboid $`[l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]]`$."
    //! @tparam DIM : Dimension of space
public:
    //! @brief constructor for a periodic domain [l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]] in \R^D (D=DIM)
    //! The domain is scattered over the mpi processes.
    //! @param global_inf : vector  {lmin[0], ..., l_min[D-1]}
    //! @param global_sup : vector {lmax[0], ..., l_max[D-1]}
    //! @param a_rad : the size of each cell in any direction should be larger than a_rad. This radius is usually the largest radius of placed spheres.
    rsa_domain(const vec_d< DIM >& global_inf, const vec_d< DIM >& global_sup, const double a_rad)
        : rsa_domain_ptr_{ new ::rsa_domain<DIM>(global_inf, global_sup, 1, a_rad) } {}


    //! @brief Accesor used to get the inf boundary
    const vec_d<DIM>& get_inf() const { return rsa_domain_ptr_->get_inf(); }
    //! @brief Accesor used to get the sup boundary
    const vec_d<DIM>& get_sup() const { return rsa_domain_ptr_->get_sup(); }
    //! @brief: get the total volume of the domain
    double get_total_volume() const { return  rsa_domain_ptr_->get_total_volume(); }

    //! @return: the total volume of spheres contained in all the mpi domains
    double compute_total_volume_of_spheres() const { return rsa_domain_ptr_->compute_total_volume_of_spheres(); }
    //! @brief : Print some general informations about the domain.
    void domain_log() const { rsa_domain_ptr_->domain_log(); }
    //! @brief : Create a VTK file with all particles over subdomains.
    void paraview() const { rsa_paraview::paraview(*rsa_domain_ptr_); }
    //! @return list of spheres inside the domain
    //! @warning only on a single MPI processe
    std::vector<rsa_sphere<DIM>> extract_spheres() const { return rsa_domain_ptr_->extract_spheres(); }

    //! @return a reference to the underlying class (not interface)
    //! @warning: this is for a finer use.
    ::rsa_domain<DIM>& get_pointed_to() { return *rsa_domain_ptr_; }


private:
    std::shared_ptr<::rsa_domain<DIM>> rsa_domain_ptr_;
};

template<int DIM>
class RadiusGenerator {
    //! @brief class used for generating successively radii of spheres to be placed
    //! Interface class for sac_de_billes::RadiusGenerator<DIM>. Acts as a reference.
    //! @tparam DIM : Dimension of space
public:
    //! @brief: constructor from a list of desired radius and volume fraction
    //! \param desired_radius_nb_phase : vectors of desired (radius of spheres, nb_spheres, phase of spheres)
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(vector<tuple<double, uint64_t, int>> desired_radius_nb_phase,
        double exclusion_distance)
        :radius_generator_prt_{ new sac_de_billes::RadiusGenerator<DIM>(desired_radius_nb_phase, exclusion_distance) } {}
    //! @brief: constructor from a list of desired R and phi
    //! \param desired_radius_volumeFraction_phase : vectors of desired (radius of spheres, volume_fraction, phase of spheres)
    //! @param volume : total volume of the shape under consideration
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
        double volume,
        double exclusion_distance)
        :radius_generator_prt_{ new sac_de_billes::RadiusGenerator<DIM>(desired_radius_volumeFraction_phase, volume, exclusion_distance) } {}

    //! @return a reference to the underlying class (not interface)
    //! @warning: this is for a finer use.
    sac_de_billes::RadiusGenerator<DIM>& get_pointed_to() { return *radius_generator_prt_; }

private:
    std::shared_ptr<sac_de_billes::RadiusGenerator<DIM>> radius_generator_prt_;
};


template<int DIM>
class rsa_algo {
    //! @brief Class for execution of the algorithm.
    //! Interface class for algorithm::rsa_algo<DIM>. Acts as a reference.
    //! @tparam DIM : Dimension of space
public:
    //! @brief class for proceeding with the rsa_algo
    //! @param a_domain : contains all the cells and spheres
    //! @param a_radius_generator : a functor for producing new radii
    //! @param a_nbshots_singledraw : nb of local draws at each mpi proc, before exchanging with other procs
    rsa_algo(rsa_domain<DIM>& a_domain, user_interface::RadiusGenerator<DIM>& a_radius_generator,
        int a_nbshots_singledraw)
        :rsa_algo_ptr_{ new algorithm::rsa_algo<DIM>(a_domain.get_pointed_to(), a_radius_generator.get_pointed_to(), a_nbshots_singledraw) } {}

    //! @brief Execute the RSA algorithm.
    //! @param seed :for the pseudo-random generator.
    void proceed(size_t seed) { rsa_algo_ptr_->template proceed<1>(seed); }

    //! @brief: Reset the radius generator, but retains the already placed spheres.\n
    //! Proceed can then be appealed once more.
    void reset_radius_generator(RadiusGenerator<DIM>& a_radius_generator) { rsa_algo_ptr_->reset_radius_generator(a_radius_generator.get_pointed_to()); }


private:
    std::shared_ptr<algorithm::rsa_algo<DIM>> rsa_algo_ptr_;
};

//! @brief Perform the RSA algorithm
//! @return The local list of spheres inside the domain (NOT gathered over the MPI processes)
//! @tparam DIM : Dimension of space
//! @param L : size of the box [0, L]^3
//! @param desired_radius_volumeFraction_phase : vector of triple { (double) radius , (double) volumeFraction, (int) phase}
//! @param seed : seed for the random generator
//! @param write_paraview : write a paraview output
template<int DIM>
vector<rsa_sphere<DIM>> throw_spheres(sac_de_billes::Point<DIM> L,
    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
    size_t seed, bool write_paraview) {
    return ::throw_spheres<DIM, 1>(L, desired_radius_volumeFraction_phase, seed, write_paraview);
}

} // namespace  user_interface
