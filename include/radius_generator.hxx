//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <basic_types.hxx>
#include <rsa_random.hxx>

namespace sac_de_billes {

using namespace std;

template<int DIM>
class RadiusGenerator {
    //! @brief class used for generating successively radii of spheres to be placed
    //! @warning: the internal representation of radii takes into account a parameter called exclusionDistance, \see RadiusGenerator::setTabRadii
public:
    //! @brief: constructor from a list of desired radius and volume fraction
    //! \param desired_radius_nb_phase : vectors of desired (radius of spheres, nb_spheres, phase of spheres)
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(vector<tuple<double, uint64_t, int>> desired_radius_nb_phase,
        double exclusion_distance = 0);
    //! @brief: constructor from a list of desired R and phi
    //! \param desired_radius_volumeFraction_phase : vectors of desired (radius of spheres, volume_fraction, phase of spheres)
    //! @param volume : total volume of the shape under consideration
    //! \param exclusionDistance : minimal distance between spheres
    RadiusGenerator(vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
        double volume,
        double exclusion_distance = 0);

    //! @return are there still radii to be generated?
    bool is_there_still_radii() const { return uint64_t(current_index) < desired_radius_nb_phase.size(); }
    //! @return the current desired number of spheres
    uint64_t get_current_number() const { return get<1>(desired_radius_nb_phase[current_index]); }
    //! \return the maximal radius
    double get_max_radius() const { return max_radius; }
    //! \return the minimal radius
    double get_min_radius() const { return min_radius; }
    //! @param a_size : number of radii generated
    //! @return : a pair< vector_of_phases, vector_of_radii>, with vector_of_phases[i] corresponding to vector_of_radii[i]
    std::tuple<vec_int, vec_double> operator()(size_t a_size);
    //! @brief update the radius generator by the knowledge of placed spheres
    //! @param nb_placed_spheres 
    void update_placed(size_t nb_placed_spheres);

private:
    vector<tuple<double, uint64_t, int>> desired_radius_nb_phase;
    double exclusion_distance;
    double max_radius;
    double min_radius;
    int current_index;

private:
    //! \return the current radius
    double get_current_radius() const { return get<0>(desired_radius_nb_phase[current_index]); }
    //! \return the current phase
    int get_current_phase() const { return get<2>(desired_radius_nb_phase[current_index]); }
    //! @brief set the current desired number of spheres
    void set_current_number(uint64_t current_number) { get<1>(desired_radius_nb_phase[current_index]) = current_number; }
    //! @brief goes to the next class of radii
    void go_to_next_radius() { current_index++; }
};


template<class FUNCTION_TYPE>
class RandomRadiusGenerator {
public:
    //! @brief : internal constructor
    //! @param transform_ : nonlinear transform so that the resulting radius generated is
    //! transform_(x) for x a uniform variable in [0, 1]
    //! @param r_min_ : minimal radius that can be generated
    //! @param r_max_ : maximal radius that can be generated
    //! @param a_gen : a random engine
    RandomRadiusGenerator(FUNCTION_TYPE transform_, double r_min_, double r_max_) :
        r_min{ r_min_ }, r_max{ r_max_ }, transform{ transform_ }, uniform_law(nullptr) {}
    //! @brief : create the law with the given engine
    //! @param a_gen : a random engine
    //! @warning : this randomGenerator cannot be used without this
    void set_random_generator(std::mt19937& a_gen) {
        uniform_law.reset(new law::uniform<double>(0., 1., a_gen));
    }
    //! @return a radius randomly generated
    double operator()() {
        return transform((*uniform_law)());
    }

private:
    double r_min;
    double r_max;
    FUNCTION_TYPE transform;
    unique_ptr<law::uniform<double>> uniform_law;
};

namespace radius_generator_auxi {
template<int DIM>
vector<tuple<double, uint64_t, int>> compute_desired_radius_nb_phase(
    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
    double volume);
} // namespace  radius_generator_auxi

} // namespace  sac_de_billes

#include <radius_generator.ixx>
