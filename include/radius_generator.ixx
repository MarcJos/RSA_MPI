//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include<link_to_auxi.hxx>

namespace sac_de_billes {

template<int DIM, class RANDOM_RADIUS_GENERATOR>
RadiusGenerator<DIM, RANDOM_RADIUS_GENERATOR>::RadiusGenerator(vector<tuple<double, uint64_t, int>> desired_radius_nb_phase_,
    double exclusion_distance_) :
    desired_radius_nb_phase(desired_radius_nb_phase_), // temporary
    exclusion_distance(exclusion_distance_),
    max_radius{ 0 }, min_radius{ 0 }, // temporary
    current_index{ 0 },
    custom_radius_generator{ nullptr } {
    if (desired_radius_nb_phase.size() == 0) {
        throw runtime_error("Empty radius generator!");
    }
    std::sort(desired_radius_nb_phase.begin(), desired_radius_nb_phase.end(), [](auto& tup1, auto& tup2) {
        return get<0>(tup1) > get<0>(tup2);
        }); // descending order
    for (size_t i = 0; i < desired_radius_nb_phase.size(); i++) {
        get<0>(desired_radius_nb_phase[0]) += 0.5 * exclusion_distance;
    }
    max_radius = get<0>(desired_radius_nb_phase[0]);
    min_radius = get<0>(desired_radius_nb_phase[desired_radius_nb_phase.size() - 1]);
    exclusion_distance = exclusion_distance_;
}

template<int DIM, class RANDOM_RADIUS_GENERATOR>
RadiusGenerator<DIM, RANDOM_RADIUS_GENERATOR>::RadiusGenerator(vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
    double volume, double exclusion_distance_) :
    RadiusGenerator(
        radius_generator_auxi::compute_desired_radius_nb_phase<DIM>(
            desired_radius_volumeFraction_phase, volume),
        exclusion_distance_) {}

template<int DIM, class RANDOM_RADIUS_GENERATOR>
RadiusGenerator<DIM, RANDOM_RADIUS_GENERATOR>::RadiusGenerator(double min_radius_, double max_radius_,
    const RANDOM_RADIUS_GENERATOR* randomRadiusGenerator,
    size_t desired_nb_spheres, double exclusion_distance_) :
    desired_radius_nb_phase({ tuple<double, uint64_t, int>(0., desired_nb_spheres, 0) }),
    exclusion_distance{ exclusion_distance_ },
    max_radius{ max_radius_ },
    min_radius{ min_radius_ },
    current_index{ 0 },
    custom_radius_generator{ randomRadiusGenerator } {}

template<int DIM, class RANDOM_RADIUS_GENERATOR>
std::tuple<vec_int, vec_double> RadiusGenerator<DIM, RANDOM_RADIUS_GENERATOR>::operator()(size_t a_size, std::mt19937& random_generator) const {
    std::tuple<vec_int, vec_double> result;
    if (not custom_radius_generator) {
        result = std::tuple<vec_int, vec_double>(vec_int(a_size, this->get_current_phase()),
            vec_double(a_size, this->get_current_radius()));
    } else {
        result = custom_radius_generator->operator()(a_size, random_generator);
    }
    return result;
}

template<int DIM, class RANDOM_RADIUS_GENERATOR>
void RadiusGenerator<DIM, RANDOM_RADIUS_GENERATOR>::update_placed(size_t nb_placed_spheres) {
    this->set_current_number(this->get_current_number() - nb_placed_spheres);
    if (this->get_current_number() == 0) {
        this->go_to_next_radius();
    }
}


template<int DIM>
vector<tuple<double, uint64_t, int>> radius_generator_auxi::compute_desired_radius_nb_phase(
    vector<tuple<double, double, int>> desired_radius_volumeFraction_phase,
    double volume) {
    vector<tuple<double, uint64_t, int>> result(desired_radius_volumeFraction_phase.size());
    for (size_t i = 0; i < desired_radius_volumeFraction_phase.size(); i++) {
        auto current_tuple = desired_radius_volumeFraction_phase[i];
        double radius = get<0>(current_tuple);
        double volume_fraction = get<1>(current_tuple);
        int phase = get<2>(current_tuple);
        uint64_t nb_spheres = static_cast<uint64_t>(volume * volume_fraction / sac_de_billes::sphereTools::volumeSphere<DIM>(radius));
        result[i] = make_tuple<double, uint64_t, int>(std::move(radius), std::move(nb_spheres),
            std::move(phase));
    }
    return result;
}

} // namespace  sac_de_billes

