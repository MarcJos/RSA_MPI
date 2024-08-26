//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include <list_of_voxels.hxx>
#include <rsa_grid_traversal.hxx>

namespace voxel_list {
using namespace sac_de_billes;
using namespace std;

template<int DIM>
class list_of_cell_voxels : public rsa_grid_traversal<DIM> {
private:
    std::vector<voxel_list::list_of_voxels<DIM>> m_list_voxel; ///< list of voxels per cell
    uint64_t number_of_voxel;

public:
    //! @brief : default constructor
    list_of_cell_voxels() : rsa_grid_traversal<DIM>(), m_list_voxel{}, number_of_voxel{ 0 }{}
    //! @brief : constructor
    list_of_cell_voxels(const double a_rad, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup, double a_max_diagonal)
        : rsa_grid_traversal<DIM>(a_rad, 0, a_inf, a_sup), m_list_voxel{}, number_of_voxel{ 0 } {
        m_list_voxel.reserve(this->size());
        //
        Point<DIM> origin = create_array<DIM>(0.);
        for (size_t i = 0; i < this->size(); i++) {
            auto disc_coord = this->get_coordinates(i);
            for (size_t d = 0; d < DIM; d++) {
                origin[d] = this->get_shift()[d] + disc_coord[d] * this->get_cell_length()[d];
            }
            m_list_voxel.push_back(list_of_voxels<DIM>(origin, this->get_cell_length(), a_max_diagonal));
        }
    };


    //! @brief remove all voxels totally covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void remove_covered(const RSA_GRID& a_rsa_grid, double a_minimal_radius) {
        for (size_t i = 0; i < this->size(); i++) {
            m_list_voxel[i].remove_covered(a_rsa_grid, a_minimal_radius);
        }
    }
    //! @brief subdivide each voxel, but retain only those that are not covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void subdivide_uncovered(const RSA_GRID& rsa_cell, double a_minimal_radius) {
        for (size_t i = 0; i < this->size(); i++) {
            m_list_voxel[i].subdivide_uncovered(rsa_cell, a_minimal_radius);
        }
    }
    //! @brief get the total number of voxels
    size_t nb_voxels() const {
        size_t nb_voxels_tot = 0;
        for (size_t i = 0; i < this->size(); i++) {
            nb_voxels_tot += m_list_voxel[i].size();
        }
        return nb_voxels_tot;
    }
    //! @return the area of all the voxels
    double total_area() const {
        double res = 0;
        for (size_t i = 0; i < this->size(); i++) {
            res += m_list_voxel[i].total_area();
        }
        return res;
    }

    // random picking

    //! @return : a_nb_pts randomly picked in the voxels (with uniform law)
    arr_vec_double<DIM> pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const {
        throw runtime_error("Not programmed yet!");
    }

};

} // namespace voxel_list