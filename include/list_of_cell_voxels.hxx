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
    std::vector<voxel_list::list_of_voxels<DIM>> m_voxel; ///< list of voxels per cell
    uint64_t number_of_voxel;

public:
    //! @brief : default constructor
    list_of_cell_voxels() : rsa_grid_traversal<DIM>(), m_voxel{}, number_of_voxel{ 0 }{}
    //! @brief : this constructor calls the function "define".
    //! @see : define.
    list_of_cell_voxels(const double a_rad, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup, double a_max_diagonal)
        : rsa_grid_traversal<DIM>(a_rad, 0, a_inf, a_sup), m_voxel{}, number_of_voxel{ 0 } {
        m_voxel.reserve(this->size());
        //
        Point<DIM> origin = create_array<DIM>(0.);
        for (size_t i = 0; i < this->size(); i++) {
            auto disc_coord = this->get_coordinates(i);
            for (size_t d = 0; d < DIM; d++) {
                origin[d] = this->get_shift()[d] + disc_coord[d] * this->get_cell_length()[d];
            }
            m_voxel.push_back(list_of_voxels<DIM>(origin, this->get_cell_length(), a_max_diagonal));
        }
    };
};

} // namespace voxel_list