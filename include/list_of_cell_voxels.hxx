//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include <random>
#include <list_of_voxels.hxx>
#include <rsa_grid_traversal.hxx>
#include <link_to_auxi.hxx>

namespace voxel_list {
using namespace sac_de_billes;
using namespace std;

template<int DIM>
class list_of_cell_voxels : public rsa_grid_traversal<DIM> {
private:
    //! @brief internal forest of voxel trees
    //! each m_list_voxel[i] correspond to the cell i of the rsa_grid_traversal
    //! @warning : all voxels are supposed to be of same area
    std::vector<voxel_list::list_of_voxels<DIM>> m_list_voxel;

public:
    //! @brief : constructor
    list_of_cell_voxels(const rsa_grid_traversal<DIM>& ref_grid, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup, double a_max_diagonal);

    //! @brief remove all voxels totally covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void remove_covered(const RSA_GRID& a_rsa_grid, double a_minimal_radius);
    //! @brief subdivide each voxel, but retain only those that are not covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void subdivide_uncovered(const RSA_GRID& rsa_cell, double a_minimal_radius);

    //! @brief get the total number of voxels
    size_t nb_voxels() const { return nb_voxels_; }
    //! @return the area of all the voxels
    double total_area() const { return total_area_; }

    // random picking
    //! @return : a_nb_pts randomly picked in the voxels (with uniform law)
    arr_vec_double<DIM> pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const;

private:
    //! @brief stored total area (should be sometimes updated)
    double total_area_;
    //! @brief stored total nb of voxels (should be sometimes updated)
    size_t nb_voxels_;
    //! @brief : list of cell_id, voxel_id
    std::vector<tuple<int64_t, int64_t>> pointer_to_voxels;

    //! @brief update all relative components of the object
    void update();
    //! @brief update the pointer_to_voxels
    void update_pointer_to_voxels();
    //! @brief update total_area_ and nb_voxels_
    void update_area_and_nb_vox();
};

} // namespace voxel_list

#include <list_of_cell_voxels.ixx>