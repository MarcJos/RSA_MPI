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
    //! assumptions:
    //! all voxels are of the same size m_length
private:
    std::vector<voxel_list::list_of_voxels<DIM>> m_list_voxel; ///< list of voxels per cell
    uint64_t number_of_voxel;

public:
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
        rebuild();
    };


    //! @brief remove all voxels totally covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void remove_covered(const RSA_GRID& a_rsa_grid, double a_minimal_radius) {
        for (size_t i = 0; i < this->size(); i++) {
            m_list_voxel[i].remove_covered(a_rsa_grid, a_minimal_radius);
        }
        rebuild();
    }
    //! @brief subdivide each voxel, but retain only those that are not covered by a sphere
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    void subdivide_uncovered(const RSA_GRID& rsa_cell, double a_minimal_radius) {
        for (size_t i = 0; i < this->size(); i++) {
            m_list_voxel[i].subdivide_uncovered(rsa_cell, a_minimal_radius);
        }
        rebuild();
    }
    //! @brief get the total number of voxels
    size_t nb_voxels() const { return nb_voxels_; }
    //! @return the area of all the voxels
    double total_area() const { return total_area_; }

    // random picking
    //! @return : a_nb_pts randomly picked in the voxels (with uniform law)
    arr_vec_double<DIM> pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const {
        if (this->nb_voxels() == 0 and a_nb_pts != 0) {
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
            throw std::runtime_error("Impossible to draw centers in no voxels");
        }
        // init ret
        arr_vec_double<DIM> ret(a_nb_pts);
        // uniform random distribution on (0, 1)
        law::uniform<double> coord_picker(0, 1, random_generator);
        // uniform random distribution for voxels
        law::uniform<int> voxel_picker(0, this->nb_voxels() - 1, random_generator);
        //
        for (int64_t i = 0; i < a_nb_pts; i++) {
            int64_t id_voxel = voxel_picker();
            Point<DIM> orig_voxel = auxi::orig_voxel(this->m_origin, *(pointer_to_voxels[id_voxel]), this->m_voxel_lengths);
            for (int d = 0; d < DIM; d++) {
                ret[i][d] = orig_voxel[d] + this->m_voxel_lengths[d] * coord_picker();
            }
        }
        return ret;
    }

private:
    double total_area_;
    size_t nb_voxels_;
    //! @brief pointer to the voxels
    std::vector<const VoxelCoordinates<DIM>*> pointer_to_voxels;

    void rebuild() {
        recompute_nb_voxels();
        recompute_total_area();
        pointer_to_voxels.reserve(nb_voxels());
        pointer_to_voxels.resize(0);
        for (const auto& list_vox : m_list_voxel) {
            const auto& vox_coord = list_vox.get_voxel_coordinates();
            for (const auto& single_vox_coord : vox_coord) {
                pointer_to_voxels.push_back(&single_vox_coord);
            }
        }
    }

    void recompute_nb_voxels() {
        nb_voxels_ = 0;
        for (size_t i = 0; i < this->size(); i++) {
            nb_voxels_ += m_list_voxel[i].size();
        }
    }

    void recompute_total_area() {
        total_area_ = 0;
        for (size_t i = 0; i < this->size(); i++) {
            total_area_ += m_list_voxel[i].total_area();
        }
    }
};

} // namespace voxel_list