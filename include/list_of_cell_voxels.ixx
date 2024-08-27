//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace voxel_list {
using namespace sac_de_billes;
using namespace std;

template<int DIM>
list_of_cell_voxels<DIM>::list_of_cell_voxels(const rsa_grid_traversal<DIM>& ref_grid, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup, double a_max_diagonal)
    : rsa_grid_traversal<DIM>(ref_grid), m_list_voxel{},
    total_area_{ 0. }, nb_voxels_{ 0 }, pointer_to_voxels{}{
    m_list_voxel.reserve(this->size());
    //
    Point<DIM> origin = create_array<DIM>(0.);
    for (size_t i = 0; i < this->size(); i++) {
        auto disc_coord = this->get_coordinates(i);
        for (size_t d = 0; d < DIM; d++) {
            origin[d] = this->get_shift()[d] + disc_coord[d] * this->get_cell_length()[d];
        }
        if (not this->is_ghost(i)) { // legitimate cell
            m_list_voxel.push_back(list_of_voxels<DIM>(origin, this->get_cell_length(), a_max_diagonal));
        } else { // never pick spheres in ghost cells
            m_list_voxel.push_back(list_of_voxels<DIM>(origin));
        }
    }
    update();
}

template<int DIM>
template<class RSA_GRID>
void list_of_cell_voxels<DIM>::remove_covered(const RSA_GRID& a_rsa_grid, double a_minimal_radius) {
    for (size_t i = 0; i < this->size(); i++) {
        m_list_voxel[i].remove_covered(a_rsa_grid, a_minimal_radius);
    }
    update();
}

template<int DIM>
template<class RSA_GRID>
void list_of_cell_voxels<DIM>::subdivide_uncovered(const RSA_GRID& rsa_cell, double a_minimal_radius) {
    for (size_t i = 0; i < this->size(); i++) {
        m_list_voxel[i].subdivide_uncovered(rsa_cell, a_minimal_radius);
    }
    update();
}


template<int DIM>
arr_vec_double<DIM> list_of_cell_voxels<DIM>::pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const {
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
        const auto& voxel_list = this->m_list_voxel[get<0>(pointer_to_voxels[id_voxel])];
        Point<DIM> orig_voxel = voxel_list.orig_voxel(get<1>(pointer_to_voxels[id_voxel]));
        for (int d = 0; d < DIM; d++) {
            ret[i][d] = orig_voxel[d] + voxel_list.get_voxel_lengths()[d] * coord_picker();
        }
    }
    return ret;
}

template<int DIM>
void  list_of_cell_voxels<DIM>::update() {
    update_area_and_nb_vox();
    update_pointer_to_voxels();
}

template<int DIM>
void list_of_cell_voxels<DIM>::update_pointer_to_voxels() {
    pointer_to_voxels.resize(nb_voxels());
    int64_t k = 0;
    for (int64_t i = 0; i < m_list_voxel.size(); i++) {
        for (int64_t j = 0; j < m_list_voxel[i].nb_voxels(); j++) {
            pointer_to_voxels[k] = tuple<int64_t, int64_t>{ i, j };
            k++;
        }
    }
}

template<int DIM>
void list_of_cell_voxels<DIM>::update_area_and_nb_vox() {
    nb_voxels_ = 0;
    for (size_t i = 0; i < this->size(); i++) {
        nb_voxels_ += m_list_voxel[i].nb_voxels();
    }
    total_area_ = 0;
    for (size_t i = 0; i < this->size(); i++) {
        total_area_ += m_list_voxel[i].total_area();
    }
}

} // namespace voxel_list