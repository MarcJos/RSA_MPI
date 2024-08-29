//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include "link_to_auxi.hxx"

namespace voxel_list {
using namespace sac_de_billes;
using namespace std;

template<int DIM>
list_of_voxels<DIM>::list_of_voxels(Point<DIM> a_global_origin, Point<DIM> a_origin, Point<DIM> a_domain_length, double a_max_diagonal) :
    list_of_voxels(a_global_origin, a_origin,
        a_domain_length, a_max_diagonal,
        list_of_voxels<DIM>::compute_voxel_length(a_domain_length, a_max_diagonal)) {}

template<int DIM>
list_of_voxels<DIM>::list_of_voxels(Point<DIM> a_global_origin, Point<DIM> a_origin, Point<DIM> a_domain_length, double a_max_diagonal, Point<DIM> a_cell_length) :
    list_of_voxels(a_global_origin, a_origin,
        list_of_voxels<DIM>::compute_shift_nb_vox(a_global_origin, a_origin, list_of_voxels<DIM>::compute_voxel_length(a_cell_length, a_max_diagonal)),
        list_of_voxels<DIM>::compute_nb_vox(a_domain_length, a_cell_length, a_max_diagonal),
        list_of_voxels<DIM>::compute_voxel_length(a_cell_length, a_max_diagonal)) {}


template<int DIM>
list_of_voxels<DIM>::list_of_voxels(Point<DIM> a_global_origin, Point<DIM> a_origin,
    vec_i<DIM> a_shift,
    vec_i<DIM> a_nb_vox, Point<DIM> a_voxel_length) :
    m_voxel_coordinates{},
    m_voxel_lengths(a_voxel_length),
    m_global_origin(a_global_origin),
    m_corners_voxel{} {
    array<array<int, 2>, DIM> limits;
    // warning : not very efficient!
    for (int d = 0; d < DIM; d++) {
        limits[d][0] = 0 + a_shift[d];
        limits[d][1] = a_nb_vox[d] + a_shift[d];
    }
    m_voxel_coordinates = sac_de_billes::getAllIndices<DIM, int64_t>(limits);
    auxi::create_corners_voxel_inplace<DIM>(m_voxel_lengths, m_corners_voxel);
}

template<int DIM>
template<class TEST>
void list_of_voxels<DIM>::remove_if(TEST test) {
    int64_t old_size = this->nb_voxels();
    BooleanVector is_voxel_covered(old_size, false);
    // tests
    for (int64_t id_vox = 0; id_vox < old_size; id_vox++) {
        is_voxel_covered[id_vox] = test(id_vox);
    }
    // remove
    for (int64_t id_vox = 0; id_vox < old_size; id_vox++) {
        if (is_voxel_covered[id_vox]) {
            std::memcpy(&(m_voxel_coordinates[id_vox]),
                &(m_voxel_coordinates[(old_size - 1)]),
                sizeof(m_voxel_coordinates[0]));
            is_voxel_covered[id_vox] = is_voxel_covered[(old_size - 1)];
            id_vox--;
            old_size--;
        }
    }
    // last coordinates are not used, because of covered_voxels
    m_voxel_coordinates.resize(old_size);
}

template<int DIM>
template<class RSA_GRID>
void list_of_voxels<DIM>::remove_covered(const RSA_GRID& a_rsa_grid, double a_minimal_radius) {
    auto test_covered = [&](size_t id_vox) {
        return auxi::is_covered<DIM>(this->orig_voxel(id_vox), this->m_voxel_lengths,
            this->get_corners_voxel(), a_rsa_grid, a_minimal_radius);
        };
    remove_if(test_covered);
}

template<int DIM>
template<class RSA_GRID>
void list_of_voxels<DIM>::subdivide_uncovered(const RSA_GRID& a_rsa_grid, double a_minimal_radius) {


    if (this->nb_voxels() == 0) { return; }
    auto new_voxel_lengths = 0.5 * m_voxel_lengths; // size of voxels divided by 2.
    //
    uint64_t old_size = this->nb_voxels();
    constexpr int nb_corners = auxi_function::puissance<DIM>(2);

    const auto& tabcorner = sac_de_billes::path::TabCorner<DIM>::get().getTab();
    auto new_corners_voxel = m_corners_voxel;
    auxi::create_corners_voxel_inplace<DIM>(new_voxel_lengths, new_corners_voxel);

    static std::vector<DiscPoint<DIM>> new_voxel_coordinates{}; // static for allocating memory not too often
    new_voxel_coordinates.resize(10);
    // Will store things for identifying spheres
    vector<pair<int, int>> list_of_cell_id_sphere_id{};
    uint64_t k = 0;
    for (uint64_t i = 0; i < old_size; i++) {
        Point<DIM> center_voxel = this->center(i);
        auxi::get_list_of_spheres<DIM>(center_voxel, m_voxel_lengths, a_rsa_grid, a_minimal_radius, list_of_cell_id_sphere_id);
        for (size_t i_corner = 0; i_corner < nb_corners; i_corner++) {
            //int64_t j = i * nb_corners + i_corner; // explicitly equal to having j++ at the end of the inner loop
            Point<DIM> origin_voxel;
            for (size_t d = 0; d < DIM; d++) {
                int64_t new_vox_coord = 2 * m_voxel_coordinates[i][d] + tabcorner[i_corner][d];
                new_voxel_coordinates[k][d] = new_vox_coord;
                origin_voxel[d] = m_global_origin[d] + new_vox_coord * new_voxel_lengths[d];
            }
            // only retain voxels that are not covered
            if (not auxi::is_covered<DIM>(origin_voxel, new_corners_voxel, a_rsa_grid, a_minimal_radius, list_of_cell_id_sphere_id)) {
                k++;
                if (k >= new_voxel_coordinates.size()) {
                    new_voxel_coordinates.resize(1.5 * new_voxel_coordinates.size());
                }
            }
        }
    }
    new_voxel_coordinates.resize(k);// last coordinates are not used, because of covered_voxels

    set_voxel_lengths(new_voxel_lengths);
    std::swap(m_voxel_coordinates, new_voxel_coordinates);
}

template<int DIM>
inline void list_of_voxels<DIM>::subdivide() {
    if (this->nb_voxels() == 0) { return; }
    m_voxel_lengths *= 0.5;
    //
    int64_t old_size = this->size();
    int64_t nb_corners = auxi_function::puissance<DIM>(2);
    std::vector<DiscPoint<DIM>> new_voxel_coordinates(nb_corners * old_size);
    for (size_t i = 0; i < old_size; i++) {
        const auto& tabcorner = sac_de_billes::path::TabCorner<DIM>::get().getTab();
        for (size_t j = 0; j < nb_corners; j++) {
            for (size_t d = 0; d < DIM; d++) {
                new_voxel_coordinates[i * nb_corners + j][d]
                    = 2 * m_voxel_coordinates[i][d] + tabcorner[j][d];
            }
        }
    }
    swap(m_voxel_coordinates, new_voxel_coordinates);
}

template<int DIM>
arr_vec_double<DIM> list_of_voxels<DIM>::pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const {
    // init ret
    arr_vec_double<DIM> ret{};
    if (this->nb_voxels() == 0 and a_nb_pts != 0) {
        std::cerr << __PRETTY_FUNCTION__ << std::endl;
        throw std::runtime_error("Impossible to draw centers in no voxels");
    }
    ret.resize(a_nb_pts);
    // uniform random distribution on (0, 1)
    law::uniform<double> coord_picker(0, 1, random_generator);
    // uniform random distribution for voxels
    law::uniform<int> voxel_picker(0, this->nb_voxels() - 1, random_generator);
    //
    for (int64_t i = 0; i < a_nb_pts; i++) {
        int64_t id_voxel = voxel_picker();
        Point<DIM> orig_voxel = this->orig_voxel(id_voxel);
        for (int d = 0; d < DIM; d++) {
            ret[i][d] = orig_voxel[d] + this->m_voxel_lengths[d] * coord_picker();
        }
    }
    return ret;
}

template<int DIM>
inline Point<DIM> list_of_voxels<DIM>::corner(int64_t a_id_voxel, int64_t a_index) const {
    Point<DIM> result = this->m_global_origin;
    const auto& discCorner = sac_de_billes::path::TabCorner<DIM>::get().getTab()[a_index];
    for (int d = 0; d < DIM; d++) {
        result[d] += (discCorner[d] + get_disc_coord(a_id_voxel)[d]) * m_voxel_lengths[d];
    }
    return result;
}

template<int DIM>
double list_of_voxels<DIM>::total_area() const {
    return this->nb_voxels() * auxi_function::productOf<double>(m_voxel_lengths);
}

template<int DIM>
inline Point<DIM> list_of_voxels<DIM>::orig_voxel(int64_t a_id_voxel) const {
    return auxi::orig_voxel<DIM>(this->m_global_origin, this->get_disc_coord(a_id_voxel), this->m_voxel_lengths);
}

template<int DIM>
Point<DIM> auxi::orig_voxel(const Point<DIM>& origin,
    const DiscPoint<DIM>& discPoint,
    const Point<DIM>& voxel_lengths) {
    Point<DIM> result = origin;
    for (int d = 0; d < DIM; d++) {
        result[d] += discPoint[d] * voxel_lengths[d];
    }
    return result;
}


template<int DIM>
inline Point<DIM> list_of_voxels<DIM>::center(int64_t a_id_voxel) const {
    Point<DIM> result = this->m_global_origin;
    for (int d = 0; d < DIM; d++) {
        result[d] += (get_disc_coord(a_id_voxel)[d] + 0.5) * m_voxel_lengths[d];
    }
    return result;
}

template<int DIM>
void list_of_voxels<DIM>::print(std::ostream& f) const {
    f << "nb of voxels " << this->size() << std::endl;
    f << "voxel length :"; auxi_function::writeVectorToString(this->m_voxel_lengths, f); f << std::endl;
    for (int64_t i = 0; i < this->size(); i++) {
        for (int d = 0; d < DIM; d++) {
            f << this->get_disc_coord(i)[d] << " ";

        }
        f << endl;
    }
}

template<int DIM>
template<class RSA_GRID>
bool list_of_voxels<DIM>::is_covered(int64_t a_id_voxel,
    const RSA_GRID& a_rsa_grid,
    double a_minimal_radius) const {
    return auxi::is_covered<DIM>(this->orig_voxel(a_id_voxel), m_voxel_lengths,
        this->get_corners_voxel(),
        a_rsa_grid, a_minimal_radius);
}

template<int DIM>
inline DiscPoint<DIM> list_of_voxels<DIM>::get_disc_coord(int64_t a_id_voxel) const {
    return m_voxel_coordinates[a_id_voxel];
}

template<int DIM>
vec_i<DIM> list_of_voxels<DIM>::compute_shift_nb_vox(Point<DIM> a_global_origin,
    Point<DIM> a_origin, Point<DIM> a_cell_length) {
    vec_i<DIM> ret;
    for (int d = 0; d < DIM; d++) {
        ret[d] = std::floor(0.5 + (a_origin[d] - a_global_origin[d]) / a_cell_length[d]);
        if (abs((a_origin[d] - a_global_origin[d]) - ret[d] * a_cell_length[d]) > 1e-9) {
            throw runtime_error("Incoherent origin and global origin");
        }
    }
    return ret;
}

template<int DIM>
vec_i<DIM> list_of_voxels<DIM>::compute_nb_vox(const Point<DIM>& a_domain_length,
    double a_max_diagonal) {
    vec_i<DIM> ret;
    for (int d = 0; d < DIM; d++) {
        ret[d] = 1 + static_cast<int>(sqrt(DIM) * a_domain_length[d] / a_max_diagonal);
    }
    return ret;
}

template<int DIM>
vec_i<DIM> list_of_voxels<DIM>::compute_nb_vox(const Point<DIM>& a_domain_length,
    const Point<DIM>& a_cell_length, double a_max_diagonal) {
    vec_i<DIM> nb_vox_cell = list_of_voxels<DIM>::compute_nb_vox(a_cell_length, a_max_diagonal);
    vec_i<DIM> multiplicator;
    vec_i<DIM> ret;
    for (int d = 0; d < DIM; d++) {
        multiplicator[d] = static_cast<int>(floor(a_domain_length[d] / a_cell_length[d] + 0.5));
        if (abs(multiplicator[d] * a_cell_length[d] - a_domain_length[d]) > 1e-9) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Incompatible cell length and domain length");
        }
        ret[d] = nb_vox_cell[d] * multiplicator[d];
    }
    return ret;
}

template<int DIM>
Point<DIM> list_of_voxels<DIM>::compute_voxel_length(const Point<DIM>& a_domain_length,
    double a_max_diagonal) {
    Point<DIM> ret;
    auto my_nb_vox = list_of_voxels<DIM>::compute_nb_vox(a_domain_length, a_max_diagonal);
    for (int d = 0; d < DIM; d++) {
        ret[d] = a_domain_length[d] / my_nb_vox[d];
    }
    return ret;
}

template<int DIM>
void list_of_voxels<DIM>::set_voxel_lengths(const Point<DIM>& voxel_lengths) {
    m_voxel_lengths = voxel_lengths;
    auxi::create_corners_voxel_inplace<DIM>(voxel_lengths, m_corners_voxel);
}

//
template<int DIM, class RSA_DATA_STORAGE>
double auxi::distance_squared(const Point<DIM>& pt0, const RSA_DATA_STORAGE& rsa_data_storage, int64_t id) {
    return geomTools::distanceCarre<DIM>(rsa_data_storage.get_center(id), pt0);
}

template<int DIM, class VectorOfPoints>
inline bool auxi::is_in_sphere(const VectorOfPoints& corner_voxels,
    const Point<DIM>& delta_center, double a_sphere_radius, double a_minRadius) {
    static_assert(std::is_same < std::array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>, VectorOfPoints>::value);
    // if this assert is false, should change the limit of the loop
    double sphereInfluenceRadiusSquared = auxi_function::puissance<2>(a_sphere_radius + a_minRadius);
    for (int i = 0; i < auxi_function::puissance<DIM>(2); i++) {
        double distanceSquared = sac_de_billes::geomTools::distanceCarre<DIM>(corner_voxels[i], delta_center);
        if (distanceSquared > sphereInfluenceRadiusSquared) {
            return false;
        }
    }
    return true;
}

template<int DIM, class RSA_GRID>
bool auxi::is_covered(const Point<DIM>& origin_voxel, const Point<DIM>& voxel_lengths,
    const array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corners_voxel,
    const RSA_GRID& rsa_grid,
    double a_minimal_radius) {
    Point<DIM> center = origin_voxel + 0.5 * voxel_lengths;
    int64_t cell_idx = rsa_grid.compute_cell_idx(center);
    auto& nc_list = rsa_grid.get_list_of_neighbor_cells();
    // corners_voxels
    for (size_t nc = 0; nc < nc_list.size(); nc++) {
        const auto& spheres = rsa_grid.get_data(cell_idx + nc_list[nc]);
        for (size_t s_id = 0; s_id < spheres.size(); s_id++) {
            Point<DIM> delta_center = spheres.get_center(s_id) - origin_voxel;
            if (auxi::is_in_sphere<DIM>(corners_voxel, delta_center, spheres.get_rad(s_id), a_minimal_radius)) {
                return true;
            }
        }
    }
    return false;
}

template<int DIM, class RSA_GRID>
bool auxi::is_covered(const Point<DIM>& origin_voxel,
    const array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corners_voxel,
    const RSA_GRID& rsa_grid,
    double a_minimal_radius,
    const vector<pair<int, int>>& list_of_cell_id_sphere_id) {
    if (list_of_cell_id_sphere_id.size() == 0) {
        return false;
    }
    //
    for (size_t i = 0; i < list_of_cell_id_sphere_id.size(); i++) {
        const auto& spheres = rsa_grid.get_data(list_of_cell_id_sphere_id[i].first);
        auto s_id = list_of_cell_id_sphere_id[i].second;
        Point<DIM> delta_center = spheres.get_center(s_id) - origin_voxel;
        double sphere_radius = spheres.get_rad(s_id);
        if (auxi::is_in_sphere<DIM>(corners_voxel, delta_center, sphere_radius, a_minimal_radius)) {
            return true;
        }
    }
    return false;
}

template<int DIM, class RSA_GRID>
void auxi::get_list_of_spheres(const Point<DIM>& center_voxel, const Point<DIM>& voxel_lengths,
    const RSA_GRID& rsa_grid, double a_minimal_radius, vector<pair<int, int>>& res) {
    int64_t cell_idx = rsa_grid.compute_cell_idx(center_voxel);
    auto& nc_list = rsa_grid.get_list_of_neighbor_cells();
    constexpr int64_t max_size = sac_de_billes::auxi_function::puissance<DIM>(2) * sac_de_billes::auxi_function::puissance<DIM>(3);
    res.resize(max_size);
    int64_t size_res = 0;
    //
    for (size_t nc = 0; nc < nc_list.size(); nc++) {
        int64_t cell_id = cell_idx + nc_list[nc];
        const auto& spheres = rsa_grid.get_data(cell_id);
        for (size_t s = 0; s < spheres.size(); s++) {
            double sphere_radius = spheres.get_rad(s);
            if (auxi::distance_squared<DIM>(center_voxel, spheres, s)
                < auxi_function::puissance<2>(sphere_radius + a_minimal_radius)) {
                res[size_res].first = cell_id;
                res[size_res].second = s;
                size_res++;
            }
        }
    }
    //
    res.resize(size_res);
}

template<int DIM>
void auxi::create_corners_voxel_inplace(const Point<DIM>& voxel_lengths,
    array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corners_voxel) {
    const auto& tabCorner = sac_de_billes::path::TabCorner<DIM>::get().getTab();
    for (int i = 0; i < auxi_function::puissance<DIM>(2); i++) {
        for (int d = 0; d < DIM; d++) {
            corners_voxel[i][d] = tabCorner[i][d] * voxel_lengths[d];
        }
    }
}

template<int DIM>
template<class FUNCTION, class TEST_ID>
void list_of_voxels<DIM>::apply_function_depending_on_absolute_cell_idx(FUNCTION function_to_be_applied, TEST_ID test,
    const auto& rsa_grid_traversal_) {
    for (int64_t id_voxel = 0; id_voxel < nb_voxels(); id_voxel++) {
        if (test(compute_absolute_cell_idx(id_voxel, rsa_grid_traversal_))) {
            function_to_be_applied(id_voxel);
        }
    }
}

template<int DIM>
template<class TEST_ID>
void list_of_voxels<DIM>::pop_voxels_depending_on_absolute_cell_idx(TEST_ID test,
    vector<VoxelCoordinates<DIM>>& selected_voxel_coordinates,
    const auto& rsa_grid_traversal_) {
    auto pop_single_voxel = [&](const auto id_voxel) {
        if (test(compute_absolute_cell_idx(id_voxel, rsa_grid_traversal_))) {
            selected_voxel_coordinates.push_back(this->m_voxel_coordinates[id_voxel]);
            return true;
        }
        return false;
        };
    remove_if(pop_single_voxel);
}

template<int DIM>
void list_of_voxels<DIM>::insert_voxels(
    const vector<VoxelCoordinates<DIM>>& voxel_absolute_coordinates) {
    this->m_voxel_coordinates.insert(this->m_voxel_coordinates.end(),
        voxel_absolute_coordinates.begin(), voxel_absolute_coordinates.end());
}

} // namespace  voxel_list
