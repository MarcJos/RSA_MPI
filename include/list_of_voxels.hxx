//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
/*
For implementing the voxel approach.
Voxels that can be subdivided.

Largely inspired from Merope
*/

#include <random>

#include "link_to_auxi.hxx"


namespace voxel_list {
using namespace sac_de_billes;
using namespace std;

//! @brief Class storing the voxels that may be not covered by existing spheres.
//! @tparam DIM : space dimension
template<int DIM>
class list_of_voxels {
    static constexpr int local_method = 1;
    // (Refers to removed code)
    // 1 better, then 0, then 2
    //! =0, test all candidate spheres when subdividing
    //! =1, test only spheres intersecting the center of the voxel (faster)
    //! =2, when suppressing and subdividing, only remove voxels by only testing their centers
    //! namely compute verify distance(center_voxel, center_sphere) + half_diagonal < radius_spheres + minimal_radius
private:
    //! @brief : storage for voxel coordinates (in 3D ) [i_1, j_1, k_1, i_2, ...] for
    //! i_m, j_m, k_m being the integer coordinates of the voxels
    //! the real coordinates are i_m * m_voxel_lengths[0], j_m * m_voxel_lengths[1], k_m * m_voxel_lengths[2]
    //! the 0 is considered to be in the domain (local coordinates, thus)
    vector<array<int64_t, DIM>> m_voxel_coordinates;
    //! @brief : common size of each cell
    Point<DIM> m_voxel_lengths;
    //! coordinates of the origin of the domain
    const Point<DIM> m_origin;
    //! corners of voxels (the origin is at 0)
    array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)> m_corners_voxel;

public:
    //! @brief constructor
    //! @param a_origin : origin of the domain
    //! @param a_domain_length : dimension of the domain
    //! @param a_minradius : all voxels should have diagonal lesser than a_max_diagonal
    list_of_voxels(Point<DIM> a_origin, Point<DIM> a_domain_length, double a_max_diagonal);
    //! @brief constructor
    //! @param a_origin : origin of the domain
    //! @param a_nb_vox : nb of voxels in each DIM direction
    //! @param a_voxel_length : voxel_length
    list_of_voxels(Point<DIM> a_origin, vec_i<DIM> a_nb_vox, Point<DIM> a_voxel_length);
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
    size_t size() const;
    //! @return the area of all the voxels
    double total_area() const;

    // random picking

    //! @return : a_nb_pts randomly picked in the voxels (with uniform law)
    arr_vec_double<DIM> pick_points(int64_t a_nb_pts, std::mt19937& random_generator) const;

    // Individual function for each voxel

    //! \return the corner coordinates
    //! \param a_index in {0, 1}^d
    Point<DIM> corner(int64_t a_id_voxel, int64_t a_index) const;
    //! coordinates of the closest corner to 0
    Point<DIM> orig_voxel(int64_t a_id_voxel) const;
    //! coordinate of the center of the voxel
    Point<DIM> center(int64_t a_id_voxel) const;

    //! @return : if the voxel is covered by spheres stored in a_rsa_data_storage and pointed-to by a_rsa_grid
    //! the voxel is covered if covered by the sphere, the radius of which is incremented by a_minimal_radius
    //! @param a_id_voxel
    //! @param a_rsa_grid : implicitly, rsa_grid<DIM>
    //! @param a_minimal_radius : minimal radius of the rsa simulation
    template<class RSA_GRID>
    bool is_covered(int64_t a_id_voxel,
        const RSA_GRID& rsa_cell, double a_minimal_radius) const;

    void print(std::ostream& f) const;

private:
    //! @return : the discrete coordinates of the voxel. [0, 0, 0] is the left lowest corner of the rsa_domain
    DiscPoint<DIM> get_disc_coord(int64_t a_id_voxel) const;
    //! @return : the ideal number of voxels is each direction, given the parameters
    //! @param a_domain_length : total length of the domain
    //! @param a_max_diagonal : each voxel should have diagonal < a_max_diagonal
    static vec_i<DIM> compute_nb_vox(const Point<DIM>& a_domain_length, double a_max_diagonal);
    //! @return : the ideal voxel length given the parameters
    //! @param a_domain_length : total length of the domain
    //! @param a_max_diagonal : each voxel should have diagonal < a_max_diagonal
    static Point<DIM> compute_voxel_length(const Point<DIM>& a_domain_length, double a_max_diagonal);
    //! @brief : const getter
    const array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& get_corners_voxel() const { return m_corners_voxel; }
    //! @brief : setter
    void set_voxel_lengths(const Point<DIM>& voxel_lengths);
};

namespace auxi {

//! @return : if the voxel is covered by spheres stored in a_rsa_data_storage and pointed-to by a_rsa_grid
//! the voxel is covered if covered by the sphere, the radius of which is incremented by a_minimal_radius
//! @param origin_voxel, voxel_lengths : voxel = [origin_voxel[0] + voxel_lengths] x [ ...]
//! @param a_rsa_grid : implicitly, rsa_grid<DIM>
//! @param a_minimal_radius : minimal radius of the rsa simulation
template<int DIM, class RSA_GRID>
bool is_covered(const Point<DIM>& origin_voxel, const Point<DIM>& voxel_lengths,
    const array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corners_voxel,
    const RSA_GRID& a_rsa_grid, double a_minimal_radius);
template<int DIM, class RSA_GRID>
bool is_covered(const Point<DIM>& origin_voxel,
    const array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corners_voxel,
    const RSA_GRID& a_rsa_grid, double a_minimal_radius,
    const vector<pair<int, int>>& list_of_cell_id_sphere_id);

//! @return squared distance between a point and the center of a spehre relative to rsa_data_storage
//! @param pt0 : point
//! @param rsa_data_storage : rsa_data_storage<DIM>, stores spheres
//! @param id : index of the sphere in the rsa_data_storage
template<int DIM, class RSA_DATA_STORAGE>
double distance_squared(const Point<DIM>& pt0, const RSA_DATA_STORAGE& rsa_data_storage, int64_t id);

//! @return whether it is totally covered by the sphere or not
//! @param origin_voxel, voxel_lengths : voxel = [origin_voxel[0] + voxel_lengths] x [ ...]
//! @param a_sphere : given sphere
//! @param a_minRadius : minimal radius of all spheres
template<int DIM, class VectorOfPoints>
bool is_in_sphere(const VectorOfPoints& corner_voxels,
    const Point<DIM>& delta_center, double a_sphere_radius, double a_minRadius);

//! @return : le identifiers (cell_id, sphere_id) of spheres that may cover the subvoxels
//! @param a_rsa_grid : implicitly, rsa_grid<DIM>
//! @param a_minimal_radius : minimal radius of the rsa simulation
template<int DIM, class RSA_GRID>
void get_list_of_spheres(const Point<DIM>& center_voxel, const Point<DIM>& voxel_lengths,
    const RSA_GRID& a_rsa_grid, double a_minimal_radius, vector<pair<int, int>>& result);
//! @brief reture in-place in corner_voxels
//! the array of the corners of a voxel
//! @param origin_voxel : origin of the voxel
//! @param voxel_lengths : length of the voxel
//! @param corner_voxels : will store the result
template<int DIM>
void create_corners_voxel_inplace(const Point<DIM>& voxel_lengths,
    array<Point<DIM>, sac_de_billes::auxi_function::puissance<DIM>(2)>& corner_voxels);
} // namespace  auxi


} // namespace  voxel_list

#include<list_of_voxels.ixx>
