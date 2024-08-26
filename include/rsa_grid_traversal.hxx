//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once


#include <basic_types.hxx>

enum class TypeCell {
    Edge, Ghost, Real, All
};

//! @brief : rsa_grid is exclusively here to store a list of cell, and make sense of traversal. It does not store any data **inside** the cell
//! pave the domain delimitated by m_shift and m_shift + m_cell_length.
//! this paving strictly encloses the domain delimitated by a_inf, a_sup
//! restricting this paving to the inner cells exactly paves the domain delimitated by a_inf, a_sup

template<int DIM>
class rsa_grid_traversal {
private:
    std::array<int, sac_de_billes::auxi_function::puissance<DIM>(3)> m_neighbor_shift; ///< shifted values
    int m_ghost_layer; ///< number of layers of 2*radius size
    vec_i<DIM> m_n; ///< number of cells i j k
    vec_d<DIM> m_shift; ///< shift to local ref
    vec_d<DIM> m_cell_length; ///< cell size
    vec_d<DIM> m_cell_inverse_length; ///< cell size
    std::vector<uint64_t> m_ghost_traversal; ///< list of ghost cells
    std::vector<uint64_t> m_edge_traversal; ///< list of edge cells
    std::vector<uint64_t> m_real_traversal; ///< list of real cells


public:
    //! @brief : default constructor
    rsa_grid_traversal() {}
    //! @brief : this constructor calls the function "define".
    //! @see : define.
    rsa_grid_traversal(const double a_rad, const int a_ghost_layer,
        const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup);

    //! @brief : get the total number of cells
    int size() const;
    //! @brief : get the list of neighbor cell shift values
    const std::array<int, sac_de_billes::auxi_function::puissance<DIM>(3)>& get_list_of_neighbor_cells() const { return m_neighbor_shift; }
    //! @brief : This function fills the interval [being,end]
    //! according to the "inside" domain length
    void get_interval_inside(const int a_dim, double& begin, double& end) const;
    //! getter
    const vec_d<DIM>& get_cell_length() const { return m_cell_length; }

    template<TypeCell typeCell>
    bool is_type_cell(int idx) const;

    //! @brief : This function return true if the cell is a ghost cell
    bool is_ghost(const vec_i<DIM>& idx) const;
    //! @brief : This function return true if the cell is a ghost cell
    bool is_ghost(const int idx) const;
    //! @brief : This function return true if the cell is a edge cell
    bool is_edge(const vec_i<DIM>& idx) const;
    //! @brief : This function return true if the cell is a real cell
    bool is_real(const vec_i<DIM>& idx) const;
    //! @brief : This function return true if the cell is a real cell
    bool is_real(vec_i<DIM>&& idx) const;
    //! @brief : This function return true if the cell is in the cell list
    bool is_cell(const vec_i<DIM>& idx) const;
    //! @brief : This function return true if the cell is in the cell list
    bool is_cell(vec_i<DIM>&& idx) const;
    //! @brief : This function return true if the cell is in the cell list
    bool is_cell(const int idx) const;

    //! @brief : This function computes the cell index from the cartesian position of the cell
    // DIM 3 : idx = i + j * nx + k * nx * ny
    uint64_t index(const vec_i<DIM>& a_ijk) const;
    //! @brief : This function reverses the index function
    vec_i<DIM> get_coordinates(const int a_idx) const;
    //! @brief : This function computes the cell index from the cartesian position of the particle
    // DIM 3 : idx -> ijk (global ref) + apply shift in local ref
    uint64_t compute_cell_idx(const vec_d<DIM>& a_pos) const;

protected:
    //! @brief  apply the Lambda for each cell index of desired type
    //! @tparam function : function depending on the index, to be applied
    //! @tparam typeCell : TypeCell for traversal
    template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
    void apply_on_cell_id(LambdaFunction& function, Arrays&... arrays) const;
    template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
    void apply_on_cell_id(LambdaFunction& function, Arrays&... arrays);

    template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
    void apply_on_chosen_cells_id(const std::vector<uint64_t>& trasversed_cells, LambdaFunction& function, Arrays&... arrays) const;
    template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
    void apply_on_chosen_cells_id(const std::vector<uint64_t>& trasversed_cells, LambdaFunction& function, Arrays&... arrays);

private:
    //! @brief : getter for the member m_ghost_traversal / m_edge_traversal  / real_traversal
    template<TypeCell typeCell>
    const std::vector<uint64_t>& get_traversal() const;

    //! @brief : construct the object
    void define(const double a_rad, const int a_ghost_layer,
        const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup);
    //! @brief : This function returns (nx * ny * nz)
    int compute_size(const vec_i<DIM>& a_n) const;
    //! @brief : computes ghost cell indexes and store them in m_ghost_traversal
    void build_ghost_traversal();
    //! @brief : computes edge cell indexes and store them in m_edge_traversal
    void build_edge_traversal();
    //! @brief : computes real cell indexes and store them in m_real_traversal
    void build_real_traversal();
};

#include<rsa_grid_traversal.ixx>
