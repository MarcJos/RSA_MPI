//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM>
rsa_grid_traversal<DIM>::rsa_grid_traversal(const double a_rad, const int a_ghost_layer,
    const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup) {
    define(a_rad, a_ghost_layer, a_inf, a_sup);
}

template<int DIM>
void rsa_grid_traversal<DIM>::define(const double a_rad, const int a_ghost_layer, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup) {
    m_ghost_layer = a_ghost_layer;

    // We define : cell size (m_cell_length, double), number of cells (m_n, int), shift = lower boundary (m_shift, double)
    for (int d = 0; d < DIM; d++) {
        double length = a_rad * 2; // 2 because the minimal distance between two centers is rad1 + rad2
        double inside_domain_width = a_sup[d] - a_inf[d];
        int nb_cells_inside = static_cast<int>(inside_domain_width / length);
        if (nb_cells_inside == 0) {
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
            throw std::runtime_error("Too small domain wrt ghost layer");
        }
        m_cell_length[d] = inside_domain_width / nb_cells_inside;
        m_cell_inverse_length[d] = 1. / m_cell_length[d];
        m_shift[d] = a_inf[d] - m_ghost_layer * m_cell_length[d];
        m_n[d] = nb_cells_inside + 2 * m_ghost_layer;
    }

    // Function used to figure out the neighbor cell shift values.
    // dim[0] = 1, dim[1] = dx, dim[2]=dx*dy ..., dim[n]=d1*d2*...*dn
    auto shift_dim = [](const int a_shift, const int a_dim, const vec_i<DIM>& a_n) -> int {
        assert(DIM > a_dim);
        int ret = a_shift;
        for (int d = a_dim - 1; d >= 0; d--) {
            ret *= a_n[d];
        }
        return ret;
        };

    // This section computes the shift values to figure out the neighbor cells index such as :
    // List of neighbor cells (idx) = {idx + m_neighbor_shift[0] , ... , idx + m_neighbor_shift[pow(3,DIM)-1]}
    std::vector<int> indices_for_m_neighbor_shift(sac_de_billes::auxi_function::puissance<DIM>(3));
    std::vector<std::array<double, DIM>> vectors_for_m_neighbor_shift(sac_de_billes::auxi_function::puissance<DIM>(3));
    for (int i = 0; i < sac_de_billes::auxi_function::puissance<DIM>(3); i++) {
        indices_for_m_neighbor_shift[i] = i;
        int tmp = i;
        int r = 0, val = 0;
        for (int d = 0; d < DIM; d++) {
            r = tmp % 3;
            val += shift_dim(r - 1, d, m_n);
            vectors_for_m_neighbor_shift[i][d] = r - 1;
            tmp = (tmp - r) / 3;
        }
        m_neighbor_shift[i] = val;
    }

    std::sort(indices_for_m_neighbor_shift.begin(), indices_for_m_neighbor_shift.end(),
        [&vectors_for_m_neighbor_shift](auto i1, auto i2) {
            return sac_de_billes::geomTools::normeCarre<DIM>(vectors_for_m_neighbor_shift[i1]) <
                sac_de_billes::geomTools::normeCarre<DIM>(vectors_for_m_neighbor_shift[i2]);
        });
    auto copy_of_m_neighbor_shift = m_neighbor_shift;
    for (int i = 0; i < sac_de_billes::auxi_function::puissance<DIM>(3); i++) {
        m_neighbor_shift[i] = copy_of_m_neighbor_shift[indices_for_m_neighbor_shift[i]];
    }
    // Warning. Better ordering for neighbors. Not perfect!

    // build traversals
    build_ghost_traversal();
    build_edge_traversal();
    build_real_traversal();
}


template<int DIM>
int rsa_grid_traversal<DIM>::compute_size(const vec_i<DIM>& a_n) const {
    int ret = 1;
    for (int d = 0; d < DIM; d++) {
        ret *= a_n[d];
    }
    return ret;
}

template<int DIM>
void rsa_grid_traversal<DIM>::build_ghost_traversal() {
    int ncells = this->size();
    for (int idx = 0; idx < ncells; idx++) {
        auto coord = get_coordinates(idx);
        if (is_ghost(coord)) {
            m_ghost_traversal.push_back(idx);
        }
    }
}

template<int DIM>
void rsa_grid_traversal<DIM>::build_edge_traversal() {
    int ncells = this->size();
    for (int idx = 0; idx < ncells; idx++) {
        auto coord = get_coordinates(idx);
        if (is_edge(coord)) {
            m_edge_traversal.push_back(idx);
        }
    }
}


template<int DIM>
template<TypeCell typeCell>
const std::vector<uint64_t>& rsa_grid_traversal<DIM>::get_traversal() const {
    if constexpr (typeCell == TypeCell::Real) {
        return m_real_traversal;
    } else if constexpr (typeCell == TypeCell::Edge) {
        return m_edge_traversal;
    } else  if constexpr (typeCell == TypeCell::Ghost) {
        return m_ghost_traversal;
    }
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
void rsa_grid_traversal<DIM>::apply_on_cell_id(LambdaFunction& function, Arrays&... arrays) {
    if constexpr (typeCell == TypeCell::All) {
        if constexpr (use_omp) {
#pragma omp parallel for
            for (int id_cell = 0; id_cell < this->size(); id_cell++) {
                function(id_cell, (arrays[id_cell])...);
            }
        } else {
            for (int id_cell = 0; id_cell < this->size(); id_cell++) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    } else {
        const auto& trasversed_cells = this->template get_traversal<typeCell>();
        apply_on_chosen_cells_id<TypeCell::All, use_omp>(trasversed_cells, function, arrays...);
    }
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
void rsa_grid_traversal<DIM>::apply_on_cell_id(LambdaFunction& function,
    Arrays&... arrays) const {
    if constexpr (typeCell == TypeCell::All) {
        if constexpr (use_omp) {
#pragma omp parallel for
            for (int id_cell = 0; id_cell < this->size(); id_cell++) {
                function(id_cell, (arrays[id_cell])...);
            }
        } else {
            for (int id_cell = 0; id_cell < this->size(); id_cell++) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    } else {
        const auto& trasversed_cells = this->template get_traversal<typeCell>();
        apply_on_chosen_cells_id<TypeCell::All, use_omp>(trasversed_cells, function, arrays...);
    }
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
void rsa_grid_traversal<DIM>::apply_on_chosen_cells_id(const std::vector<uint64_t>& trasversed_cells,
    LambdaFunction& function, Arrays&... arrays) const {
    const int size = trasversed_cells.size();
    if constexpr (use_omp) {
#pragma omp parallel for
        for (int i_id_cell = 0; i_id_cell < size; i_id_cell++) {
            int id_cell = trasversed_cells[i_id_cell];
            if constexpr (typeCell == TypeCell::All) {
                function(id_cell, (arrays[id_cell])...);
            } else if (is_type_cell<typeCell>(id_cell)) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    } else {
        for (int i_id_cell = 0; i_id_cell < size; i_id_cell++) {
            int id_cell = trasversed_cells[i_id_cell];
            if constexpr (typeCell == TypeCell::All) {
                function(id_cell, (arrays[id_cell])...);
            } else if (is_type_cell<typeCell>(id_cell)) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    }
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction, class ...Arrays>
void rsa_grid_traversal<DIM>::apply_on_chosen_cells_id(const std::vector<uint64_t>& trasversed_cells,
    LambdaFunction& function, Arrays&... arrays) {
    const int size = trasversed_cells.size();
    if constexpr (use_omp) {
#pragma omp parallel for
        for (int i_id_cell = 0; i_id_cell < size; i_id_cell++) {
            int id_cell = trasversed_cells[i_id_cell];
            if constexpr (typeCell == TypeCell::All) {
                function(id_cell, (arrays[id_cell])...);
            } else if (is_type_cell<typeCell>(id_cell)) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    } else {
        for (int i_id_cell = 0; i_id_cell < size; i_id_cell++) {
            int id_cell = trasversed_cells[i_id_cell];
            if constexpr (typeCell == TypeCell::All) {
                function(id_cell, (arrays[id_cell])...);
            } else if (is_type_cell<typeCell>(id_cell)) {
                function(id_cell, (arrays[id_cell])...);
            }
        }
    }
}

template<int DIM>
void rsa_grid_traversal<DIM>::build_real_traversal() {
    int ncells = this->size();
    for (int idx = 0; idx < ncells; idx++) {
        auto coord = get_coordinates(idx);
        if (!is_ghost(coord)) {
            m_real_traversal.push_back(idx);
        }
    }

#ifdef PRINT_ALL_DEBUG
    if (DIM == 2) {
        std::cout << " X " << m_n[0] << " Y " << m_n[1] << std::endl;
        for (auto it : m_real_traversal)
            std::cout << it << " ";
        std::cout << std::endl;
    }
#endif
}

template<int DIM>
template<TypeCell typeCell>
bool rsa_grid_traversal<DIM>::is_type_cell(int idx) const {
    if constexpr (typeCell == TypeCell::All) {
        return true;
    }
    //
    auto coord = this->get_coordinates(idx);
    if constexpr (typeCell == TypeCell::Ghost) {
        return is_ghost(coord);
    } else if constexpr (typeCell == TypeCell::Edge) {
        return is_edge(coord);
    } else if constexpr (typeCell == TypeCell::Real) {
        return is_real(coord);
    } else {
        static_assert(typeCell == TypeCell::All);
    }
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_ghost(const int idx) const {
    auto coord = this->get_coordinates(idx);
    return is_ghost(coord);
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_ghost(const vec_i<DIM>& idx) const {
    if (!is_cell(idx)) {
        return false;
    }
    for (int d = 0; d < DIM; d++) {
        if (idx[d] < m_ghost_layer) return true;
        if (idx[d] > m_n[d] - 1 - m_ghost_layer) return true;
    }
    return false;
}

template<int DIM>
int rsa_grid_traversal<DIM>::size() const {
    return compute_size(m_n);
}

template<int DIM>
void rsa_grid_traversal<DIM>::get_interval_inside(const int a_dim, double& begin, double& end) const {
    const int layer = 1 + m_ghost_layer; // (edge+ghost)
    begin = m_shift[a_dim] + layer * m_cell_length[a_dim];
    end = m_shift[a_dim] + (m_n[a_dim] - layer) * m_cell_length[a_dim];
    std::cout << " dim " << a_dim << " begin " << begin << " end " << end << std::endl;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_edge(const vec_i<DIM>& idx) const {
    for (int d = 0; d < DIM; d++) {
        if (idx[d] == m_ghost_layer) return true; // lower edge side
        if (idx[d] == m_n[d] - 1 - m_ghost_layer) return true; // upper edge side
    }
    return false;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_real(const vec_i<DIM>& idx) const {
    for (int d = 0; d < DIM; d++) {
        if (idx[d] < m_ghost_layer) return false; // lower edge side
        if (idx[d] > m_n[d] - 1 - m_ghost_layer) return false; // upper edge side
    }
    return true;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_real(vec_i<DIM>&& idx) const {
    for (int d = 0; d < DIM; d++) {
        if (idx[d] < m_ghost_layer) return false; // lower edge side
        if (idx[d] > m_n[d] - 1 - m_ghost_layer) return false; // upper edge side
    }
    return true;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_cell(const int idx) const {
    if (idx < 0) return false; // lower side
    if (idx >= this->size()) return false; // upper side
    return true;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_cell(const vec_i<DIM>& idx) const {
    for (int d = 0; d < DIM; d++) {
        if (idx[d] < 0) return false; // lower side
        if (idx[d] >= m_n[d]) return false; // upper side
    }
    return true;
}

template<int DIM>
bool rsa_grid_traversal<DIM>::is_cell(vec_i<DIM>&& idx) const {
    for (int d = 0; d < DIM; d++) {
        if (idx[d] < 0) return false; // lower side
        if (idx[d] >= m_n[d]) return false; // upper side
    }
    return true;
}

template<int DIM>
uint64_t rsa_grid_traversal<DIM>::index(const vec_i<DIM>& a_ijk) const {
    uint64_t ret = 0;
    for (int d = DIM - 1; d >= 0; d--) {
        ret = ret * m_n[d] + a_ijk[d];
    }
    return ret;
}

template <int DIM>
int mult(const int N, const vec_i<DIM>& a_in) {
    assert(N < DIM);
    int ret = 1;
    for (int i = 0; i < N; i++) {
        ret *= a_in[i];
    }
    return ret;
}

template<int DIM>
void compute_coordinates(const int a_dim, const int a_idx, vec_i<DIM>& a_res, const vec_i<DIM>& a_sizes) {
    // could be optimized
    if (a_dim < 0) return;
    assert(a_idx >= 0);
    int size_dim = mult<DIM>(a_dim, a_sizes);
    a_res[a_dim] = a_idx / size_dim;
    compute_coordinates<DIM>(a_dim - 1, a_idx - a_res[a_dim] * size_dim, a_res, a_sizes);
}

template<int DIM>
vec_i<DIM> rsa_grid_traversal<DIM>::get_coordinates(const int a_idx) const {
    vec_i<DIM> ret;
    compute_coordinates<DIM>(int(DIM) - 1, a_idx, ret, m_n);
    return ret;
}

template<int DIM>
uint64_t rsa_grid_traversal<DIM>::compute_cell_idx(const vec_d<DIM>& a_pos) const {
    uint64_t ret = 0;
    for (int d = DIM - 1; d >= 0; d--) {
        int64_t ijk_d = (a_pos[d] - m_shift[d]) * m_cell_inverse_length[d];
        assert(ijk_d >= 0);
        ret = ret * m_n[d] + ijk_d;
    }
    return ret;
}
