//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM>
template<typename... Args>
bool rsa_ghost_area<DIM>::add_particle(Args&&... a_args) {
    if (is_ghost(std::forward<Args>(a_args)...)) {
        m_data.add_sphere(std::forward<Args>(a_args)...);
        return true;
    }
    return false;
}

template<int DIM>
template<typename ...Args>
bool rsa_ghost_area<DIM>::is_ghost(Args && ...a_args) {
    bool ret = check_interval(get_inf(), get_sup(), get_center<DIM>(std::forward<Args>(a_args)...));
    return ret;
}

template<int DIM>
void rsa_ghost_area<DIM>::apply_periodicity() {
    if (!m_periodic) return;
    for (int dim = 0; dim < DIM; dim++) {
        const double shift = m_shift[dim];
        if (std::abs(shift) < 1E-16) continue;
        m_data.shift_centers(dim, shift);
    }
}
