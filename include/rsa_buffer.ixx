//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<typename... Args>
rsa_buffer<Args...>::rsa_buffer(int a_proc_id, int a_size) {
    assert(a_proc_id >= 0);
    assert(a_size >= 0);
    set_proc_id(a_proc_id);
    set_size(a_size);
}

template<typename... Args>
void rsa_buffer<Args...>::set_size(size_t a_size) {
    resize(a_size);
}

template<typename... Args>
void rsa_buffer<Args...>::resize(size_t a_size) {
    m_size = a_size;
    int bytesize = m_size * args_size();
    m_data.resize(bytesize);
}

template<typename... Args>
void rsa_buffer<Args...>::purge() {
    m_size = 0;
    m_data.clear();
}

//! pack

template<typename... Args>
void rsa_buffer<Args...>::pack(int a_first, int a_last, const Args*... m_args) {
    int shift = 0;
    auto& info = get_data();
    assert(a_last - a_first <= size() && " this buffer is not correctly resized");
    pack_impl(info.data(), a_first, a_last, shift, m_args...);
}

template<typename... Args>
void rsa_buffer<Args...>::pack(int a_size, const Args*... m_args) {
    pack(0, a_size, m_args...);
}


template<typename... Args>
template<typename U>
void rsa_buffer<Args...>::pack_impl(char* a_data, int a_first, int a_last, int a_shift, const U* a_ptr) {
    const char* cast_ptr = reinterpret_cast<const char*>(a_ptr);
    int size = (a_last - a_first) * sizeof(U);
    std::copy(cast_ptr, cast_ptr + size, a_data + a_shift);
}

template<typename... Args>
template<typename U, typename... miniArgs>
void rsa_buffer<Args...>::pack_impl(char* a_data, int a_first, int a_last, int a_shift, const U* a_ptr, const miniArgs*... a_args) {
    pack_impl(a_data, a_first, a_last, a_shift, a_ptr);
    int shift = a_shift + size() * sizeof(U);
    pack_impl(a_data, a_first, a_last, shift, a_args...);
}

// unpacking step
template<typename... Args>
void rsa_buffer<Args...>::unpack(int a_first, int a_last, Args*... m_args) const {
    int shift = 0;
    const auto& info = get_data();
    const char* ptr = info.data();
    assert(a_last - a_first <= size());
    unpack_impl(ptr, a_first, a_last, shift, m_args...);
}

template<typename... Args>
void rsa_buffer<Args...>::unpack(int a_size, Args*... m_args) const {
    unpack(0, a_size, m_args...);
}

template<typename... Args>
template<typename U>
void rsa_buffer<Args...>::unpack_impl(const char* a_data, int a_first, int a_last, int a_shift, U* a_ptr) const {
    char* cast_ptr = reinterpret_cast<char*>(a_ptr);
    std::copy(a_data + a_shift + a_first * sizeof(U), a_data + a_shift + a_last * sizeof(U), cast_ptr);
}


template<typename... Args>
template<typename U, typename... miniArgs>
void rsa_buffer<Args...>::unpack_impl(const char* a_data, int a_first, int a_last, int a_shift, U* a_ptr, miniArgs*... a_args) const {
    unpack_impl(a_data, a_first, a_last, a_shift, a_ptr);
    int shift = a_shift + size() * sizeof(U);
    unpack_impl(a_data, a_first, a_last, shift, a_args...);
}


template<typename... Args>
template<typename ... MpiInfo>
void rsa_buffer<Args...>::send_size(MpiInfo&&... a_info) {
    assert(get_proc_id() >= 0);
    MPI_Isend(&m_size, 1, MPI_INT, m_proc_id, a_info ...);
}


template<typename... Args>
template<typename ... MpiInfo>
void rsa_buffer<Args...>::recv_size(MpiInfo&&... a_info) {
    assert(get_proc_id() >= 0);
    //			std::cout << " recv_size("<<m_proc_id<<","<<m_size<<")" <<std::endl;
    MPI_Irecv(&m_size, 1, MPI_INT, m_proc_id, a_info ...);
}

template<typename... Args>
template<typename ... MpiInfo>
void rsa_buffer<Args...>::send_data(MpiInfo&&... a_info) {
    assert(get_proc_id() >= 0);
    char* ptr = m_data.data();
    int size = m_data.size();
    MPI_Isend(ptr, size, MPI_CHAR, m_proc_id, a_info ...);
}

template<typename... Args>
template<typename ... MpiInfo>
void rsa_buffer<Args...>::recv_data(MpiInfo&&... a_info) {
    assert(get_proc_id() >= 0);
    resize(m_size);
    char* ptr = m_data.data();
    int size = m_data.size();
    MPI_Irecv(ptr, size, MPI_CHAR, m_proc_id, a_info ...);
}

//! buffer_for_spheres

template<int DIM>
template<class RSA_DATA_STORAGE>
inline void buffer_for_spheres<DIM>::pack(const RSA_DATA_STORAGE& a_rsa_data_storage) {
    buffer_auxi::buffer_for_spheres_<DIM>::pack(0, a_rsa_data_storage.size(),
        a_rsa_data_storage.get_priority_ptr(),
        a_rsa_data_storage.get_phase_ptr(),
        a_rsa_data_storage.get_rad_ptr(),
        a_rsa_data_storage.get_position_ptr());
}

template<int DIM>
template<class RSA_DATA_STORAGE>
inline void buffer_for_spheres<DIM>::unpack(RSA_DATA_STORAGE& rsa_data_storage) {
    uint64_t shift = rsa_data_storage.size();
    uint64_t size_unpack = this->size();
    rsa_data_storage.resize(shift + size_unpack);
    buffer_auxi::buffer_for_spheres_<DIM>::unpack(0, size_unpack,
        rsa_data_storage.get_priority_ptr() + shift,
        rsa_data_storage.get_phase_ptr() + shift,
        rsa_data_storage.get_rad_ptr() + shift,
        rsa_data_storage.get_position_ptr() + shift);
}
