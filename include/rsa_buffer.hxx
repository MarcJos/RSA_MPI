//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <mpi.h>
#include <vector>

template<typename... Args>
class rsa_buffer {
public:
	//! @brief constructor
	//! @param a_proc_id : id of the processor
	//! @param a_size : size of the data in terms of nb of objects to send
	rsa_buffer(int a_proc_id, int a_size);
	//! @brief default constructor
	rsa_buffer() : rsa_buffer(0, 0) {}

private:
	int m_proc_id;
	int m_size;
	std::vector<char> m_data;

public://! @brief getter
	uint64_t size() const { return m_size; }
	//! @brief getter
	int get_proc_id() const { return m_proc_id; }
	//! @brief getter
	std::vector<char>& get_data() { return m_data; }
	//! @brief getter
	const std::vector<char>& get_data() const { return m_data; }

	//! @brief resize the buffer to contain a_size objects
	//! @param a_size
	void set_size(size_t a_size);
	//! @brief : setter
	void set_proc_id(size_t a_proc_id) { m_proc_id = a_proc_id; }


	//! @brief remove all the objects in the buffer
	void purge();

	//! @brief : pack all the data contained in the pointers m_args into m_data as chars
	//! @param a_first, a_last : m_arg[a_first ... a_last[ should be packed
	//! @param ...m_args : variadic quantities to be packed
	void pack(int a_first, int a_last, const Args*... m_args);
	//! @see : pack(0, a_size, m_args...);
	void pack(int a_size, const Args*... m_args);

	// unpacking step

	//! @brief unpack all the data contained in m_data into the pointers m_args
	//! @param a_first, a_last  : the m_args should be filled at indices [a_first, ..., a_last[
	//! @param ...m_args : variadic pointers to be filled into
	void unpack(int a_first, int a_last, Args*... m_args) const;
	//! @see : unpack(0, a_size, m_args...);
	void unpack(int a_size, Args*... m_args) const;

	//! @brief MJ : warning ne comprends pas
	//! @param ...a_info
	template<typename ... MpiInfo>
	void send_size(MpiInfo&&... a_info);
	//! @brief MJ : warning ne comprends pas
	//! @param ...a_info
	template<typename ... MpiInfo>
	void recv_size(MpiInfo&&... a_info);
	//! @brief MJ : warning ne comprends pas
	//! @param ...a_info
	template<typename ... MpiInfo>
	void send_data(MpiInfo&&... a_info);
	//! @brief MJ : warning ne comprends pas
	//! @param ...a_info
	template<typename ... MpiInfo>
	void recv_data(MpiInfo&&... a_info);

private:
	//! @brief resize the buffer to contain a_size objects
	//! @param a_size
	void resize(size_t a_size);
	//! @return size of a unique object
	constexpr size_t args_size() const { return (sizeof(Args) + ...); }


	template<typename U>
	void pack_impl(char* a_data, int a_first, int a_last, int a_shift, const U* a_ptr);
	template<typename U, typename... miniArgs>
	void pack_impl(char* a_data, int a_first, int a_last, int a_shift, const U* a_ptr, const miniArgs*... a_args);

	template<typename U, typename... miniArgs>
	void unpack_impl(const char* a_data, int a_first, int a_last, int a_shift, U* a_ptr, miniArgs*... a_args) const;
	template<typename U>
	void unpack_impl(const char* a_data, int a_first, int a_last, int a_shift, U* a_ptr) const;
};

namespace buffer_auxi {
template<int DIM>
using buffer_for_spheres_ = rsa_buffer<uint64_t, uint64_t, double, vec_d<DIM>>;
}// namespace  buffer_auxi

//! @brief the last double are here to store the center of the sphere
//! @tparam DIM : dimension of the space
template<int DIM>
class buffer_for_spheres : public buffer_auxi::buffer_for_spheres_<DIM> {
public:
	buffer_for_spheres(int a_proc_id, int a_size) : buffer_auxi::buffer_for_spheres_<DIM>(a_proc_id, a_size) {}
	buffer_for_spheres() : buffer_auxi::buffer_for_spheres_<DIM>() {}

	template<class RSA_DATA_STORAGE>
	void pack(const RSA_DATA_STORAGE& rsa_data_storage);
	template<class RSA_DATA_STORAGE>
	void unpack(RSA_DATA_STORAGE& rsa_data_storage);
private:
	using buffer_auxi::buffer_for_spheres_<DIM>::unpack;
	using buffer_auxi::buffer_for_spheres_<DIM>::pack;
};




#include<rsa_buffer.ixx>
