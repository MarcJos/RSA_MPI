//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
/**
 * \file rsa_data_storage.hxx
 * \brief particle storage in a SOA data structure.
 * \author Raphael Prat
 * \version 0.1
 * \date 22 avril 2023
 */


#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>

#include <helper.hxx>
#include <rsa_decoration.hxx>
#include <basic_types.hxx>
#include <rsa_buffer.hxx>

 /**
  *  rsa_data_storage is class to manage particle data.
  *  DIM is the position dimension of the particle.
  */

template<int DIM>
class rsa_data_storage {
	friend class buffer_for_spheres<DIM>;

private:
	vec_int m_priority;  ///< priority
	vec_int m_phase; ///<  global identifier
	vec_double m_rad; ///< sphere radius
	arr_vec_double<DIM> m_pos; ///< position in DIM-D

	// Accessors
public:
	size_t size() const {
		size_t ret = m_rad.size();
		assert(ret == m_pos.size() && "error size rx");
		assert(ret == m_priority.size() && "error size id");
		assert(ret == m_phase.size() && "error size id");
		return ret;
	}

	const vec_d<DIM>& get_center(int a_idx) const { return m_pos[a_idx]; }
	double get_center(int a_idx, int d) const { return get_center(a_idx)[d]; }
	double get_rad(const int a_idx) const { return m_rad[a_idx]; }
	uint64_t get_priority(const int a_idx) const { return m_priority[a_idx]; }
	uint64_t get_phase(const int a_idx) const { return m_phase[a_idx]; }

	//! @return : create an object sphere from the SoA
	//! @param a_index : index of the sphere
	rsa_sphere<DIM> create_sphere(const size_t a_index) const;

	// mutators
public:
	//! @brief add a new sphere
	template<typename... Args>
	void add_sphere(const uint64_t a_priority, const uint64_t a_id, const double a_rad, const Args... a_args);
	//! @brief add a new sphere
	void add_sphere(const uint64_t a_priority, const uint64_t a_id, const double a_rad, const vec_d<DIM>& a_center);
	//! @brief add a new sphere
	void add_sphere(const uint64_t a_priority, const rsa_sphere<DIM>& a_sphere);
	//! @brief add a new sphere
	//! the center is located at a_center[.][a_index_center]
	void add_sphere(const uint64_t a_priority, const uint64_t a_id, const double a_rad,
		int  a_index_center, const arr_vec_double<DIM>& a_centers);
	//! @brief add a new sphere
	//! the center is located at a_center[.][a_index_center]
	void add_sphere(const int a_idx, const rsa_data_storage<DIM>& a_data);
	//! @brief add a new spheres
	template<class PriorityVector>
	void add_spheres(const PriorityVector& a_priority, const vec_int& a_phases, const std::vector<double>& a_rad,
		const arr_vec_double<DIM>& a_centers);
	//! @brief merge another rsa_data_storage into the current one
	void add_spheres(const rsa_data_storage<DIM>& a_in);
	//! @brief : removes the sphere situated at index a_index
	void remove_sphere(const size_t a_index);
	//! @brief mutator
	//! @param a_idx : index
	//! @param a_in : input for phase
	void set_phase(const int a_idx, const uint64_t a_in) { m_phase[a_idx] = a_in; }
	//! @brief mutator
	//! @param a_idx : index
	//! @param a_in : input radius
	void set_rad(const int a_idx, const double a_in) { return m_rad[a_idx] = a_in; }
	//! @brief mutator
	//! @param a_dim : dimension
	//! @param a_value : input value of the shift
	void shift_centers(const double a_dim, const double a_value);

	//! @brief  : This function clear all data in the data storage
	void purge();


	// IO functions
public:
	void print_spheres() const;

	void print_particle(const int a_idx) const;

	void write_spheres(std::string a_name = "spheres_serial.vtk") const;

	void print_log() const;

	// Check functions
	//! @param i, j : indices of particles
	//! @return whether particle i and j intersect each other
	bool check_collision(int i, int j) const;
	//! @return whether no placed particle intersect another one
	bool check_particles() const;
	//! @return whether a given sphere a_sphere intersects an already placed particle
	//! @param a_sphere : a given sphere
	bool check_intersection(const rsa_sphere<DIM>& a_sphere) const;
	//! @return whether a given sphere a_sphere intersects an already placed particle
	//! @param a_sphere : a given sphere
	bool check_radius() const;

protected:
	//! @brief resize all the storages with the new size a_size
	void resize(const int a_size);

private:
	// ptr accessor
	const vec_d<DIM>* get_position_ptr() const { return m_pos.data(); }
	const double* get_rad_ptr() const { return m_rad.data(); }
	const uint64_t* get_priority_ptr() const { return m_priority.data(); }
	const uint64_t* get_phase_ptr() const { return m_phase.data(); }

	// non-const ptr accessors
	vec_d<DIM>* get_position_ptr() { return m_pos.data(); }
	double* get_rad_ptr() { return m_rad.data(); }
	uint64_t* get_priority_ptr() { return m_priority.data(); }
	uint64_t* get_phase_ptr() { return m_phase.data(); }
};

//! @brief qualify the status of spheres
enum class TypeTag {
	Validated, // definitely validated
	Undecided, // can become Validated, Conflict, or Skip
	Conflict, // in conflict with another sphere
	Skip,	// definitely throwed away
	// others are combinations
	Any, // any sphere
	Skip_Or_Undecided, //
	Not_Validated //
};


//! warning : this heritage is highly fragile.
template<int DIM>
class rsa_data_storage_with_tags : public rsa_data_storage<DIM> {
private:
	//! vector of labels (do skip) per cell use in detection
	std::vector<TypeTag> m_tags;

public:
	//! @brief : removes the sphere situated at index a_index
	void remove_sphere(const size_t a_index);
	//! @brief add a new sphere
	//! the center is located at a_center[.][a_index_center]
	void add_sphere(const int a_idx, const rsa_data_storage_with_tags<DIM>& a_data);
	void add_sphere(const int a_idx, const rsa_data_storage<DIM>& a_data, TypeTag typeTag);
	//! @brief : remove all particles tagged by do remove
	template<TypeTag typeTag = TypeTag::Any>
	void remove_tagged_spheres();

	//! @brief does the i-th sphere satisfy the given typeTag?
	template<TypeTag typeTag>
	bool is_tag(int a_idx) const;
	TypeTag get_tag(int a_idx) const { return m_tags[a_idx]; }
	template<TypeTag typeTag>
	void set_tag(int a_idx) { m_tags[a_idx] = typeTag; }

private:
	void remove_tag(const size_t a_index);
};

template<int DIM, class my_class>
const vec_d<DIM>& get_center(int i, const my_class& a_rsa_data_storage) {
	return a_rsa_data_storage.get_center(i);
}



#include <rsa_data_storage.ixx>
