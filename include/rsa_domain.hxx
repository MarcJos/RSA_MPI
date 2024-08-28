//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <mpi.h>

#include <basic_types.hxx>
#include <rsa_buffer.hxx>
#include <rsa_grid.hxx>
#include <rsa_data_storage.hxx>
#include <rsa_ghost_area.hxx>
#include <operator_cells_and_data_storage.hxx>
#include <operator_ghost.hxx>
//#include <operator_algorithm.hxx>

using namespace sac_de_billes;

template<int DIM>
class rsa_domain {
	//! @brief: interface class for ::rsa_domain<DIM>. Acts as a reference.
	//! Pave a cuboid $`[l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]]`$."
	using Buffer = buffer_for_spheres<DIM>;

private:
	int m_id;  ///< domain id
	vec_d<DIM> m_inf; ///< lower boundary
	vec_d<DIM> m_sup; ///< upper boundary
	vec_d<DIM> m_global_inf; ///< global lower boundary
	vec_d<DIM> m_global_sup; ///< global upper boundary
	rsa_grid<DIM> m_grid; ///< particles data storage
	std::vector<rsa_ghost_area<DIM>> m_ghost; ///< list of ghost areas
	int m_ghost_layer; ///< number of ghost layers
	std::vector<Buffer> m_send; ///< list of buffers used to send particles
	std::vector<Buffer> m_recv; ///< list of buffers used to receive particles
	double m_rad; ///< stores the a_rad given in constructor

	//! @brief : default constructor is private, do not use it.
	rsa_domain() {}
	//! @brief check if the domain is sufficiently large to have cells large enough
	//! @param a_rad : maximal radius of a sphere
	void check_sufficiently_large(double a_rad);

public:

	//! @brief constructor for a periodic domain [l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]] in \R^D (D=DIM)
	//! The domain is scattered over the mpi processes.
	//! @param global_inf : vector  {lmin[0], ..., l_min[D-1]}
	//! @param global_sup : vector {lmax[0], ..., l_max[D-1]}
	//! @param a_ghost_layer : set the number of ghost layers.
	//! @param a_rad : the size of each cell in any direction should be larger than a_rad. This radius is usually the largest radius of placed spheres.
	//! @brief: interface class for ::rsa_domain<DIM>. Acts as a reference.
	//! Pave a cuboid $`[l_min[0], l_max[0]] x ... [l_min[D-1], l_max[D-1]]`$."
	rsa_domain(const vec_d<DIM>& global_inf, const vec_d<DIM>& global_sup,
		const int a_ghost_layer = 0, const double a_rad = 0.05);

	//! @brief: accessor
	double get_m_rad() const { return m_rad; }
	//! @brief : Accessor used to get the data storage (getter)
	rsa_grid<DIM>& get_grid() { return m_grid; }
	const rsa_grid<DIM>& get_grid() const { return m_grid; }
	//! @brief Accesor used to get the inf boundary
	const vec_d<DIM>& get_inf() const { return m_inf; }
	//! @brief Accesor used to get the global inf boundary
	const vec_d<DIM>& get_global_inf() const { return m_global_inf; }
	//! @brief Accesor used to get the sup boundary
	const vec_d<DIM>& get_sup() const { return m_sup; }
	//! @brief Accesor used to get the number of host_layer
	int get_ghost_layer() const { return m_ghost_layer; }
	//! @brief : Accesor used to get the different ghost areas
	std::vector<rsa_ghost_area<DIM>>& get_ghost_areas() { return m_ghost; }
	//! @brief Accesor used to get the send buffers
	std::vector<Buffer>& get_send_buffers() { return m_send; }
	//! @brief Accesor used to get the receive buffers
	std::vector<Buffer>& get_recv_buffers() { return m_recv; }
	//! @brief : const getter
	int get_m_id() const { return m_id; }

	//! @brief : This function computes lower boundaries and upper boundaries
	//! for every subdomains according to their mpi ranks.
	//! @warning : Modify it for having voxels instead of slices.
	void lb(const int a_id, const int a_n_mpi,
		const vec_d<DIM>& global_inf, const vec_d<DIM>& global_sup, const double a_rad);
	//! @brief : This function updates the ghost area with new particles between two draws.
	void update_ghost(rsa_data_storage<DIM>& a_ghost_data);
	//!
	void add(rsa_data_storage<DIM>& a_data) { build_grid(m_grid, a_data, TypeTag::Validated); }
	//! @brief : Print some general informations about the domain.
	void domain_log() const;

	//! @brief: get the total volume of the domain
	double get_total_volume() const;
	//! @return: the total volume of spheres contained in all the mpi domains
	double compute_total_volume_of_spheres() const;
	//! @return list of real spheres inside the domain
	vector<rsa_sphere<DIM>> extract_spheres() const { return this->get_grid().extract_spheres(); }
};

namespace auxi_rsa_domain {
template<int DIM>
void create_layer_bounds_inplace(const vec_d<DIM>& domain_inf, const vec_d<DIM>& domain_sup,
	sac_de_billes::DiscPoint<DIM> relative_position, vec_d<DIM>& layer_inf, vec_d<DIM>& layer_sup, double a_rad);
} // namespace  auxi_rsa_domain

#include<rsa_domain.ixx>
