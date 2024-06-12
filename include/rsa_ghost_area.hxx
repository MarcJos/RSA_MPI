//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <helper.hxx>
#include <rsa_data_storage.hxx>

//! @brief class for storing spheres in zones that should be exchanged.
//! Each domain possesses ghost areas around it.
template<int DIM>
class rsa_ghost_area {
private:
	static constexpr rsa_helper::helper_check_interval<DIM> check_interval = {};

	// members
	bool m_periodic; ///< if this section is
	vec_d<DIM> m_inf; ///< lower boundary
	vec_d<DIM> m_sup; ///< upper boundary
	rsa_data_storage<DIM> m_data; ///< data storage
	vec_d<DIM> m_shift;

public:
	//! @brief constructor
	rsa_ghost_area(const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup, const bool a_periodic, const vec_d<DIM>& a_shift) :
		m_periodic(a_periodic), m_inf(a_inf), m_sup(a_sup), m_data(), m_shift(a_shift) {}

	//! @brief  : This function returns the number of particles in the ghost area
	int size() { return m_data.size(); }
	//! @brief  : This function gets the data storage (m_data)
	rsa_data_storage<DIM>& get_data() { return m_data; }
	const rsa_data_storage<DIM>& get_data() const { return m_data; }
	//! @brief Accesor used to get the inf boundary
	const vec_d<DIM>& get_inf() const { return m_inf; }
	//! @brief Accesor used to get the sup boundary
	const vec_d<DIM>& get_sup() const { return m_sup; }

	//! @param ...a_args : argument necessary for adding a particle
	//! @return : whether the particle has really been added or not
	template<typename... Args>
	bool add_particle(Args&&... a_args);
	//! @return : whether a particle is a ghost or not
	//! @param ...a_args : arguments for the center of the particle
	template<typename... Args>
	bool is_ghost(Args&&... a_args);



	//! @brief  : This function purges the data storage between two draws
	void purge() { m_data.purge(); }
	//! @brief  : This function applies a shift value (m_shift)
	void apply_periodicity();
};

#include<rsa_ghost_area.ixx>
