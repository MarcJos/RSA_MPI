//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <basic_types.hxx>
#include <rsa_data_storage.hxx>
#include <rsa_grid_traversal.hxx>

//! @brief : rsa_grid stores the list of cells
//! pave the domain delimitated by m_shift and m_shift + m_cell_length.
//! this paving strictly encloses the domain delimitated by a_inf, a_sup
//! restricting this paving to the inner cells exactly paves the domain delimitated by a_inf, a_sup
template<int DIM>
class rsa_grid : public rsa_grid_traversal<DIM> {
private:
	std::vector<rsa_data_storage_with_tags<DIM>> m_data; ///< list of data per cell
	uint64_t number_of_spheres_fast;

public:
	//! @brief : default constructor
	rsa_grid() : rsa_grid_traversal<DIM>(), number_of_spheres_fast{ 0 } {}
	//! @brief : this constructor calls the function "define".
	//! @see : define.
	rsa_grid(const double a_rad, const int a_ghost_layer, const vec_d<DIM>& a_inf, const vec_d<DIM>& a_sup)
		: rsa_grid_traversal<DIM>(a_rad, a_ghost_layer, a_inf, a_sup), number_of_spheres_fast{ 0 } {
		m_data.resize(this->size());
	};

	//! extract the data into a single rsa_data_storage
	template<TypeCell typeCell, TypeTag typeTag = TypeTag::Any>
	rsa_data_storage<DIM> extract_data() const;

	//! @brief : remove all particles tagged by do remove
	template<TypeCell typeCell, TypeTag typeTag>
	void remove_tagged_spheres();
	template<TypeCell typeCell, TypeTag typeTag>
	void remove_tagged_spheres(std::vector<uint64_t>& restrict_to_cells_id);
	//! @brief : purge cells
	template<TypeCell typeCell = TypeCell::All>
	void purge();
	//! @brief copy the particles of a_from inside it
	//! @param
	template<TypeCell typeCell>
	void add_spheres(const rsa_grid<DIM>& a_from);
	//! @brief add a single sphere to the grid
	//! @param i : index of the sphere
	//! @param a_spheres : position of the sphere inside the storage a_spheres
	void add_single_sphere(int i, const rsa_data_storage<DIM>& a_spheres, TypeTag typeTag);

	//! @brief : This function sets to true the a_idx-th sphere of the cell a_cell
	template<TypeTag typeTag>
	void set_tag(int a_cell, const int a_idx) { m_data[a_cell].template set_tag<typeTag>(a_idx); }
	//! @brief : return true is the a_idx-th sphere of the cell a_cell is tagged/labeled
	template<TypeTag typeTag>
	bool is_tag(int a_cell, const int a_idx) const { return m_data[a_cell].template is_tag<typeTag>(a_idx); }

	//! @brief : get the number of spheres
	//! @warning : too expensive !
	template<TypeCell typeCell = TypeCell::Real>
	uint64_t get_number_of_spheres() const;
	//! @brief : setter
	void set_number_of_spheres_fast(uint64_t number_of_spheres_fast_) { number_of_spheres_fast = number_of_spheres_fast_; }
	//! @brief : const getter
	uint64_t get_number_of_spheres_fast() const;
	//! @brief : total volume of spheres *locally* on the process
	double local_volume_of_spheres() const;
	//! @return list of real spheres inside the domain
	std::vector<rsa_sphere<DIM>> extract_spheres() const;

	//! @brief : get the data storage of the cell a_idx
	const rsa_data_storage_with_tags<DIM>& get_data(const int a_idx) const { return m_data[a_idx]; }
	rsa_data_storage_with_tags<DIM>& get_data(const int a_idx) { return m_data[a_idx]; }

	//! @return all the spheres that are tagged as "conflict"
	template<TypeCell typeCell = TypeCell::All>
	rsa_data_storage<DIM> get_conflict_spheres() const;

	//! @brief check whether particles are intersected
	//! @warning : do not use this function outside of a assert function (memory footprint, slow)
	bool check_particles() const;
	//! @return : how many data of given tag
	template<TypeCell typeCell, TypeTag typeTag>
	int64_t how_many(std::vector<uint64_t>& cells_id) const;
	template<TypeCell typeCell, TypeTag typeTag>
	int64_t how_many() const;

	//! @brief  apply the Lambda for each cell index of desired type
	//! @tparam function : function depending on the index, to be applied
	//! @tparam typeCell : TypeCell for traversal
	template<TypeCell typeCell, bool use_omp, class LambdaFunction>
	void apply_on_cells(LambdaFunction& function) const;
	template<TypeCell typeCell, bool use_omp, class LambdaFunction>
	void apply_on_cells(LambdaFunction& function);
	//
	template<TypeCell typeCell, bool use_omp, class LambdaFunction>
	void apply_on_cells(std::vector<uint64_t>& cells_id, LambdaFunction& function) const;
	template<TypeCell typeCell, bool use_omp, class LambdaFunction>
	void apply_on_cells(std::vector<uint64_t>& cells_id, LambdaFunction& function);
};

#include<rsa_grid.ixx>
