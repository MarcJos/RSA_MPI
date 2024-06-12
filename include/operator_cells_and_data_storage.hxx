//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <basic_types.hxx>
#include <rsa_data_storage.hxx>
#include <helper.hxx>
#include <rsa_random.hxx>


//! @return wWether the 2 spheres collide
//! @param a_rad1, a_rad2 : radius of sphere 1, sphere 2
//! @param a_pos1, a_pos2 : centers of sphere 1, sphere 2
template<int DIM>
bool is_there_collision(double a_rad1, const vec_d<DIM>& a_pos1,
	double a_rad2, const vec_d<DIM>& a_pos2);

//! @brief This function determines if a sphere, defined by its id, radius, and position, collides with at least one sphere within a grid. This function can also take into account priority and whether the sphere has been tagged.
//! @return Whether there is a collision or not. If there is a collision, this function returns the cell id and the position of the particle in contact. Otherwise, it returns false, -1, -1.
//! @param a_id The id of the sphere, with the value corresponding to the priority.
//! @param a_rad The radius of the sphere.
//! @param a_pos The center of the sphere.
//! @param a_cell_id The id of the cell in the grid a_grid.
//! @param a_list_of_cells The list of offset cells for the grid a_grid.
//! @param a_grid The grid of cells containing spheres used for collision detection.
template<int DIM, bool PRIORITY, bool LABEL, bool ONLY_VALIDATED, class VEC>
std::tuple<bool, int, int> detection(uint64_t a_id, double a_rad, const vec_d<DIM>& a_pos, const int a_cell_id, const VEC& a_list_of_cells, const rsa_grid<DIM>& a_grid);

//! @brief This function removes collisions between new spheres and old spheres within a grid.
template<int DIM, bool fast>
void remove_collisions_with_validated_spheres(rsa_grid<DIM>& a_old_cells, std::vector<uint64_t>& a_cells_with_conflicts);

//! @brief Finds and resolves conflicts (if possible) within a grid of cells.
//! @param a_grid : The grid of cells to analyze for conflicts. Output : All conflicts resolved and unresolved are removed
//! @return A data storage containing the unresolved conflicts.
template<int DIM>
void find_solve_conflicts(rsa_grid<DIM>& a_grid, std::vector<uint64_t>& a_cells_with_conflicts);

//! @brief Add spheres in a of grid of cells..
//! @param a_cells : The grid of cells that receives spheres.
//! @param a_data_storage : The data storage containing the spheres.
template<int DIM>
void build_grid(rsa_grid<DIM>& a_cells, const rsa_data_storage_with_tags<DIM>& a_data_storage);
template<int DIM>
void build_grid(rsa_grid<DIM>& a_cells, const rsa_data_storage<DIM>& a_data_storage, TypeTag typeTag);

//! @brief Removes ghost spheres from the grid.
//! @param a_grid : The grid of cells from which ghost particles will be removed.
template<int DIM>
void remove_ghost(rsa_grid<DIM>& a_grid);

//! @brief Copies data from one data storage to another one.
//! @param a_from : The source data storage containing particles to be copied.
//! @param a_dest : The destination data storage to store the copied particles.
template<int DIM>
void copy_particles(const rsa_data_storage<DIM>& a_from, rsa_data_storage<DIM>& a_dest);

//! @brief Copies particles from one grid to another one.
//! @param a_from : The source grid containing particles to be copied.
//! @param a_to : The destination grid to store the copied particles.
template<int DIM>
void copy_particles(const rsa_grid<DIM>& a_from, rsa_grid<DIM>& a_to);

#include<operator_cells_and_data_storage.ixx>
