//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <rsa_grid.hxx>
#include <rsa_data_storage.hxx>

//! @brief Checks if there are no spheres in the ghost areas of the grid of cells.
//! @param a_grid : The grid of cells to check for spheres in ghost areas.
template<int DIM>
bool check_no_sphere_in_ghost(const rsa_grid<DIM>& a_grid);

//! @brief Checks if all particles are present and accounted for in the cells of the grid.
//! @param a_grid : The grid of cells to check for particles.
template<int DIM>
bool check_particles_in_cells(const rsa_grid<DIM>& a_grid);

//! @brief Checks if there are no duplicate spheres in the cells of the grid.
//! @param a_grid : The grid of cells to check for duplicate spheres.
template<int DIM>
bool check_no_doublon(const rsa_grid<DIM>& a_grid);

//! @brief Removes any duplicate spheres from the cells of the grid.
//! @param a_grid : The grid of cells from which duplicate spheres will be removed.
template<int DIM>
void remove_doublons(rsa_grid<DIM>& a_grid);

#include <operator_check_cell.ixx>
