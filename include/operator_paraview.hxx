//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

namespace rsa_paraview {
//! @brief Exports data from the provided data storage to a Paraview file for debugging and visualization purposes.
//! @param a_name : The name of the output Paraview file.
//! @param a_data : The data storage containing the information to be exported.
template<int DIM>
void debug_paraview(std::string a_name, rsa_data_storage<DIM>& a_data);

//! @brief Exports data from a grid of cells to a Paraview file for debugging and visualization purposes.
//! @param a_name : The name of the output Paraview file.
//! @param a_grid : The grid of cells containing the information to be exported.
template<int DIM, TypeCell>
void debug_paraview(std::string a_name, rsa_grid<DIM>& a_grid);

//! @brief Exports grid data, a draw grid, and conflict data to a Paraview file for debugging and visualization purposes.
//! @param a_name : The name of the output Paraview file.
//! @param a_grid : The grid of cells containing the information to be exported.
//! @param a_draw_grid : The grid used for drawing.
//! @param a_conflicts : The data storage containing conflict information to be exported.
template<int DIM>
void debug_paraview(std::string a_name, rsa_grid<DIM>& a_grid, rsa_grid<DIM>& a_draw_grid, rsa_data_storage<DIM>& a_conflicts);

//! @brief Exports grid data to Paraview files visualization purposes.
//! @param a_name : The name of the output Paraview file.
//! @param a_grid : The grid of cells containing the information to be exported.
template<int DIM>
void paraview(std::string a_name, const rsa_grid<DIM>& a_grid);

//! @brief Exports grid data to Paraview files visualization purposes.
//! @param a_domain : The domain containing the information to be exported.
template<int DIM>
void paraview(const rsa_domain<DIM>& a_domain, std::string a_name = "ParaviewOutput");
} // namespace rsa_paraview

#include <operator_paraview.ixx>
