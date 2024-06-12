//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <cstring>
#include <rsa_grid.hxx>
#include <rsa_ghost_area.hxx>
#include <operator_cells_and_data_storage.hxx>
#include <helper.hxx>

template<int DIM>
using Buffers = std::vector<buffer_for_spheres<DIM>>;
template<int DIM>
using GhostAreas = std::vector<rsa_ghost_area<DIM>>;

//! @brief : Fill the different ghost areas with a list of particles (a_position)
template<int DIM, TypeTag typeTag = TypeTag::Any>
void fill_ghost(const rsa_grid<DIM>& a_grid, std::vector<rsa_ghost_area<DIM>>& a_ghost);

//! @brief : Fill the buffer with informations in a buffer.
template<int DIM>
void fill_buffer_with_ghost(rsa_ghost_area<DIM>& a_ghost, buffer_for_spheres<DIM>& a_buffer);

//! @brief Updates ghost particles by exchanging data between processes.
//! @param a_recv : The vector of buffer_for_spheres storing particles received from neighboring processes.
//! @param a_send : The vector of buffer_for_spheres storing particles to be sent to neighboring processes.
//! @param a_ghost : The vector of rsa_ghost_area containing ghost particle information.
//! This vector is consumed after the process.
//! @param a_data : The data storage to be updated with exchanged ghost particles.
template<int DIM>
void update_ghost_impl(
  std::vector<buffer_for_spheres<DIM>>& a_recv,
  std::vector<buffer_for_spheres<DIM>>& a_send,
  std::vector<rsa_ghost_area<DIM>>& a_ghost,
  rsa_data_storage<DIM>& a_data);

//! @brief Upate the a_grid with ghosts from a_ghost_areas, accross MPI
//! @param a_grid : a given rsa_grid
//! @param a_recv, a_send : buffer for receiving and send
//! @param a_ghost_areas : ghost areas that should be transferred accross MPI
template<int DIM>
void put_ghost_into(rsa_grid<DIM>& a_grid, GhostAreas<DIM>& a_ghost_areas,
  Buffers<DIM>& a_recv, Buffers<DIM>& a_send, TypeTag typeTag);

#include<operator_ghost.ixx>
