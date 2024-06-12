//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

#include <rsa_grid.hxx>
#include <rsa_data_storage.hxx>

template<int DIM>
bool check_no_sphere_in_ghost(const rsa_grid<DIM>& a_grid) {
	auto& ghosts = a_grid.get_ghost_traversal();
	for (int cell = 0; cell < ghosts.size(); cell++) {
		const auto& spheres = a_grid.get_data(ghosts[cell]);
		if (spheres.size() != 0) {
			return false;
		}
	}
	return true;
}

template<int DIM>
bool check_particles_in_cells(const rsa_grid<DIM>& a_grid) {
	rsa_helper::helper_extractor<DIM> extract;
	vec_d<DIM> position;
	for (int cell = 0; cell < a_grid.size(); cell++) {
		const auto& spheres = a_grid.get_data(cell);
		for (int s1 = 0; s1 < spheres.size(); s1++) {
			extract(spheres.get_position(), s1, position);
			const int cell_id = a_grid.compute_cell_idx(position);
			if (cell_id != cell) return false;
		}
	}
	return true;
}

template<int DIM>
bool check_no_doublon(const rsa_grid<DIM>& a_grid) {
	int number_of_doublons = 0;
	rsa_helper::helper_extractor<DIM> extract;
	vec_d<DIM> position_1;
	vec_d<DIM> position_2;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int cell = 0; cell < a_grid.size(); cell++) {
		const auto& spheres = a_grid.get_data(cell);
		for (int s1 = 0; s1 < spheres.size(); s1++) {
			for (int s2 = s1 + 1; s2 < spheres.size(); s2++) {
				const bool same_id = spheres.get_priority(s1) == spheres.get_priority(s2);
				position_1 = spheres.get_center(s1);
				position_2 = spheres.get_center(s2);
				bool same_position = true;
				for (int dim = 0; dim < DIM; dim++) {
					//if (position_1[dim] != position_2[dim]) {
					if (std::fabs(position_1[dim] - position_2[dim]) >= 1e-16) {
						same_position = false;
						break;
					}
				}
				if (same_position) {
					std::cerr << "Rank[" << rank << "] doublon[" << number_of_doublons << "] "
						<< " index1 " << s1
						<< " index2 " << s2
						<< " id1 " << spheres.get_priority(s1)
						<< " id2 " << spheres.get_priority(s2)
						<< " posx " << position_1[0]
						<< " posy " << position_1[1] << std::endl;
					number_of_doublons++;
				}
			}
		}
	}
	if (number_of_doublons > 0) {
		std::cerr << "Rank[" << rank << "] contains " << number_of_doublons << std::endl;
		return false;
	}
	return true;
}

template<int DIM>
void remove_doublons(rsa_grid<DIM>& a_grid) {
	rsa_helper::helper_extractor<DIM> extract;
	vec_d<DIM> position_1;
	vec_d<DIM> position_2;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int compt = 0;
	for (int cell = 0; cell < a_grid.size(); cell++) {
		auto& spheres = a_grid.get_data(cell);
		for (int s1 = 0; s1 < spheres.size(); s1++) {
			for (int s2 = spheres.size() - 1; s2 > s1; s2--) {
				const bool same_id = spheres.get_priority(s1) == spheres.get_priority(s2);
				extract(spheres.get_position(), s1, position_1);
				extract(spheres.get_position(), s2, position_2);
				bool same_position = true;
				for (int dim = 0; dim < DIM; dim++) {
					if (std::fabs(position_1[dim] - position_2[dim]) >= 1e-16) {
						same_position = false;
						break;
					}
				}
				if (same_position) {
					compt++;
					spheres.remove_sphere(s2);
				}
			}
		}
	}
	if (compt > 0) {
		std::cerr << "Rank[" << rank << "] remove " << compt << " doublon(s)" << std::endl;
	}
}
