//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM>
bool is_there_collision(double a_rad1, const vec_d<DIM>& a_pos1,
	double a_rad2, const vec_d<DIM>& a_pos2) {
	double distanceCarre = 0;
	for (int d = 0; d < DIM; d++) {
		const double dr = a_pos1[d] - a_pos2[d];
		distanceCarre += dr * dr;
	}
	return distanceCarre < sac_de_billes::auxi_function::puissance<2>(a_rad1 + a_rad2);
}

template<int DIM, bool PRIORITY, bool SKIP, bool ONLY_VALIDATED, class VEC>
std::tuple<bool, int, int> detection(uint64_t a_id, double a_rad, const vec_d<DIM>& a_pos,
	const int a_cell_id, const VEC& a_list_of_cells, const rsa_grid<DIM>& a_grid) {
	for (size_t ncell = 0; ncell < a_list_of_cells.size(); ncell++) {
		const int cell_id = a_cell_id + a_list_of_cells[ncell];
		const auto& spheres = a_grid.get_data(cell_id);
		for (size_t s = 0; s < spheres.size(); s++) {
			if constexpr (ONLY_VALIDATED) {
				if (not a_grid.template is_tag<TypeTag::Validated>(cell_id, s)) {
					continue;
				}
			}
			if constexpr (SKIP) // if the sphere is tagged as skipped
			{
				if (a_grid.template is_tag<TypeTag::Skip>(cell_id, s)) {
					continue;
				}
			}

			const double rad = spheres.get_rad(s);
			const auto& pos = spheres.get_center(s);

			if constexpr (PRIORITY) { // if priority between spheres is to be considered
				if (spheres.get_priority(s) < a_id) {
					if (is_there_collision<DIM>(a_rad, a_pos, rad, pos)) {
						return { true, cell_id, s };
					}
				}
			} else {
				if (is_there_collision<DIM>(a_rad, a_pos, rad, pos)) {
					return { true, cell_id, s };
				}
			}
		}
	}
	return { false, -1, -1 };
}

template<int DIM>
void remove_collisions_with_validated_spheres(rsa_grid<DIM>& a_old_grid,
	std::vector<uint64_t>& a_cells_with_conflicts) {
	constexpr bool use_priority = false;
	constexpr bool use_skip = false;
	constexpr bool use_only_validated = true;
	auto& nc_list = a_old_grid.get_list_of_neighbor_cells();
	// iterate over cells of new spheres
	auto remove_collision_loc = [&a_old_grid, &nc_list, use_priority,
		use_skip, use_only_validated](int cell_id, auto& new_spheres) {
		for (int s = new_spheres.size() - 1; s >= 0; s--) {
			assert((new_spheres.template is_tag<TypeTag::Undecided>(s) or new_spheres.template is_tag<TypeTag::Validated>(s)));
			if (new_spheres.template is_tag<TypeTag::Undecided>(s)) {
				const double rad = new_spheres.get_rad(s);
				const auto& pos = new_spheres.get_center(s);
				auto [boom, n, i] = detection<DIM, use_priority, use_skip, use_only_validated>(-1, rad, pos, cell_id, nc_list, a_old_grid);
				if (boom) new_spheres. template set_tag<TypeTag::Skip>(s);
			}
		}
		};
	a_old_grid.template apply_on_cells<TypeCell::Real, true>(a_cells_with_conflicts, remove_collision_loc);
	a_old_grid.template remove_tagged_spheres<TypeCell::Real, TypeTag::Skip>(a_cells_with_conflicts);
}


template<int DIM>
void find_solve_conflicts(rsa_grid<DIM>& a_grid,
	std::vector<uint64_t>& a_cells_with_conflicts) {

	using info = std::pair<int, int>;
	constexpr bool use_priority(true), use_skip(true);

	std::vector<info> ind = {};
	const auto& neigh_cells_list = a_grid.get_list_of_neighbor_cells();
	for (uint64_t id = 0; id < a_cells_with_conflicts.size(); id++) {
		uint64_t cell = a_cells_with_conflicts[id];
		info add;
		add.first = cell;
		auto& spheres = a_grid.get_data(cell);
		for (size_t s = 0; s < spheres.size(); s++) {
			if (a_grid.template is_tag<TypeTag::Undecided>(cell, s)) {
				add.second = s;
				ind.push_back(add);
			}
		}
	}

	sort(ind.begin(), ind.end(),
		[&a_grid](const std::pair<int, int>& a, const std::pair<int, int>& b) {
			int i0 = a_grid.get_data(a.first).get_priority(a.second);
			int i1 = a_grid.get_data(b.first).get_priority(b.second);
			if (i0 != i1) {
				return i0 < i1;
			} else {
				std::uniform_real_distribution<> distrib(-1., 1.);
				std::mt19937 random_generator{};
				double d_pos_00 = a_grid.get_data(a.first).get_center(a.second, 0);
				double d_pos_11 = a_grid.get_data(b.first).get_center(b.second, 0);
				d_pos_00 /= (d_pos_00 + d_pos_11) * 0.5 * std::numeric_limits<int>::max();
				d_pos_11 /= (d_pos_00 + d_pos_11) * 0.5 * std::numeric_limits<int>::max();
				int pos_00 = static_cast<int>(d_pos_00);
				int pos_10 = static_cast<int>(d_pos_11);
				std::minmax(pos_00, pos_10);
				std::seed_seq seq{ pos_00 , pos_10 };
				random_generator.seed(seq);
				double criterion = distrib(random_generator);
				while (criterion == 0) {
					criterion = distrib(random_generator);
				}
				return criterion < 0.;
			}
		});

	for (size_t index = 0; index < ind.size(); index++) {
		auto [cell_id, sphere_id] = ind[index];
		if (a_grid.is_ghost(cell_id)) continue;
		const auto& spheres = a_grid.get_data(cell_id);
		const double rad = spheres.get_rad(sphere_id);
		const uint64_t id = spheres.get_priority(sphere_id);
		const auto& pos = spheres.get_center(sphere_id);
		auto [boom, cell_id_2, sid] = detection<DIM, use_priority, use_skip, false>(id, rad, pos, cell_id, neigh_cells_list, a_grid);
		if (boom) {
			if (a_grid.is_ghost(cell_id_2) or a_grid.template is_tag<TypeTag::Conflict>(cell_id_2, sid)) {
				a_grid.template set_tag<TypeTag::Conflict>(cell_id, sphere_id);
			} else {
				a_grid.template set_tag<TypeTag::Skip>(cell_id, sphere_id);
			}
		}
	}
}

template<int DIM>
void build_grid(rsa_grid<DIM>& a_grid, const rsa_data_storage_with_tags<DIM>& a_data) {
	const int size = a_data.size();
	for (int i = 0; i < size; i++) {
		a_grid.add_single_sphere(i, a_data);
	}
}

template<int DIM>
void build_grid(rsa_grid<DIM>& a_grid, const rsa_data_storage<DIM>& a_data, TypeTag typeTag) {
	const int size = a_data.size();
	for (int i = 0; i < size; i++) {
		a_grid.add_single_sphere(i, a_data, typeTag);
	}
}


template<int DIM>
void remove_ghost(rsa_grid<DIM>& a_grid) {
	a_grid.template purge<TypeCell::Ghost>();
}

template<int DIM>
void copy_particles(const rsa_data_storage<DIM>& a_from, rsa_data_storage<DIM>& a_dest) {
	a_dest.add_spheres(a_from);
}

template<int DIM>
void copy_particles(const rsa_grid<DIM>& a_from, rsa_grid<DIM>& a_to) {
	a_to.add_spheres(a_from);
}

