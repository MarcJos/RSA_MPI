//! Copyright : Apache 2.0, see LICENSE 
//! 
//! Copyright : see license.txt
//!
//! \brief
#pragma once

#include<median_of_medians.hxx>
#include <chrono>
#include <ctime>
#include <cmath>

namespace algorithm {

template<int DIM>
rsa_algo<DIM>::rsa_algo(rsa_domain<DIM>& a_domain, RadiusGenerator<DIM>& a_radius_generator,
	int a_nbshots_singledraw, int a_n_draw) :
	m_nbshots_singledraw{ std::move(a_nbshots_singledraw) },
	m_n_draw{ a_n_draw },
	m_radius_generator{ a_radius_generator },
	m_domain{ a_domain },
	m_ghost_data{},
	m_r_max{ a_radius_generator.get_max_radius() },
	miss_rate{ auxi::magical_default_miss_rate<DIM>() } {
	this->get_grid() = rsa_grid<DIM>(a_radius_generator.get_max_radius(), a_domain.get_ghost_layer(),
		a_domain.get_inf(), a_domain.get_sup());
	if (a_radius_generator.get_max_radius() > a_domain.get_m_rad()) {
		throw runtime_error("Impossible to have the maximal radius of radius generator larger than domain implicit radius");
	}
}


template<int DIM>
void rsa_algo<DIM>::proceed_naive(std::mt19937& random_generator) {
	auto ulaw = auxi::create_random_law(m_domain, random_generator);
	// get data storage
	auto& ghost_areas = m_domain.get_ghost_areas();
	auto center_generator = [&ulaw](int size) {
		return auxi::generate_sphere_positions<DIM>(ulaw.data(), size);
		};
	//!
	auto priority_generator = [&random_generator](int size) {
		return generate_priority<int>(size, random_generator);
		};
	//!
	const auto& radius_gen = this->m_radius_generator;
	auto radius_generator = [&radius_gen, &random_generator](int size) {
		return radius_gen(size, random_generator);
		};

	// This loop iterates over draw
	for (int draw = 0; draw < m_n_draw; draw++) {
		auto nb_shots = m_nbshots_singledraw;
		int64_t total_nb_shots = nb_shots * rsa_mpi::get_number_of_mpi_processes();
		bool may_outreach_nb_spheres = (total_nb_shots >= m_radius_generator.get_current_number());
		uint64_t nb_spheres_total_max = m_radius_generator.get_current_number();
		int64_t nb_new_spheres_loc = single_draw(center_generator, priority_generator, radius_generator,
			nb_shots, nb_spheres_total_max, may_outreach_nb_spheres);
		int64_t nb_new_spheres_glob = rsa_mpi::compute_mpi_sum(nb_new_spheres_loc);
		m_radius_generator.update_placed(nb_new_spheres_glob);
		bool should_continue = m_radius_generator.is_there_still_radii();
		if (not should_continue) {
			rsa_mpi::message("End of generation, due to lack of new radii");
			break;
		}
	}
}



template<int DIM>
void rsa_algo<DIM>::check_vox_time() {
	voxel_list::list_of_voxels<DIM> uncovered_voxels(m_domain.get_inf(), m_domain.get_sup() - m_domain.get_inf(),
		(m_radius_generator.get_max_radius() + m_radius_generator.get_min_radius()));
	for (int draw = 0; true; draw++) {
		int64_t total_nb_vox = rsa_mpi::compute_mpi_sum(static_cast<int64_t>(uncovered_voxels.size()));
		if (total_nb_vox == 0) {
			rsa_mpi::message("End of check");
			break;
		}
		update_covered_voxels(uncovered_voxels, 1.);
	}
}

template<int DIM>
template<bool plog>
void rsa_algo<DIM>::proceed_voxel(std::mt19937& random_generator) {
	// create voxels for generating spheres
	voxel_list::list_of_voxels<DIM> uncovered_voxels(m_domain.get_inf(), m_domain.get_sup() - m_domain.get_inf(),
		(m_r_max + m_radius_generator.get_min_radius()));
	auto center_generator = [&uncovered_voxels, &random_generator](int size) {
		return uncovered_voxels.pick_points(size, random_generator);
		};
	//!
	auto priority_generator = [&random_generator](int size) {
		return generate_priority<int>(size, random_generator);
		};
	//!
	const auto& radius_gen = this->m_radius_generator;
	auto radius_generator = [&radius_gen, &random_generator](int size) {
		return radius_gen(size, random_generator);
		};

	// This loop iterates over draw
	uint64_t old_nb_spheres = 0;

	int64_t final_nb_miss = 0, final_nb_shots = 0;

	constexpr bool timelog = false;
	std::chrono::time_point<std::chrono::system_clock> start_time;
	std::chrono::time_point<std::chrono::system_clock> current_time;
	if constexpr (timelog) {
		start_time = std::chrono::system_clock::now();
		rsa_mpi::message("Time:", 0, "ms");
	}

	for (int draw = 0; true; draw++) {
		//! warning : assumes all the rsa_domains are all equally subdivided with same voxel_lengths
		double intensity = compute_intensity_poisson(uncovered_voxels.size());
		double maximum_intensity = rsa_mpi::compute_mpi_max(intensity);

		uint64_t total_nb_vox = rsa_mpi::compute_mpi_sum(static_cast<int64_t>(uncovered_voxels.size()));
		if (total_nb_vox == 0) {
			rsa_mpi::message("[log rsa_mpi] End of algorithm : fully packed");
			break;
		}

		uint64_t nb_shots = compute_nb_shots_voxel(random_generator, intensity, maximum_intensity);
		uint64_t total_nb_shots = rsa_mpi::compute_mpi_sum(nb_shots);
		if (total_nb_shots == 0) {
			continue;
		}

		bool may_outreach_nb_spheres = (total_nb_shots >= m_radius_generator.get_current_number());


		uint64_t nb_spheres_total_max = m_radius_generator.get_current_number();
		uint64_t nb_new_spheres_loc = single_draw(center_generator, priority_generator, radius_generator,
			nb_shots, nb_spheres_total_max, may_outreach_nb_spheres);
		uint64_t nb_new_spheres_glob = rsa_mpi::compute_mpi_sum(nb_new_spheres_loc);
		m_radius_generator.update_placed(nb_new_spheres_glob);
		bool should_continue = m_radius_generator.is_there_still_radii();

		// measure efficiency
		int new_nb_spheres = this->get_grid().template get_number_of_spheres_fast();
		assert(new_nb_spheres > old_nb_spheres);
		assert(nb_shots > new_nb_spheres - old_nb_spheres);
		uint64_t nb_miss = nb_shots - (new_nb_spheres - old_nb_spheres);
		old_nb_spheres = new_nb_spheres;
		uint64_t total_nb_miss = rsa_mpi::compute_mpi_sum(nb_miss); // could fuse these reduction operations with nb_shots?
		double miss_rate = total_nb_miss / (0.000001 + total_nb_shots);
		if constexpr (plog) { // print output
			rsa_mpi::message("--> Miss rate = " + to_string(int(100 * miss_rate)) + "%"
				+ "               ; Nb shots = " + to_string(total_nb_shots)
				+ " ; Nb miss = ", to_string(total_nb_miss)
			);
		}
		final_nb_miss += total_nb_miss;
		final_nb_shots += total_nb_shots;

		//! update covered voxels
		update_covered_voxels<plog>(uncovered_voxels, miss_rate); // optimize_efficiency

		if (not should_continue) {
			rsa_mpi::message("End of generation, due to lack of new radii");
			break;
		}
		if constexpr (timelog) {
			current_time = std::chrono::system_clock::now();
			const double duration = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
			rsa_mpi::message("Time:", duration, "ms");
		}
	}

	{ // Final message
		if constexpr (timelog) {
			current_time = std::chrono::system_clock::now();
			const double duration = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
			rsa_mpi::message("Time:", duration, "ms");
		}

		double miss_rate = (1. * final_nb_miss) / final_nb_shots;
		rsa_mpi::message("[log rsa_mpi] Final Miss rate = " + to_string(int(100 * miss_rate)) + "%"
			+ "               ; Final Nb shots = " + to_string(final_nb_shots)
			+ " ; Final Nb miss = " + to_string(final_nb_miss)
			+ " ; \n[log rsa_mpi] Total Nb spheres = " + to_string(final_nb_shots - final_nb_miss));
	}
}

template<int DIM>
template<bool plog>
void  rsa_algo<DIM>::update_covered_voxels(voxel_list::list_of_voxels<DIM>& uncovered_voxels,
	double miss_rate) {
	if (miss_rate > desired_miss_rate()) {
		uncovered_voxels.remove_covered(this->get_grid(), this->m_radius_generator.get_min_radius());
		uncovered_voxels.subdivide_uncovered(this->get_grid(), this->m_radius_generator.get_min_radius());
	}

	if constexpr (plog) {
		int64_t nb_voxels = uncovered_voxels.size();
		int64_t max_nb_voxels = rsa_mpi::compute_mpi_max(nb_voxels);
		int64_t total_nb_voxels = rsa_mpi::compute_mpi_sum(nb_voxels);
		rsa_mpi::message("Total nb voxels: " + std::to_string(total_nb_voxels));
		rsa_mpi::message("Max nb voxels  : " + std::to_string(max_nb_voxels));
	}
}

template<int DIM>
int rsa_algo<DIM>::compute_nb_shots_voxel(std::mt19937& random_generator,
	double intensity, double maximum_intensity) {
	double mutliplicator_intensity = m_nbshots_singledraw / maximum_intensity;
	if (mutliplicator_intensity > 1) {
		mutliplicator_intensity = 1.;
	}
	std::poisson_distribution<> poisson_distrib(mutliplicator_intensity * intensity);
	int nb_shots = poisson_distrib(random_generator);
	return nb_shots;
}

template<int DIM>
double rsa_algo<DIM>::compute_intensity_poisson(int a_number_of_voxels) {
	return m_intensity_factor * a_number_of_voxels;
}

template<int DIM>
template<int method, bool plog>
void rsa_algo<DIM>::proceed(size_t seed) {
	static_assert(method == 0 or method == 1);
	std::mt19937 random_generator = law::create_random_generator(seed);
	if constexpr (method == 0) {
		this->proceed_naive(random_generator);
	}
	if constexpr (method == 1) {
		this->proceed_voxel<plog>(random_generator);
	}
}

template<int DIM>
void rsa_algo<DIM>::reset_radius_generator(RadiusGenerator<DIM>& a_radius_generator) {
	auto old_r_max = m_radius_generator.get_max_radius();
	auto new_r_max = a_radius_generator.get_max_radius();
	if (old_r_max < new_r_max) {
		rsa_mpi::message(("Impossible to impose a smaller maximal radius after initializing the algorithm"
			"\n"
			"This is due to the unrderliying ghost data, which rely on a given r_max"));
		throw runtime_error("Incorrect maximal radius");
	}
	m_r_max = std::max(m_r_max, new_r_max);
	m_radius_generator = std::move(a_radius_generator);
}

template<int DIM>
template<class CenterGenerator, class PriorityGenerator, class RadiusGenerator>
int64_t rsa_algo<DIM>::single_draw(CenterGenerator& center_generator,
	PriorityGenerator& priority_generator,
	RadiusGenerator& a_radius_generator,
	int nb_shots,
	uint64_t nb_spheres_total_max,
	bool may_outreach_nb_spheres) {
	// get data storage
	auto& ghost_areas = m_domain.get_ghost_areas();
	// This function generates a_size spheres and try to store them in the data_storage
	int64_t nb_added_spheres = generate_spheres<DIM>(
		this->get_grid(),
		m_domain.get_recv_buffers(), m_domain.get_send_buffers(), ghost_areas,
		center_generator, a_radius_generator, priority_generator,
		m_ghost_data, nb_shots, nb_spheres_total_max, may_outreach_nb_spheres);
	// some checks used for debugging
	assert(this->get_grid().check_particles() && " collision detected");
	assert(check_no_doublon(this->get_grid()));
	//
	return nb_added_spheres;
}

template<int DIM, int method, bool plog>
rsa_algo<DIM> uniform_generate(rsa_domain<DIM>& a_domain, double a_rad,
	int a_size, int a_n_draw, size_t seed) {
	RadiusGenerator<DIM> radius_generator(
		vector<tuple<double, double, int>>{ {a_rad, 1., 0} }, a_domain.get_total_volume()
	);
	assert(radius_generator.get_current_number() > 0);
	return uniform_generate<DIM, method, plog>(a_domain, radius_generator, a_size, a_n_draw, seed);
}

template<int DIM, int method, bool plog>
rsa_algo<DIM> uniform_generate(rsa_domain<DIM>& a_domain,
	sac_de_billes::RadiusGenerator<DIM>& radius_generator,
	int a_size, int a_n_draw, size_t seed) {
	rsa_algo<DIM> rsaalgo(a_domain, radius_generator, a_size, a_n_draw);
	rsaalgo.template proceed<method, plog>(seed);
	return rsaalgo;
}


template<int DIM>
void check_time_unavoidable_vox(rsa_domain<DIM>& a_domain, double min_radius, double exclusion_distance) {
	sac_de_billes::RadiusGenerator<DIM> radius_generator(vector<tuple<double, uint64_t, int>> {tuple<double, uint64_t, int>{min_radius, 1, 0}}, exclusion_distance);
	rsa_algo<DIM> rsaalgo(a_domain, std::move(radius_generator), 0, 0);
	rsaalgo.check_vox_time();
}

template<int DIM, typename FuncLaw>
arr_vec_double<DIM> auxi::generate_sphere_positions(FuncLaw* a_func_law, int a_size) {
	arr_vec_double<DIM> pos(a_size);
	//! warning, incorrect sense
	for (int d = 0; d < DIM; d++) {
		auto& law = a_func_law[d];
		for (int it = 0; it < a_size; it++) {
			pos[it][d] = law();
		}
	}
	return pos;
}

template<int DIM>
std::vector<law::uniform<double>> auxi::create_random_law(const rsa_domain<DIM>& a_domain,
	std::mt19937& random_generator) {
	std::vector<law::uniform<double>> ulaw{};

	const auto inf = a_domain.get_inf();
	const auto sup = a_domain.get_sup();

	// law are defined for every dimensions
	for (int i = 0; i < DIM; i++) {
		ulaw.emplace_back(law::uniform<double>(inf[i], sup[i], random_generator)); // fill law, using boundaries
	}

	return ulaw;
}

template<int DIM>
vector<uint64_t> auxi::compute_conflict_cells(const rsa_data_storage<DIM>& a_spheres, const rsa_grid<DIM>& a_grid) {
	vector<uint64_t> cells_with_conflicts = {};
	for (size_t i = 0; i < a_spheres.size(); i++) {
		const int cell_id = a_grid.compute_cell_idx(a_spheres.get_center(i));
		cells_with_conflicts.push_back(cell_id);
	}
	std::sort(cells_with_conflicts.begin(), cells_with_conflicts.end());
	auto last = std::unique(cells_with_conflicts.begin(), cells_with_conflicts.end());
	cells_with_conflicts.erase(last, cells_with_conflicts.end());
	return cells_with_conflicts;
}

template<int DIM>
double auxi::magical_default_miss_rate() {
	if constexpr (DIM <= 3) {
		return 0.95; // magical
	} else {
		return 1 - 0.05 * auxi_function::puissance<DIM - 3>(0.1); // magical
	}
}

template<int DIM>
void auxi::recompute_conflict_cells(const rsa_grid<DIM>& a_grid, std::vector<uint64_t>& a_cells_with_conflicts) {
	int current_size = a_cells_with_conflicts.size();
	for (int64_t i = a_cells_with_conflicts.size() - 1; i >= 0; i--) {
		bool is_conflict = false;
		const auto& spheres = a_grid.get_data(a_cells_with_conflicts[i]);
		for (size_t s = 0; s < spheres.size(); s++) {
			if (spheres.template is_tag<TypeTag::Not_Validated>(s)) {
				is_conflict = true;
				break;
			}
		}
		if (not is_conflict) {
			a_cells_with_conflicts[i] = a_cells_with_conflicts[current_size - 1];
			current_size--;
		}
	}
	a_cells_with_conflicts.resize(current_size);
}

template<class INT_TYPE>
std::vector<uint64_t> generate_priority(int a_size, std::mt19937& a_random_generator) {
	auto ulaw = law::uniform<INT_TYPE>(0, std::numeric_limits<INT_TYPE>::max(), a_random_generator);
	std::vector<uint64_t> ret(a_size);
	for (int id = 0; id < a_size; id++) {
		ret[id] = ulaw();
	}
	return ret;
}

template<int DIM>
bool add_to_sample(
	rsa_grid<DIM>& a_grid,
	std::vector<buffer_for_spheres<DIM>>& a_recv, std::vector<buffer_for_spheres<DIM>>& a_send,
	std::vector<rsa_ghost_area<DIM>>& a_ghost_areas,
	std::vector<uint64_t>& a_cells_with_conflicts,
	uint64_t a_maximal_nb_spheres,
	int64_t& local_nb_new_spheres,
	bool may_outreach_nb_spheres) {
	assert((a_grid.template how_many<TypeCell::Real, TypeTag::Skip>() == 0));
	assert((a_grid.template how_many<TypeCell::Real, TypeTag::Conflict>() == 0));
	assert((a_grid.template how_many<TypeCell::Ghost, TypeTag::Not_Validated>() == 0));
	// We test if the particles build_grid with the particles inside the data_storage (we use cells). 
	remove_collisions_with_validated_spheres<DIM>(a_grid, a_cells_with_conflicts);
	auxi::recompute_conflict_cells<DIM>(a_grid, a_cells_with_conflicts);

	// Update data and ghost paritcles with draw edges particles
	fill_ghost<DIM, TypeTag::Undecided>(a_grid, a_ghost_areas);
	put_ghost_into(a_grid, a_ghost_areas, a_recv, a_send, TypeTag::Undecided);

	// detect conflicts (unsolvable conflict are sotred in conflicts) 
	find_solve_conflicts<DIM>(a_grid, a_cells_with_conflicts);

	//! section for ensuring taking only the first cells.
	uint64_t nb_undecided = a_grid.template how_many<TypeCell::Real, TypeTag::Undecided>(a_cells_with_conflicts);
	bool has_outreached_nb_spheres = false;
	uint64_t total_nb_undecided = 0;
	if (may_outreach_nb_spheres) {
		total_nb_undecided = rsa_mpi::compute_mpi_sum(nb_undecided);
		has_outreached_nb_spheres = (total_nb_undecided > a_maximal_nb_spheres);
	}

	//! test if would add less than desired number spheres
	if (has_outreached_nb_spheres) { // if more
		int64_t final_nb_undecided = auxi::retain_only_first_elements<DIM>(a_grid, a_cells_with_conflicts, a_maximal_nb_spheres);
		auxi::update_sphere_tags<DIM>(a_grid, a_recv, a_send, a_ghost_areas, a_cells_with_conflicts, final_nb_undecided);
		local_nb_new_spheres += final_nb_undecided;
		a_maximal_nb_spheres = 0;
		return false;
	} else { // if less
		auxi::update_sphere_tags<DIM>(a_grid, a_recv, a_send, a_ghost_areas, a_cells_with_conflicts, nb_undecided);
		local_nb_new_spheres += nb_undecided;
		a_maximal_nb_spheres -= total_nb_undecided;
		// in case no more spheres should be placed
		if (a_maximal_nb_spheres == 0) {
			a_grid.template remove_tagged_spheres<TypeCell::All, TypeTag::Not_Validated>(a_cells_with_conflicts);
			return false;
		} else {
			// test if there is at least one conflict
			int64_t local_nb_conflict = a_grid.template how_many<TypeCell::Real, TypeTag::Not_Validated>(a_cells_with_conflicts);
			int64_t global_nb_conflicts = rsa_mpi::compute_mpi_sum(local_nb_conflict);
			bool is_conflict = (global_nb_conflicts != 0);
			if (is_conflict) {
				auxi::recompute_conflict_cells<DIM>(a_grid, a_cells_with_conflicts);
				return true;
			} else {
				return false;
			}
		}
	}
}


template<int DIM>
void auxi::update_sphere_tags(rsa_grid<DIM>& a_grid,
	std::vector<buffer_for_spheres<DIM>>& a_recv, std::vector<buffer_for_spheres<DIM>>& a_send,
	std::vector<rsa_ghost_area<DIM>>& a_ghost_areas,
	std::vector<uint64_t>& a_cells_with_conflicts,
	uint64_t nb_undecided) {
	// Update data and ghost paritcles with draw edges particles
	fill_ghost<DIM, TypeTag::Undecided>(a_grid, a_ghost_areas);
	put_ghost_into(a_grid, a_ghost_areas, a_recv, a_send, TypeTag::Validated);
	a_grid.template remove_tagged_spheres<TypeCell::Ghost, TypeTag::Not_Validated>();
	// update data
	a_grid.template remove_tagged_spheres<TypeCell::Real, TypeTag::Skip>(a_cells_with_conflicts);

	uint64_t nb_spheres = a_grid.get_number_of_spheres_fast();
	a_grid.set_number_of_spheres_fast(nb_spheres + nb_undecided);
	auto local_change = [](int cell_id, auto& cell) {
		for (size_t s = 0; s < cell.size(); s++) {
			if (cell.template is_tag<TypeTag::Undecided>(s)) {
				cell.template set_tag<TypeTag::Validated>(s);
			} else if (cell.template is_tag<TypeTag::Conflict>(s)) {
				cell.template set_tag<TypeTag::Undecided>(s);
			}
		}
		};
	a_grid.template apply_on_cells<TypeCell::Real, true>(a_cells_with_conflicts, local_change);
	// Check
	assert(check_no_doublon(a_grid));
}

template<int DIM>
int64_t auxi::retain_only_first_elements(rsa_grid<DIM>& a_grid,
	std::vector<uint64_t>& a_cells_with_conflicts,
	const uint64_t& a_maximal_nb_spheres) {

	//! reconstruct list of litigious spheres
	vec_int priority{};
	std::vector<tuple<int64_t, int>> cell_id_sphere_id{};
	auto add_cell_id_sphere_id = [&priority, &cell_id_sphere_id](int cell_id, auto& cell) {
		for (size_t s = 0; s < cell.size(); s++) {
			if (cell.template is_tag<TypeTag::Undecided>(s)) {
				priority.push_back(cell.get_priority(s));
				tuple<int64_t, int> new_ids(cell_id, s);
				cell_id_sphere_id.push_back(new_ids);
			}
		}
		};
	a_grid.template apply_on_cells<TypeCell::Real, false>(a_cells_with_conflicts, add_cell_id_sphere_id);
	//! get the pivotal priority
	double pivot = median_of_medians::find_pivot_for_first_elements(priority, a_maximal_nb_spheres, std::less<uint64_t>{});
	//!
	//! update a_grid, according to it
	//! 1) skip unchosen spheres
	int64_t nb_undecided = 0;
	for (size_t i = 0; i < cell_id_sphere_id.size(); i++) {
		if (priority[i] > pivot) {
			auto cell_id = get<0>(cell_id_sphere_id[i]);
			int s_id = get<1>(cell_id_sphere_id[i]);
			a_grid.template set_tag<TypeTag::Skip>(cell_id, s_id);
		} else {
			nb_undecided++;
		}
	}
	//! 2) mark all sphere not undecided or validated as Skip
	auto tag_unchosen_spheres = [](int cell_id, auto& cell) {
		for (size_t s = 0; s < cell.size(); s++) {
			if (not cell.template is_tag<TypeTag::Undecided>(s) and not cell.template is_tag<TypeTag::Validated>(s)) {
				cell.template set_tag<TypeTag::Skip>(s);
			}
		}
		};
	a_grid.template apply_on_cells<TypeCell::Real, true>(a_cells_with_conflicts, tag_unchosen_spheres);
	//!
	return nb_undecided;
}

void print_draw_size(int a_draw_size) {
	int global = 0;
	MPI_Allreduce(&a_draw_size, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	rsa_mpi::message("--> Number of spheres in the sample", global);
}

void print_log(int number_of_conflicts, int number_of_particles) {
	int global[2] = { 0,0 };
	int local[2] = { number_of_conflicts,number_of_particles };
	MPI_Allreduce(local, global, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	rsa_mpi::message("---> Number of spheres (including ghost)", global[1], "- Number of conflicts:", global[0]);
}

template<int DIM, class CenterGenerator, class PriorityGenerator, class RadiusGenerator>
int64_t generate_spheres(rsa_grid<DIM>& a_grid,
	Buffers<DIM>& a_recv, Buffers<DIM>& a_send,
	GhostAreas<DIM>& a_ghost_areas,
	CenterGenerator& a_center_generator,
	RadiusGenerator& a_radius_generator,
	PriorityGenerator& a_priority_generator,
	rsa_data_storage<DIM>& a_ghost_data,
	int a_size,
	uint64_t nb_spheres_total_max,
	bool may_outreach_nb_spheres) {
	rsa_data_storage<DIM> spheres = generate_candidates<DIM>(
		a_center_generator, a_radius_generator, a_priority_generator, a_size);

	int64_t local_nb_new_spheres = 0;

	//
	build_grid(a_grid, spheres, TypeTag::Undecided);

	// identify the cells containing spheres with conflicts
	vector<uint64_t> cells_with_conflicts = auxi::compute_conflict_cells<DIM>(spheres, a_grid);

	while (add_to_sample(a_grid,
		a_recv, a_send,
		a_ghost_areas, cells_with_conflicts,
		nb_spheres_total_max, local_nb_new_spheres, may_outreach_nb_spheres)) {
	}

	assert((a_grid.template how_many<TypeCell::All, TypeTag::Not_Validated>() == 0));
	assert(a_grid.check_particles() && " collision detected");

	return local_nb_new_spheres;
}

template<int DIM, class CenterGenerator, class PriorityGenerator, class RadiusGenerator>
rsa_data_storage<DIM> generate_candidates(
	CenterGenerator& a_center_generator,
	RadiusGenerator& a_radius_generator,
	PriorityGenerator& a_priority_generator,
	int a_size) {

	rsa_data_storage<DIM> spheres;
	auto pos = a_center_generator(a_size);
	auto priority = a_priority_generator(a_size);
	auto phases_radii = a_radius_generator(a_size);
	spheres.add_spheres(priority, std::get<0>(phases_radii), std::get<1>(phases_radii), pos);

	return spheres;
}

} // namespace algorithm
