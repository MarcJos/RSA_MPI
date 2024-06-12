//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction>
void rsa_grid<DIM>::apply_on_cells(LambdaFunction& function) const {
	this->template apply_on_cell_id<typeCell, use_omp>(function, this->m_data);
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction>
void rsa_grid<DIM>::apply_on_cells(LambdaFunction& function) {
	this->template apply_on_cell_id<typeCell, use_omp>(function, this->m_data);
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction>
void rsa_grid<DIM>::apply_on_cells(std::vector<uint64_t>& cells_id, LambdaFunction& function) const {
	this->template apply_on_chosen_cells_id<typeCell, use_omp>(cells_id, function, this->m_data);
}

template<int DIM>
template<TypeCell typeCell, bool use_omp, class LambdaFunction>
void rsa_grid<DIM>::apply_on_cells(std::vector<uint64_t>& cells_id, LambdaFunction& function) {
	this->template apply_on_chosen_cells_id<typeCell, use_omp>(cells_id, function, this->m_data);
}

template<int DIM>
template<TypeCell typeCell, TypeTag typeTag>
int64_t rsa_grid<DIM>::how_many(std::vector<uint64_t>& cells_id) const {
	int64_t local = 0;
	auto counter = [&local](int cell_id, const auto& cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			if (cell.template is_tag<typeTag>(i)) {
#pragma omp atomic
				local += 1;
			}
		}
		};
	this->template apply_on_cells<typeCell, true>(cells_id, counter);
	return local;
}

template<int DIM>
template<TypeCell typeCell, TypeTag typeTag>
int64_t rsa_grid<DIM>::how_many() const {
	int64_t local = 0;
	auto counter = [&local](int cell_id, const auto& cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			if (cell.template is_tag<typeTag>(i)) {
#pragma omp atomic
				local += 1;
			}
		}
		};
	this->template apply_on_cells<typeCell, true>(counter);
	return local;
}

template<int DIM>
template<TypeCell typeCell>
uint64_t rsa_grid<DIM>::get_number_of_spheres() const {
	int res = 0;
	auto add_sphere_cell = [&res](int id_cell, auto& cell_data) {
#pragma omp atomic
		res += cell_data.size();
		};
	this->template apply_on_cells<typeCell, true>(add_sphere_cell);
	return res;
}

template<int DIM>
double rsa_grid<DIM>::local_volume_of_spheres() const {
	double res = 0.;
	auto add_sphere_cell = [&res](int id_cell, auto& cell_data) {
		double local_res = 0.;
		for (size_t i = 0; i < cell_data.size(); i++) {
			local_res += sac_de_billes::sphereTools::volumeSphere<DIM>(cell_data.get_rad(i));
		}
#pragma omp atomic
		res += local_res;
		};
	this->template apply_on_cells<TypeCell::Real, true>(add_sphere_cell);
	return res;
}

template<int DIM>
std::vector<rsa_sphere<DIM>> rsa_grid<DIM>::extract_spheres() const {
	std::vector<rsa_sphere<DIM>> res{};
	auto extract_spheres_into_res = [&res](int64_t cell_id, const auto& cell_data) {
		for (size_t i = 0; i < cell_data.size(); i++) {
			res.push_back(cell_data.create_sphere(i));
		}
		};
	this->apply_on_cells<TypeCell::Real, false>(extract_spheres_into_res);
	return res;
}


template<int DIM>
uint64_t rsa_grid<DIM>::get_number_of_spheres_fast() const {
	// just for validation
	[[maybe_unused]] auto for_test = [this]() {
		int res = 0;
		auto add_sphere_cell = [&res](int id_cell, auto& cell_data) {
			for (size_t i = 0; i < cell_data.size(); i++) {
				if (cell_data.template is_tag<TypeTag::Validated>(i)) {
#pragma omp atomic
					res++;
				}
			}
			};
		this->template apply_on_cells<TypeCell::Real, true>(add_sphere_cell);
		return res;
		};
	assert(number_of_spheres_fast == for_test());
	// just for validation

	return number_of_spheres_fast;
}


template<int DIM>
template<TypeCell typeCell>
rsa_data_storage<DIM> rsa_grid<DIM>::get_conflict_spheres() const {
	rsa_data_storage<DIM> conflicts{};
	auto add_sphere_in_conflict_cell = [&conflicts](auto id, const auto& cell) {
		for (size_t i = 0; i < cell.size(); i++) {
			if (cell.is_conflict(i)) {
				conflicts.add_sphere(i, cell);
			}
		}
		};
	this->template apply_on_cells<typeCell>(add_sphere_in_conflict_cell);
	return conflicts;
}

template<int DIM>
template<TypeCell typeCell, TypeTag typeTag>
rsa_data_storage<DIM> rsa_grid<DIM>::extract_data() const {
	rsa_data_storage<DIM> res;
	if constexpr (typeTag == TypeTag::Any) {
		auto extract_data_cell = [&res, this](int id_cell, const auto& spheres) {
			copy_particles(spheres, res);
			};
		this->template apply_on_cells<typeCell, false>(extract_data_cell);
		return res;
	} else {
		auto extract_data_cell = [&res, this](int id_cell, const auto& spheres) {
			for (size_t i = 0; i < spheres.size(); i++) {
				if (spheres.template is_tag<typeTag>(i)) {
					res.add_sphere(i, spheres);
				}
			}
			};
		this->template apply_on_cells<typeCell>(extract_data_cell);
		return res;
	}
}

template<int DIM>
template<TypeCell typeCell>
void rsa_grid<DIM>::purge() {
	auto purge_cell = [](int cell_id, auto& spheres) {
		spheres.purge();
		};
	this->template apply_on_cells<typeCell, true>(purge_cell);
}

template<int DIM>
template<TypeCell typeCell>
void rsa_grid<DIM>::add_spheres(const rsa_grid<DIM>& a_from) {
	assert(a_from.size() == this->size());
	auto add_sphere_cell = [&a_from](int cell_id, auto& to) {
		const auto& from = a_from.get_data(cell_id);
		copy_particles(from, to);
		};
	this->template apply_on_cells<typeCell, false>(add_sphere_cell);
}

template<int DIM>
void rsa_grid<DIM>::add_single_sphere(int i, const rsa_data_storage<DIM>& a_spheres, TypeTag typeTag) {
	size_t cell_id = this->compute_cell_idx(a_spheres.get_center(i));
	auto& cell_spheres = this->m_data[cell_id];
	cell_spheres.add_sphere(i, a_spheres, typeTag);
}

template<int DIM>
template<TypeCell typeCell, TypeTag typeTag>
void rsa_grid<DIM>::remove_tagged_spheres() {
	auto remove_tagged_spheres_locally = [](int id, auto& cell_data) {
		cell_data.template remove_tagged_spheres<typeTag>();
		};
	this->template apply_on_cells<typeCell, true>(remove_tagged_spheres_locally);
}

template<int DIM>
template<TypeCell typeCell, TypeTag typeTag>
void rsa_grid<DIM>::remove_tagged_spheres(std::vector<uint64_t>& restrict_to_cells_id) {
	auto remove_tagged_spheres_locally = [](int id, auto& cell_data) {
		cell_data.template remove_tagged_spheres<typeTag>();
		};
	this->template apply_on_cells<typeCell, true>(restrict_to_cells_id, remove_tagged_spheres_locally);
}

template<int DIM>
bool rsa_grid<DIM>::check_particles() const {
	auto data = this->template extract_data<TypeCell::All>();
	return data.check_particles();
}

