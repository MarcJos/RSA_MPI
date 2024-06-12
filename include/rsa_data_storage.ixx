//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM>
rsa_sphere<DIM> rsa_data_storage<DIM>::create_sphere(const size_t a_index) const {
	return rsa_sphere<DIM>(this->get_center(a_index), get_rad(a_index), get_phase(a_index));
}

template<int DIM>
template<typename... Args>
void rsa_data_storage<DIM>::add_sphere(const uint64_t a_priority, const uint64_t a_id, const double a_rad, const Args... a_args) {
	// resize
	int old_size = size();
	resize(old_size + 1);
	// old_size is the index of the new element
	rsa_helper::helper_setter<DIM> setter;
	m_priority[old_size] = a_priority;
	m_phase[old_size] = a_id;
	m_rad[old_size] = a_rad;
	setter(old_size, m_pos, a_args...);
}

template<int DIM>
void rsa_data_storage<DIM>::add_sphere(const uint64_t a_priority, const uint64_t a_id, const double a_rad, const vec_d<DIM>& a_center) {
	// resize
	int old_size = size();
	resize(old_size + 1);
	// old_size is the index of the new element
	m_priority[old_size] = a_priority;
	m_phase[old_size] = a_id;
	m_rad[old_size] = a_rad;
	m_pos[old_size] = a_center;
}

template<int DIM>
void rsa_data_storage<DIM>::add_spheres(const rsa_data_storage<DIM>& a_in) {
	auto self = this;
	const int dst_size = self->size();
	const int src_size = a_in.size();
	if (src_size == 0) {
		return;
	} // nothing to do

	self->resize(dst_size + src_size);

	// define dst
	uint64_t* const dst_id = self->get_priority_ptr();
	double* const dst_rad = self->get_rad_ptr();
	auto* const dst_pos = self->get_position_ptr();

	// define src
	const uint64_t* const src_id = a_in.get_priority_ptr();
	const double* const src_rad = a_in.get_rad_ptr();
	const auto* src_pos = a_in.get_position_ptr();

	// do copy
	std::memcpy(dst_id + dst_size, src_id, src_size * sizeof(dst_id[0]));
	std::memcpy(dst_rad + dst_size, src_rad, src_size * sizeof(dst_rad[0]));
	std::memcpy(dst_pos + dst_size, src_pos, src_size * sizeof(dst_pos[0]));
}
template<int DIM>
void rsa_data_storage<DIM>::add_sphere(const uint64_t a_priority, const rsa_sphere<DIM>& a_sphere) {
	add_sphere(
		a_priority,
		a_sphere.phase,
		a_sphere.radius,
		a_sphere.center
	);
}

template<int DIM>
void rsa_data_storage<DIM>::add_sphere(const uint64_t a_prio, const uint64_t a_id, const double a_rad,
	int a_index_center, const arr_vec_double<DIM>& a_centers) {
	this->add_sphere(a_prio, a_id, a_rad, a_centers[a_index_center]);
}

template<int DIM>
void rsa_data_storage<DIM>::add_sphere(const int a_index, const rsa_data_storage<DIM>& a_data) {
	add_sphere(a_data.get_priority(a_index), a_data.get_phase(a_index), a_data.get_rad(a_index), a_data.get_position_ptr()[a_index]);
}

template<int DIM>
template<class PriorityVector>
void rsa_data_storage<DIM>::add_spheres(const PriorityVector& priority,
	const vec_int& a_phase,
	const std::vector<double>& a_rad,
	const arr_vec_double<DIM>& a_centers) {
	int number_of_spheres = size();
	int size = a_centers.size();
	if (size == 0) { return; } // nothing to do
	if (priority.size() != size or a_phase.size() != size or a_rad.size() != size) {
		std::cerr << __PRETTY_FUNCTION__ << std::endl;
		throw std::runtime_error("Incompatible size");
	}

	auto copy_data = [number_of_spheres, size](auto dest_ptr, auto from_ptr) {
		std::memcpy(dest_ptr + number_of_spheres, from_ptr, size * sizeof(dest_ptr[0]));
		};

	resize(number_of_spheres + size);
	copy_data(get_priority_ptr(), priority.data());
	copy_data(get_phase_ptr(), a_phase.data());
	copy_data(get_rad_ptr(), a_rad.data());
	copy_data(get_position_ptr(), a_centers.data());
}

template<int DIM>
void rsa_data_storage<DIM>::remove_sphere(const size_t a_index) {
	const int new_size = size() - 1;
	m_priority[a_index] = m_priority[new_size];
	m_phase[a_index] = m_phase[new_size];
	m_rad[a_index] = m_rad[new_size];
	m_pos[a_index] = m_pos[new_size];
	resize(new_size);
}

template<int DIM>
void rsa_data_storage<DIM>::shift_centers(const double a_dim, const double a_value) {
	assert(a_dim < DIM);
	// no vectorization here
	const int size = this->size();
	auto centers = this->get_position_ptr();
	for (int idx = 0; idx < size; idx++) {
		centers[idx][a_dim] += a_value;
	}
}

template<int DIM>
void rsa_data_storage<DIM>::purge() {
	m_priority.clear();
	m_phase.clear();
	m_rad.clear();
	m_pos.clear();
}

template<int DIM>
void rsa_data_storage<DIM>::resize(const int a_size) {
	m_priority.resize(a_size);
	m_phase.resize(a_size);
	m_rad.resize(a_size);
	m_pos.resize(a_size);
}

template<int DIM>
void rsa_data_storage<DIM>::print_spheres() const {
	auto number_of_spheres = size();
	assert(size() < 1E9);
	for (size_t i = 0; i < number_of_spheres; i++) {
		auto positions = this->get_position();
		if constexpr (DIM == 2)	rsa_mpi::message(i, positions[0][i], positions[1][i]);
		if constexpr (DIM == 3)	rsa_mpi::message(i, positions[0][i], positions[1][i], positions[2][i]);
	}
}

template<int DIM>
void rsa_data_storage<DIM>::print_particle(const int a_idx) const {
	std::string line = " particle " + std::to_string(a_idx) + " (pos:";
	for (int d = 0; d < DIM; d++) {
		line += std::to_string(m_pos[a_idx][d]) + " ";
	}
	line += "-- rad:" + std::to_string(m_rad[a_idx]) + ")";
	std::cerr << line << std::endl;
}

template<int DIM>
void rsa_data_storage<DIM>::write_spheres(std::string a_name) const {
	if constexpr (DIM < 2 or DIM > 3) {
		throw std::runtime_error("Not possible to vizualize in paraview!");
	}

	std::ofstream outFile(a_name);
	if (!outFile) {
		std::cerr << "Erreur : impossible de créer le fichier de sortie !" << std::endl;
		return;
	}
	outFile << "# vtk DataFile Version 3.0" << std::endl;
	outFile << "Spheres" << std::endl;
	outFile << "ASCII" << std::endl;
	outFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	auto numSpheres = size();
	outFile << "POINTS " << numSpheres << " float" << std::endl;
	for (int i = 0; i < numSpheres; i++) {
		auto sphere = create_sphere(i);
		sac_de_billes::auxi_function::writeVectorToString(sphere.center, outFile, " ");
		if constexpr (DIM == 2) {
			outFile << " 0";
		}
		outFile << std::endl;
	}

	outFile << std::endl;
	outFile << "POINT_DATA " << numSpheres << std::endl;
	outFile << "SCALARS Id int 1" << std::endl;
	outFile << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < numSpheres; i++) {
		outFile << get_phase(i) << std::endl;
	}

	outFile << std::endl;
	outFile << "SCALARS radius double 1" << std::endl;
	outFile << "LOOKUP_TABLE default" << std::endl;
	for (int i = 0; i < numSpheres; i++) {
		outFile << get_rad(i) << std::endl;
	}
	outFile.close();
}

template<int DIM>
void rsa_data_storage<DIM>::print_log() const {
	rsa_mpi::message(" -- start grid log -- ");
	rsa_mpi::message("   --> number of spheres:", size());
	rsa_mpi::message(" --   end grid log -- ");
}

template<int DIM>
bool rsa_data_storage<DIM>::check_collision(int i, int j) const {
	double acc = 0;
	for (int d = 0; d < DIM; d++) {
		double tmp = m_pos[i][d] - m_pos[j][d];
		acc += tmp * tmp;
	}
	double rad = m_rad[i] + m_rad[j];
	bool ret = acc < rad * rad;
	return ret;
}

template<int DIM>
bool rsa_data_storage<DIM>::check_particles() const {
	int size_ = size();
	for (int i = 0; i < size_; i++) {
		for (int j = i + 1; j < size_; j++) {
			bool detect = check_collision(i, j);
			if (detect) {
				std::cout << " Error particle " << i << " and particle " << j << std::endl;
				print_particle(i);
				print_particle(j);
				return false;
			}
		}
	}
	return true;
}

template<int DIM>
bool rsa_data_storage<DIM>::check_intersection(const rsa_sphere<DIM>& a_sphere) const {
	for (int i = 0; i < size(); i++) {
		auto sphere = create_sphere(i);
		if (sac_de_billes::areIntersected<DIM>(a_sphere, sphere)) {
			return false;
		}
	}
	return true;
}

template<int DIM>
bool rsa_data_storage<DIM>::check_radius() const {
	for (int i = 0; i < size(); i++) {
		if (get_rad(i) < 1e-16) return false;
	}
	return true;
}

template<int DIM>
template<TypeTag typeTag>
bool rsa_data_storage_with_tags<DIM>::is_tag(int i) const {
	if constexpr (typeTag == TypeTag::Skip or
		typeTag == TypeTag::Validated or
		typeTag == TypeTag::Undecided or
		typeTag == TypeTag::Conflict) {
		return this->m_tags[i] == typeTag;
	}
	if constexpr (typeTag == TypeTag::Any) {
		return true;
	}
	if constexpr (typeTag == TypeTag::Skip_Or_Undecided) {
		return (is_tag<TypeTag::Skip>(i) or is_tag<TypeTag::Undecided>(i));
	}
	if constexpr (typeTag == TypeTag::Not_Validated) {
		return not is_tag<TypeTag::Validated>(i);
	}
}

template<int DIM>
void rsa_data_storage_with_tags<DIM>::remove_sphere(const size_t a_index) {
	rsa_data_storage<DIM>::remove_sphere(a_index);
	this->remove_tag(a_index);
}


template<int DIM>
void rsa_data_storage_with_tags<DIM>::add_sphere(const int a_idx,
	const rsa_data_storage_with_tags<DIM>& a_data) {
	rsa_data_storage<DIM>::add_sphere(a_idx, a_data);
	this->m_tags.push_back(a_data.get_tag(a_idx));
}

template<int DIM>
void rsa_data_storage_with_tags<DIM>::add_sphere(const int a_idx,
	const rsa_data_storage<DIM>& a_data, TypeTag typeTag) {
	rsa_data_storage<DIM>::add_sphere(a_idx, a_data);
	this->m_tags.push_back(typeTag);
}

template<int DIM>
template<TypeTag typeTag>
void rsa_data_storage_with_tags<DIM>::remove_tagged_spheres() {
	const int size = this->size();
	assert(this->m_tags.size() == size);
	for (int i = size - 1; i >= 0; i--) {
		if (is_tag<typeTag>(i)) {
			this->remove_sphere(i);
		}
	}
}

template<int DIM>
void rsa_data_storage_with_tags<DIM>::remove_tag(const size_t a_index) {
	const int new_size = m_tags.size() - 1;
	m_tags[a_index] = m_tags[new_size];
	m_tags.resize(new_size);
}
