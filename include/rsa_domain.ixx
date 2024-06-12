//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
#include "link_to_auxi.hxx"

template<int DIM>
rsa_domain<DIM>::rsa_domain(const vec_d<DIM>& global_inf,
	const vec_d<DIM>& global_sup, const int a_ghost_layer, const double a_rad)
	:m_global_inf(global_inf), m_global_sup(global_sup), m_rad{ a_rad } {
	int mpi_rank = -1;
	int mpi_size = -1;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	assert(mpi_size != -1);
	assert(mpi_rank != -1);

	// define id
	this->m_id = mpi_rank;
	// define ghost layer
	this->m_ghost_layer = a_ghost_layer;

	// load balancing 
	lb(mpi_rank, mpi_size, global_inf, global_sup, a_rad);

	// define data storage
	this->m_grid = rsa_grid<DIM>(a_rad, a_ghost_layer, this->get_inf(), this->get_sup());
	this->check_sufficiently_large(a_rad);
}


template<int DIM>
void rsa_domain<DIM>::lb(
	const int a_id, const int a_n_mpi, const vec_d<DIM>& global_inf, const vec_d<DIM>& global_sup, const double a_rad) {
	int ndims[DIM];
	int periods[DIM];
	MPI_Comm MPI_COMM_CART;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for (int dim = 0; dim < DIM; dim++) {
		periods[dim] = 1;
		ndims[dim] = 0; // do not remove it
	}
	MPI_Dims_create(a_n_mpi, DIM, ndims);
	MPI_Cart_create(MPI_COMM_WORLD, DIM, ndims, periods, true, &MPI_COMM_CART);

	int coord[DIM];
	MPI_Cart_coords(MPI_COMM_CART, rank, DIM, coord);


	vec_d<DIM> domain_size;
	for (int dim = 0; dim < DIM; dim++) {
		domain_size[dim] = global_sup[dim] - global_inf[dim];
		double size = domain_size[dim] / ndims[dim];
		m_inf[dim] = global_inf[dim] + coord[dim] * size;
		m_sup[dim] = global_inf[dim] + (coord[dim] + 1) * size;
	}

	int coord_neig[DIM];
	vec_d<DIM> shift;
	// j'aurais du tout capturer ...
	auto fill = [&shift, &coord_neig, &domain_size, &coord, &periods, &ndims](const auto& relative_pos, const vec_d<DIM>& inf, const vec_d<DIM>& sup, bool& periodic)->bool {
		for (int dim = 0; dim < DIM; dim++) {
			coord_neig[dim] = coord[dim] + relative_pos[dim];
			shift[dim] = 0;
			if (coord_neig[dim] == -1) {
				if (periods[dim] == 1) {
					coord_neig[dim] += ndims[dim];
					shift[dim] = domain_size[dim];
					periodic = true;
				} else {
					return false;
				}
			} else if (coord_neig[dim] == ndims[dim]) {
				if (periods[dim] == 1) {
					coord_neig[dim] -= ndims[dim];
					shift[dim] = -domain_size[dim];
					periodic = true;
				} else {
					return false;
				}
			}
		}
		return true;
		};

	array<array<int, 2>, DIM> limits;
	for (int d = 0; d < DIM; d++) {
		limits[d][0] = -1;
		limits[d][1] = 2;
	}
	auto list_relative_positions = sac_de_billes::getAllIndices<DIM, long>(limits);
	assert(list_relative_positions.size() == sac_de_billes::auxi_function::puissance<DIM>(3));

	for (const auto& relativePosition : list_relative_positions) {
		// check not zero
		bool should_continue = true;
		for (int dim = 0; dim < DIM; dim++) {
			if (relativePosition[dim] != 0) should_continue = false;
		}
		if (should_continue) continue;
		//
		bool not_exist = false;
		bool periodic = false;
		vec_d<DIM> inf, sup;
		auxi_rsa_domain::create_layer_bounds_inplace<DIM>(get_inf(), get_sup(), relativePosition, inf, sup, a_rad);
		not_exist = fill(relativePosition, inf, sup, periodic); // -1 <= r-1 <= 1 // fill coord_neig,inf,sup
		if (not_exist == false) continue;
		int neig;
		MPI_Cart_rank(MPI_COMM_CART, coord_neig, &neig);
		m_send.push_back(Buffer(neig, 0));
		m_recv.push_back(Buffer(neig, 0));
		rsa_ghost_area<DIM> ghost(inf, sup, periodic, shift);
		m_ghost.push_back(ghost);
	}
}

template<int DIM>
void auxi_rsa_domain::create_layer_bounds_inplace(const vec_d<DIM>& domain_inf, const vec_d<DIM>& domain_sup,
	sac_de_billes::DiscPoint<DIM> relative_position, vec_d<DIM>& layer_inf, vec_d<DIM>& layer_sup, double a_rad) {
	for (int dim = 0; dim < DIM; dim++) {
		int pos = relative_position[dim];
		if (pos == -1) {
			layer_inf[dim] = domain_inf[dim];
			layer_sup[dim] = domain_inf[dim] + 2 * a_rad;
		} else if (pos == 0) {
			layer_inf[dim] = domain_inf[dim];
			layer_sup[dim] = domain_sup[dim];
		} else {
			layer_inf[dim] = domain_sup[dim] - 2 * a_rad;
			layer_sup[dim] = domain_sup[dim];
		}
	}
}

template<int DIM>
void rsa_domain<DIM>::check_sufficiently_large(double a_rad) {
	for (int dim = 1; dim < DIM; dim++) {
		if (get_sup()[dim] - get_inf()[dim] < 2 * a_rad) {
			std::cerr << "proc" << rsa_mpi::get_my_rank() << "---->" << __PRETTY_FUNCTION__ << std::endl;
			throw runtime_error("The domain is not sufficiently large!");
		}
	}
}

template<int DIM>
void rsa_domain<DIM>::update_ghost(rsa_data_storage<DIM>& a_ghost_data) {
	update_ghost_impl(this->m_recv, this->m_send, this->m_ghost, a_ghost_data);
}

template<int DIM>
void rsa_domain<DIM>::domain_log() const {
	rsa_mpi::message(" ==== log: domain [", m_id, "] ==== ");
	if (rsa_mpi::is_master_rank()) {
		std::cout << " == m_inf : ";
		sac_de_billes::auxi_function::writeVectorToString(get_inf(), std::cout); std::cout << std::endl;
		std::cout << " == m_sup : ";
		sac_de_billes::auxi_function::writeVectorToString(get_sup(), std::cout); std::cout << std::endl;
	}
	rsa_mpi::message(" == sbuffer size : ", m_send.size());
	rsa_mpi::message(" == rbuffer size : ", m_recv.size());
	rsa_mpi::message(" ================================== ");
}

template<int DIM>
double rsa_domain<DIM>::get_total_volume() const {
	auto lengths = m_global_sup - m_global_inf;
	return sac_de_billes::auxi_function::productOf<double>(lengths);
}

template<int DIM>
double rsa_domain<DIM>::compute_total_volume_of_spheres() const {
	double local_volume = this->get_grid().local_volume_of_spheres();
	double global_volume = rsa_mpi::compute_mpi_sum(local_volume);
	return  global_volume;
}
