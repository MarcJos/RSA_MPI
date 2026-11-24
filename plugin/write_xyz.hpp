#pragma once

#include <rsa_data_storage.hxx>
#include <rsa_buffer.hxx>
#include <array>
#include <fstream>
#include <sstream>
#include <string>

template<int DIM>
rsa_data_storage<DIM> gather_data(rsa_data_storage<DIM>& a_data) {
	using Buffer = buffer_for_spheres<DIM>;
	// mpi stuff
	int mpi_size;
	int mpi_rank;
	auto comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);

	rsa_data_storage<DIM> data;

	// Particles are gathered on the master mpi process : 0
	// Three steps :
	// 	--> master rank (1) get the number of particles  (2) get particles (3) unpack particles in the data_storage
	// 	--> slave ranks (1) send the number of particles (2) pack particles in a buffer (3) send particles 
	if (mpi_rank != 0) // slaves
	{
		MPI_Request request[2];
		int tag_size = mpi_rank;
		int tag_data = mpi_size + mpi_rank;

		// (1)
		Buffer send(0, int(a_data.size()));
		send.send_size(tag_size, comm, &request[0]);
		// (2)
		send.pack(a_data);
		// (3)
		send.send_data(tag_data, comm, &request[1]);

		// not really useful
		MPI_Waitall(2, request, MPI_STATUSES_IGNORE);
	} else // master
	{
		copy_particles(a_data, data);
		std::vector<Buffer> recvs;
		int current_size = int(a_data.size()); // current size
		std::vector<MPI_Request> request(mpi_size - 1);
		recvs.resize(mpi_size - 1);

		// (1)
		// recv size step
		for (int from = 1; from < mpi_size; from++) {
			Buffer& recv = recvs[from - 1];
			int tag = from;
			recv.set_size(0);
			recv.set_proc_id(from);
			recv.recv_size(tag, comm, &request[from - 1]);
		}
		MPI_Waitall(mpi_size - 1, request.data(), MPI_STATUSES_IGNORE);
		int total_size = current_size;
		for (auto it : recvs) total_size += it.size();
		// (2)
		// recv data step
		for (int from = 1; from < mpi_size; from++) {
			auto& recv = recvs[from - 1];
			int tag = mpi_size + from;
			recv.recv_data(tag, comm, &request[from - 1]);
		}

		MPI_Waitall(mpi_size - 1, request.data(), MPI_STATUSES_IGNORE);

		// (3)
		// The third step unpacks particles in the data storage
		for (int from = 1; from < mpi_size; from++) {
			auto& recv = recvs[from - 1];
			recv.unpack(data);
		}
	}
	return data;
}

template<int DIM>
void write_xyz(std::string a_name, rsa_data_storage<DIM>& spheres, const std::array<double, DIM>& a_sup) {
	int mpi_rank;
	auto comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &mpi_rank);
	const int type = 0;
	auto all_spheres = gather_data(spheres);
	if (mpi_rank == 0) {
		std::ofstream outFile(a_name);
		if (!outFile) {
			std::cerr << "Erreur : impossible de créer le fichier de sortie " << a_name << " !" << std::endl;
			return;
		}
		// XYZ files expect 3 coordinates: pad with 0 when DIM < 3.
		constexpr int n_coords = DIM < 3 ? 3 : DIM;
		auto write_coords = [&outFile](const auto& a_vec) {
			for (int d = 0; d < n_coords; d++) {
				outFile << " " << (d < DIM ? a_vec[d] : 0.0);
			}
			outFile << "\n";
		};
		outFile << all_spheres.size() << "\n";
		std::ostringstream header;
		for (int d = 0; d < n_coords; d++) {
			header << (d == 0 ? "" : " ") << (d < DIM ? a_sup[d] : 0.0);
		}
		outFile << header.str() << "\n";
		for (size_t i = 0; i < all_spheres.size(); i++) {
			const auto center = all_spheres.get_center(i);
			outFile << type;
			write_coords(center);
		}
	}
}
