//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once
#include <filesystem>
#include <iostream>
#include <string>

namespace rsa_paraview {

static int step = 0;

void write_pvtp(std::string basedir, std::string basename, size_t number_of_files) {
	std::string name = basedir + "/" + basename + ".pvtp";
	std::ofstream outFile(name);
	if (!outFile) {
		std::cerr << "Erreur : impossible de créer le fichier de sortie suivant: " << name << std::endl;
		return;
	}
	outFile << " <VTKFile type=\"PPolyData\"> " << std::endl;
	outFile << "   <PPolyData GhostLevel=\"0\">" << std::endl;
	outFile << "    <PPointData Scalar=\"id, tag\">" << std::endl;
	outFile << "      <PDataArray type=\"Int64\" Name=\"Id\"/>" << std::endl;
	outFile << "      <PDataArray type=\"Int32\" Name=\"Tag\"/>" << std::endl;
	outFile << "      <PDataArray type=\"Float64\" Name=\"Radius\"/>" << std::endl;
	outFile << "    </PPointData>" << std::endl;
	outFile << "     <PPoints>" << std::endl;
	outFile << "       <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << std::endl;
	outFile << "     </PPoints> " << std::endl;
	for (size_t i = 0; i < number_of_files; i++) {
		std::string subfile = basename + "/" + basename + "_" + std::to_string(i) + ".vtp";
		outFile << "     <Piece Source=\"" << subfile << "\"/>" << std::endl;
	}
	outFile << "   </PPolyData>" << std::endl;
	outFile << " </VTKFile>" << std::endl;
}

void write_vtp(std::string name, size_t size,
	std::stringstream& buff_r, std::stringstream& buff_id, std::stringstream& buff_tag,
	std::stringstream& buff_radius) {
	name = name + ".vtp";
	std::ofstream outFile(name);
	if (!outFile) {
		std::cerr << "Erreur : impossible de créer le fichier de sortie suivant: " << name << std::endl;
		return;
	}

	outFile << "<VTKFile type=\"PolyData\">" << std::endl;
	outFile << " <PolyData>" << std::endl;
	outFile << "   <Piece NumberOfPoints=\"" << size << "\" NumberOfPolys=\"" << 0 << "\">" << std::endl;
	outFile << "   <PointData>" << std::endl;
	outFile << "     <DataArray type=\"Int64\" Name=\"Id\" format=\"ascii\">" << std::endl;
	outFile << buff_id.rdbuf() << std::endl;
	outFile << "     </DataArray>" << std::endl;
	outFile << "     <DataArray type=\"Int32\" Name=\"Tag\" format=\"ascii\">" << std::endl;
	outFile << buff_tag.rdbuf() << std::endl;
	outFile << "     </DataArray>" << std::endl;
	outFile << "     <DataArray type=\"Float64\" Name=\"Radius\" format=\"ascii\">" << std::endl;
	outFile << buff_radius.rdbuf() << std::endl;
	outFile << "     </DataArray>" << std::endl;
	outFile << "   </PointData>" << std::endl;
	outFile << "   <Points>" << std::endl;
	outFile << "     <DataArray type=\"Float64\" Name=\"\"  NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
	outFile << buff_r.rdbuf() << std::endl;
	outFile << "     </DataArray>" << std::endl;
	outFile << "   </Points>" << std::endl;
	outFile << "  </Piece>" << std::endl;
	outFile << " </PolyData>" << std::endl;
	outFile << "</VTKFile>" << std::endl;
}

template<int DIM>
void paraview(std::string a_name, const rsa_grid<DIM>& a_grid) {
	namespace fs = std::filesystem;
	// mpi stuff
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Define file and directory names
	std::string basename = "rsa_mpi";
	std::string directory = a_name + "/" + basename + "_" + std::to_string(step);
	std::string filename = directory + "/" + basename + "_" + std::to_string(step) + "_" + std::to_string(rank);

	// Create directories
	if (rank == 0) {
		fs::create_directory(a_name);
		fs::create_directory(directory);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	size_t count(0);
	std::stringstream buff_r; // store vertices
	std::stringstream buff_id; // store faces
	std::stringstream buff_tag; // store face offsets
	std::stringstream buff_radius; // store face offsets

	// lambda function to fill stringstream
	auto fill_buffers = [&count, &buff_r, &buff_id, &buff_tag, &buff_radius](const size_t cell_id, auto& cell_data) {
		size_t n_spheres = cell_data.size();
		for (size_t s = 0; s < n_spheres; s++) {
			auto r = cell_data.get_center(s);
			buff_r << r[0] << " " << r[1] << " " << r[2] << " ";
			buff_id << cell_data.get_phase(s) << " ";
			//buff_tag << int8_t(cell_data.get_tag(s)) << " ";
			buff_tag << int(cell_data.get_tag(s)) << " ";
			buff_radius << cell_data.get_rad(s) << " ";
		}
		count += n_spheres;
		};

	// false -> no OpenMP
	a_grid.template apply_on_cells<TypeCell::Real, false>(fill_buffers);

	if (rank == 0) {
		std::string dir = a_name;
		std::string name = basename + "_" + std::to_string(step);
		write_pvtp(dir, name, size);
	}
	write_vtp(filename, count, buff_r, buff_id, buff_tag, buff_radius);

	// just increment step
	step++;
}

template<int DIM>
void paraview(const rsa_domain<DIM>& a_domain, std::string a_name) {
	rsa_paraview::paraview(a_name, a_domain.get_grid());
}

template<int DIM>
void debug_paraview(std::string a_name, rsa_data_storage<DIM>& a_data) {
	using Buffer = buffer_for_spheres<DIM>;
	// mpi stuff
	int mpi_size;
	int mpi_rank;
	auto comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);

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
			recv.unpack(a_data);
		}

		// Particles are written in a VTK file
		a_data.write_spheres(a_name);
	}
}

template<int DIM, TypeCell Tc>
void debug_paraview(std::string a_name, rsa_grid<DIM>& a_grid) {
	rsa_data_storage<DIM> data = a_grid.template extract_data<Tc>();
	debug_paraview<DIM>(a_name, data);
}


template<int DIM>
void debug_paraview(std::string a_name, rsa_grid<DIM>& a_grid, rsa_grid<DIM>& a_draw_grid, rsa_data_storage<DIM>& a_conflicts) {
	// mpi stuff
	int mpi_size;
	int mpi_rank;
	auto comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);

	rsa_data_storage<DIM> old_data = a_grid.template extract_data<TypeCell::Real>();
	rsa_data_storage<DIM> new_data = a_draw_grid.template extract_data<TypeCell::Real>();
	rsa_data_storage<DIM> data;

	// add particles
	{
		int size = old_data.size();
		int shift = 0;
		copy_particles(old_data, data);
		for (int i = 0; i < size; i++) data.set_phase(i, 0); // 
		shift += size;
		size = new_data.size();
		copy_particles(new_data, data);
		for (int i = shift; i < shift + size; i++) data.set_phase(i, 1); // 
		shift += size;
		size = a_conflicts.size();
		copy_particles(a_conflicts, data);
		for (int i = shift; i < shift + size; i++) data.set_phase(i, 2); // 
	}

	assert(data.size() == (a_grid.get_number_of_spheres() + a_draw_grid.get_number_of_spheres() + a_conflicts.size()));
	debug_paraview(a_name, data);
}
}
