//! Copyright : Apache 2.0, see LICENSE 
//! 
#pragma once

template<int DIM, TypeTag typeTag>
void fill_ghost(const rsa_grid<DIM>& a_grid,
    std::vector<rsa_ghost_area<DIM>>& a_ghost) {
    const int ghost_size = a_ghost.size();

    auto fill_ghost_cell = [ghost_size, &a_ghost](int cell_id, const auto& spheres) {
        for (size_t id = 0; id < spheres.size(); id++) {
            if constexpr (typeTag == TypeTag::Any) {
                for (int g_idx = 0; g_idx < ghost_size; g_idx++) {
                    rsa_ghost_area<DIM>& ghost = a_ghost[g_idx];
                    ghost.add_particle(id, spheres);
                    // A particle could be in different ghost areas that's why we have removed the following line
                }
            } else if (spheres.template is_tag<typeTag>(id)) {
                for (int g_idx = 0; g_idx < ghost_size; g_idx++) {
                    rsa_ghost_area<DIM>& ghost = a_ghost[g_idx];
                    ghost.add_particle(id, spheres);
                }
            }
        }
        };
    a_grid.template apply_on_cells<TypeCell::Edge, false>(fill_ghost_cell);
}

template<int DIM>
void fill_buffer_with_ghost(rsa_ghost_area<DIM>& a_ghost, buffer_for_spheres<DIM>& a_buffer) {
    a_ghost.apply_periodicity();
    auto& ghost_data = a_ghost.get_data();
    const int size = ghost_data.size();
    a_buffer.set_size(size);
    a_buffer.pack(ghost_data);
    a_ghost.purge();
    assert(a_ghost.size() == 0);
}

template<int DIM>
void update_ghost_impl(
    std::vector<buffer_for_spheres<DIM>>& a_recv,
    std::vector<buffer_for_spheres<DIM>>& a_send,
    std::vector<rsa_ghost_area<DIM>>& a_ghost,
    rsa_data_storage<DIM>& a_data) {
    // mpi stuff
    int mpi_size = -1;
    int mpi_rank = -1;
    auto comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    int tag_size = 0;
    int tag_data = mpi_size;
    std::vector<MPI_Request> request(2 * (a_send.size() + a_recv.size()));

    // Update ghost are done in 6 steps:
    // (1) send the number of particles per ghost area.
    // (2) receive these numbers
    // (3) pack data in buffers.
    // (4) send particles.
    // (5) receive particles in buffers.
    // (6) unpack particles in buffers in the data storage.

    // 
    int acc = 0;
    // step (1)
    for (size_t idx = 0; idx < a_send.size(); idx++) {
        const int size = a_ghost[idx].size();
        a_send[idx].set_size(size);
        const int tag = tag_size + mpi_rank;
        a_send[idx].send_size(tag, comm, &request[acc++]);
    }

    // step (2)
    for (auto& it : a_recv) {
        int from = it.get_proc_id();
        const int tag = tag_size + from;
        it.recv_size(tag, comm, &request[acc++]);
    }
    assert(acc == a_send.size() + a_recv.size());
    int tmp = acc;
    MPI_Waitall(tmp, request.data(), MPI_STATUSES_IGNORE);

    for (size_t idx = 0; idx < a_send.size(); idx++) {
        // step (3)
        if (a_send[idx].size() == 0) continue;
        // m_ghost data are purged after this step
        fill_buffer_with_ghost<DIM>(a_ghost[idx], a_send[idx]);
        // step (4)
        const int tag = tag_data + mpi_rank;
        a_send[idx].send_data(tag, comm, &request[acc++]);
    }

    // step (5)
    for (size_t idx = 0; idx < a_recv.size(); idx++) {
        if (a_recv[idx].size() == 0) continue;
        const int tag = tag_data + a_recv[idx].get_proc_id();
        a_recv[idx].recv_data(tag, comm, &request[acc++]);
    }

    MPI_Waitall(acc - tmp, request.data() + tmp, MPI_STATUSES_IGNORE);

    int new_particles = 0;
    for (size_t idx = 0; idx < a_recv.size(); idx++) {
        new_particles += a_recv[idx].size();
    }

    // step (6)
    for (size_t idx = 0; idx < a_recv.size(); idx++) {
        a_recv[idx].unpack(a_data);
    }
    assert(a_data.check_radius());
}

template<int DIM>
void put_ghost_into(rsa_grid<DIM>& a_grid, GhostAreas<DIM>& a_ghost_areas,
    Buffers<DIM>& a_recv, Buffers<DIM>& a_send, TypeTag typeTag) {
    rsa_data_storage<DIM> ghost;
    update_ghost_impl(a_recv, a_send, a_ghost_areas, ghost);
    build_grid(a_grid, ghost, typeTag);
}
