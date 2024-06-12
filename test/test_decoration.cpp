//! Copyright : Apache 2.0, see LICENSE 
//! 
#include <rsa_decoration.hxx>

int main(int argc, char** argv) {
	rsa_mpi::init(&argc, &argv);
	rsa_mpi::finalize();
	return MPI_SUCCESS;
}
