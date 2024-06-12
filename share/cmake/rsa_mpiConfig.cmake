if(NOT DEFINED rsa_mpi_DIR)
	set(rsa_mpi_DIR $ENV{RSA_MPI_DIR})
	if(NOT DEFINED rsa_mpi_DIR) #Env is empty
	message(FATAL_ERROR "rsa_mpi_DIR environnement variable is undefined, please defined this variable or use -Drsa_mpi_DIR=path_to_rsa_mpi") 
	endif()
endif()


if(DEFINED rsa_mpi_DIR)
	set(rsa_mpi_INCLUDE_DIRS ${rsa_mpi_DIR}/include)
endif()

# Vérifier si les variables requises sont définies
if (NOT EXISTS ${rsa_mpi_INCLUDE_DIRS}/rsa_decoration.hxx)
    message(FATAL_ERROR "Le fichier RSAMPIConfig.cmake does not define correctly the variable rsa_mpi_INCLUDE_DIRS. please define rsa_mpi_DIR")
endif()
