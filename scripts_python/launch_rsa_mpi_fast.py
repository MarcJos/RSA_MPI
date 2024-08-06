# mpiexec -n 4 python3 launch_rsa_mpi.py 

import mpi4py
mpi4py.rc.thread_level = "funneled"
from mpi4py  import  MPI
import time

import rsa_mpi_py as rsa_mpi


rsa_mpi.set_nb_threads(1)

L = [1., 1., 1.]
desired_radius_volumeFraction_phase=[[0.05, 0.25, 0], [0.025, 1, 0]]
seed = 0
see_mpi_paraview = False

the_spheres = rsa_mpi.throw_spheres_3D(L, desired_radius_volumeFraction_phase, seed, see_mpi_paraview)



# Example:  see one sphere
comm = MPI.COMM_WORLD
if (comm.Get_rank()==0):
    my_sphere = the_spheres[0]
    print(my_sphere)
    print(my_sphere.center)
    print(my_sphere.radius)
    print(my_sphere.phase)


