# mpiexec -n 4 python3 launch_rsa_mpi.py 

import mpi4py
mpi4py.rc.thread_level = "funneled"
from mpi4py  import  MPI
import time

import rsa_mpi_py as rsa_mpi


rsa_mpi.set_nb_threads(1)

l = [0., 0., 0.]
L = [1., 1., 1.]
desired_radius_volumeFraction_phase=[[0.05, 0.25, 0], [0.025, 1, 0]]
rad_max = max(rvp[0] for rvp in desired_radius_volumeFraction_phase)
seed = 0
distance_min = 0.
nb_shots_single_draw = 6000 # number of shots before MPI communications

domain = rsa_mpi.rsa_domain_3D(l, L, rad_max)
radius_generator = rsa_mpi.radius_generator_3D(desired_radius_volumeFraction_phase, domain.get_total_volume(), distance_min)
rsaalgo = rsa_mpi.rsa_algo_3D(domain, radius_generator, nb_shots_single_draw)
rsaalgo.proceed(seed)

the_spheres = domain.extract_spheres()
domain.paraview()





