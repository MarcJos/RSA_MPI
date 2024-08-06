import mpi4py
mpi4py.rc.thread_level = "funneled"
from mpi4py  import  MPI
import time

import rsa_mpi_py as rsa_mpi

rsa_mpi.set_nb_threads(1)

l_min = [0., 0., 0.]
l_max = [1., 1., 1.]
ghost_layer = 1
rad = 0.49/2
nb_shots_single_draw = 500
seed = 10
nb_subdivision_rad = 10
multiplier_rad = 0.75
see_paraview = 1

domain =  rsa_mpi.rsa_domain_3D(l_min, l_max, rad)
domain.domain_log()

def radgen(rad_):
    return rsa_mpi.radius_generator_3D([[rad_, 1, 0]], domain.get_total_volume(), 0.)

def print_here():
    local_volume = 0
    local_volume = domain.compute_total_volume_of_spheres()
    print(local_volume)

radius_generator = radgen(rad)
rsaalgo = rsa_mpi.rsa_algo_3D(domain, radius_generator, nb_shots_single_draw)
rsaalgo.proceed(seed)
print_here()

current_rad = rad
for j in range(nb_subdivision_rad):
    current_rad *= multiplier_rad
    radius_generator = radgen(current_rad)
    rsaalgo.reset_radius_generator(radius_generator)
    rsaalgo.proceed(seed)
    print_here()
    if (see_paraview == 1):
        domain.paraview()


