# Principles

- The python script can be run in mpi mode.
- Each method and class is ended by subscript "_(N)D" for indicating the dimension N. For example `throw_spheres_3D`. By convention, everything is here indicated with a final subscript `_3D`. Nevertheless, these can be replaced by `_2D`, etc.

# Header

For enabling mpi execution, the .py file should have as header:  
`import mpi4py`  
`mpi4py.rc.threaded = True`  
`mpi4py.rc.thread_level = "funneled" `  
`from mpi4py  import  MPI`  
Then, import the python interface:  
`import rsa_mpi_py as rsa_mpi`


# Fast methods
    
For throwing spheres inside a cube of length `L` of radii $`r_1, ..., r_n`$, with prescribed volume fractions $`\phi_1, ..., \phi_n`$ and phase (=identifier) $`p_1, ..., p_n`$, use:  
- `throw_spheres_3D(L, desired_radius_volumeFraction_phase, seed, write_paraview)`
    - `L` : size of the box $`[0, L]^3`$
    - `desired_radius_volumeFraction_phase` : vector of triples [(double) radius , (double) volumeFraction, (int) phase]
    - `seed` indicates the random seed for the pseudo-random generator
    - `write_paraview` : bool write a paraview output"
    - return The local list of spheres inside the domain :warning: only on a single the MPI processes



# Classes

The principle is the following: a domain is defined, then, the polydispersity of the spheres is chosen, then, the algorithm is launched. Each task is abstracted by a single class.

- `rsa_domain_3D`: Data storage for the spheres.
Pave a cuboid $`[l_{\min}[0], l_{\max}[0]] x ... [l_{\min}[D], l_{\max}[D]]`$.
Inside a domain, spheres are partitionned between different cells (=small cuboids).  
Methods:
    - `rsa_domain_3D([l_min[0], ... l_min[D]], [l_max[0], ... l_max[D]], a_ghost_layer, a_rad)`. Constructor.
        - `a_rad` : the size of each cell in any direction should be larger than `a_rad`
        - `a_ghost_layer`: Number of ghost layers around the cell. Should be taken as 1.
    - `domain_log`: Display.
    - `get_total_volume` : Return the total volume of the (global) cuboid.
    - `paraview` : Print paraview output.
    - `compute_total_volume_of_spheres` : Compute the cumulated volume of all spheres inside the whole domains
    - `extract_spheres` : return list of spheres inside the domain. :warning: only on a single MPI processe
- `radius_generator_3D`: Store the prescribed radius distribution fixed by the user.  
Methods:
    - `radius_generator_3D([[rad[0], N[0], p[0]], ...], exclusion_distance)`: First place `N[0]` spheres of radius `rad[0]`, denoted by phase `p[0]`, then `N[1]` spheres... All spheres are separated by a distance `exclusion_distance`.
    - `radius_generator_3D([[rad[0], phi[0], p[0]], ...], volume_total, exclusion_distance)`: First place `N` spheres of radius `rad[0]` denoted by phase `p[0]`  such that the resulting volume fraction is `phi[0]`, then `N` spheres... All spheres are separated by a distance `exclusion_distance`. The total volume of the domain is `volume_total`.
- `rsa_algo_3D`: Class for execution of the algorithm.  
Methods:
    - `rsa_algo_3D(rsa_domain, radius_generator)`: Constructor
    - `proceed(seed)`: Execute the RSA algorithm, with the seed for the pseudo-random generator.
    - `reset_radius_generator(a_radius_generator)`: Reset the radius generator, but retains the already placed spheres. Proceed can then be appealed once more.

