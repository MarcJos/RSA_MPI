=================
The RSA algorithm
=================

Principle
---------

This software provides an HPC implementation **equivalent** to the
following classical RSA algorithm: 

1. The desired radii of the spheres are sorted in nonincreasing order. 
2. For each radius R, the following happens: 

- if no more sphere can be placed ``break``, 
- a candidate sphere S of radius R is randomly chosen inside the global shape, 
- if S intersects an already placed sphere, then, go to the previous line. Otherwise, it is accepted. 

3. The algorithm terminates if all the desired spheres has been placed, or if no more sphere can be placed (then, the configuration of spheres is guaranteed to be *packed*).

.. note::

  Although it provably provides the same output (in terms of statistics) as the above algorithm, this software uses a more educated strategy, see “Parallel and bias-free RSA algorithm for maximal Poisson-sphere sampling, Josien & Prat (in preparation)”.

Special features
----------------

-  **Inputs**

   -  *Geometry* : The default geometry is a periodic cuboid, of
      prescribed lengths.
   -  *Polydispersity* : There exist various ways to impose
      polydispersity.

      -  Impose an explicit list of radii as a final configuration.
      -  Impose a list of radii :math:`r_1> ... > r_n`, with prescribed
         volume fractions :math:`\phi_1, ..., \phi_n` as a final
         configuration.

-  **Behaviour**:

   -  *Packed configuration* : One can prescribe a list of radii
      :math:`r_1> ... > r_n`, with prescribed volume fractions
      :math:`\phi_1, ..., \phi_n` as a final configuration. Three cases
      are possible:

      -  All the volumes fractions :math:`\phi_j` for :math:`j \leq n`
         can be reached. Then, the algorithm *terminates normally*, and
         returns a message
         ``End of generation, due to lack of new radii``.
      -  All the volumes fractions :math:`\phi_j` for :math:`j` with
         :math:`n > j` can be reached but the last :math:`\phi_n`. Then,
         the algorithm *terminates normally*, and returns a message
         ``End of algorithm : fully packed``. This guarantees that no
         more sphere of radius :math:`r\geq r_n` can be added. (To
         ensure a final packed configuration, impose
         :math:`\phi_n=1.0`.)
      -  At step :math:`j` for :math:`n>j`, it is not possible to place
         more spheres to reach the volume fraction :math:`\phi_j` with
         spheres of radius :math:`r_j`. It leads to an **infinite loop**
         :warning:. (The reason is that the algorithm is not allowed to
         consider radii :math:`r_{j+1}` before reaching :math:`\phi_j`.)
         To avoid this, consider *re-iterating*.

   -  *Re-iterate* : Starting from a previous configuration, it is
      possible to add more spheres with different radius. For example,
      this may be used to get highly heterogenous and packed
      configurations of spheres of radii :math:`r_1 > ... > r_n`, while
      guaranteeing at each step that the configuration is packed.
   -  *Randomness* : The produced configurations are pseudo-random.
      Hence, we use a pseudo-random generator with a random seed. This
      guarantees the following behaviour: for a fixed number of MPI
      processes, and for a fixed seed value, the programm returns the
      same configuration.

 
.. warning::

  This guarantee is for a fixed machine, with a fixed  software environment (compiler, MPI library, …).

.. warning:: 

  Changing the number of MPI processes changes the final configuration.

-  **Outputs**

   -  *Data* : List of spheres can be obtained directly in C++/python,
      on each MPI process.
   -  *.vtk files* : .vtk files (read by
      `Paraview <https://www.paraview.org/>`__).
