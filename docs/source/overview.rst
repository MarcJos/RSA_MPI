================
RSA MPI Overview
================

RSA MPI in a nutshell
=====================

``RSA MPI`` is a HPC library dedicated to generating random configurations
non-intersecting balls, with an unbiased the RSA algorithm. The
sequential strategy based on `Ebeida et al,
2012 <https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-8659.2012.03059.x>`__.
This strategy has been extended to the MPI framework, as decribed in
this paper “Parallel and bias-free RSA algorithm for~maximal
Poisson-sphere sampling, Josien & Prat (in preparation)”.

Our implementation successfully generated more than 12 billions of
spheres over 131,072 ``MPI`` processes in 16 seconds in dimension d=3.

This project is written in ``c++17``.

Contributors
============

- Raphaël Prat (CEA/DES/IRESNE/DEC/SESC/LDOP, HPC)
- Marc Josien (CEA/DES/IRESNE/DEC/SESC/LMCP, Mathematics)
