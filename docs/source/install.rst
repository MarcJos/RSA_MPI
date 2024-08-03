============
Installation
============

.. warning:: 

  These commands should be done in the root directory.

Installation with Cmake
=======================

1. Get pybind11


.. code-block:: bash

   bash Installation/Pre-install.sh

2. Building RSA MPI


.. code-block:: bash

   bash Installation/Install.sh

3. Testing

.. code-block:: bash

   source Env.sh
   cd build
   ctest
   cd ../

4. Generate **doxygen** documentation (html in \`./doc/’)

.. code-block:: bash

   cd include
   doxygen DoxyDoc
   cd ../

Installation with Spack
=======================

.. warning::

  These commands should be done in the root directory.

Add repository that contains the spack package.

.. code-block:: bash

   spack repo add spack_repo

.. code-block:: bash

   spack install rsampi

The executable can be appealed by ``rsa`` after writing:

.. code-block:: bash

   export rsa=$(spack location -i rsampi)/bin/rsa

How to Run RSA MPI
==================

Source the environment
----------------------

.. code-block:: bash

   source Env.sh

Then, three possibilities can be considered, depending on the complexity
of the experiment: 

- run with a command line in the standard case, 
- use a python script for more refined use, 
- build a ``C++`` script for fine-tuned use.

Command line
------------

Command lines are available with the command:

.. code-block:: bash

   rsa --help

Example of usage for 3D RSA with radius 0.05 inside a cube :math:`[0, 1]^3`.

.. code-block:: bash

   mpirun -n 4 rsa --dim 3 --radius 0.05 --seed 0 --size 10000 --inf 0, 0, 0 --sup 1, 1, 1 --paraview 1

C++ Script Example
------------------

Users may only make use of classes and functions of the namespace user_interface. See the **doxygen** documentation for the C++ manual.

**Examples of script:** 

- `main_script <scripts/main_script.cpp>`__ : The  simplest ``C++`` script. 

- `main_max_frac <scripts/main_max_frac.cpp>`__ : Example of script using re-iteration.

Using Python API
----------------

The Python API is actually a wrapper for the ``C++`` API throught the `pybind11 library <https://github.com/pybind/pybind11/>`__. See the `python manual <doc/pythond_manual.md>`__. For examples, see in `scripts_python <scripts_python>`__, to be launched as:

.. code-block:: bash

   mpiexec -np 2 python3 scripts_python/launch_rsa_mpi.py

Miscellaneous
=============

-  For the license, see the file LICENSE.
-  Please kindly report bugs and issues through the gitlab interface.
-  If you use this software, please consider citing “Parallel and
   bias-free RSA algorithm for~maximal Poisson-sphere sampling, Josien &
   Prat (in preparation)”.
