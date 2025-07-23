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

Add repository that contains the spack package.

.. code-block:: bash

   git clone https://github.com/Collab4exaNBody/spack-repos.git
   spack repo add spack-repos

.. code-block:: bash

   spack external find openmpi
   spack install rsampi

The executable can be appealed by ``rsa`` after writing:

.. code-block:: bash

   export rsa=`spack location -i rsampi`/bin/rsa

