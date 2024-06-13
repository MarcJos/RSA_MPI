from spack import *

class Rsampi(CMakePackage):
    """C++ library : HPC implementation of RSA algorithm for maximal Poisson-sphere sampling 
		"""

    homepage = "https://github.com/MarcJos/RSA_MPI"
    git = "https://github.com/MarcJos/RSA_MPI.git"

    version("1.0", commit="0c496ad70586b4fbf81d5e0af94ec4898af17ac9")

    depends_on("cmake")
    depends_on("openmpi")
    depends_on("py-pybind11")

    build_system("cmake", default="cmake")


    def cmake_args(self):
        args = []
        return args
