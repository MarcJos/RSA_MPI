from spack import *

class Rsampi(CMakePackage):
    """C++ library : HPC implementation of RSA algorithm for maximal Poisson-sphere sampling 
		"""

    homepage = "https://github.com/MarcJos/RSA_MPI"
    git = "https://github.com/MarcJos/RSA_MPI.git"

    version("1.0", commit="3caa41dd77c1880802e9f38e6f4399b380ea995d")

    depends_on("cmake@3.26.3")
    depends_on("openmpi")
    build_system("cmake", default="cmake@3.26.3")

    @run_before("cmake")
    def pre_install(self):
         bash = which("bash")
         bash("./Installation/Pre-install.sh")

    def cmake_args(self):
        args = []
        return args
