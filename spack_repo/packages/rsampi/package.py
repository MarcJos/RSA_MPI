from spack import *

class Rsampi(CMakePackage):
    """C++ library : HPC implementation of RSA algorithm for maximal Poisson-sphere sampling 
		"""

    homepage = "https://github.com/MarcJos/RSA_MPI"
    git = "https://github.com/MarcJos/RSA_MPI.git"

    version("1.0", commit="3c21d84bc448f7cebea5ab3a96074784c0dbe51e")

    depends_on("cmake")
    depends_on("openmpi")
    depends_on("py-pybind11")

    build_system("cmake", default="cmake")

    @run_before("cmake")
    def pre_install(self):
#        with working_dir(self.stage.source_path):
         bash = which("bash")
         bash("./Installation/Pre-install.sh")

    def cmake_args(self):
        args = []
        return args
