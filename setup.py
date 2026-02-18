from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

openblas_include = "C:/OpenBLAS/include"  # Windows
openblas_lib = "C:/OpenBLAS/lib"          # Windows

ext = Extension(
    "superADMM.superADMM",
    sources=["superADMM/pysuperADMM.pyx", "src/superADMM.c", "src/csparse.c", "src/ldl.c"],
    include_dirs=["src", np.get_include(),openblas_include], # 
    library_dirs=["src", openblas_lib],
    libraries=["libopenblas"],
    extra_compile_args=[
        "/O2", "/fp:fast", "/arch:AVX2", "/DNDEBUG"
    ],
)

setup(
    name="superADMM",
    version="0.7.0",
    packages=["superADMM"],
    package_data={"superADMM": ["libopenblas.dll"]}, # bundle DLL 
    include_package_data=True,
    install_requires=[ "numpy>=1.26", "scipy>=1.10"],
    ext_modules=cythonize(
        ext,
        compiler_directives={"language_level": "3"},
    ),
)
