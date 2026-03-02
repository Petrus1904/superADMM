from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import platform
import struct
from pathlib import Path
import shutil
import subprocess
import scipy_openblas64 as sob


print("OS:", platform.system())
print("OS version:", platform.release())
print("Platform details:", platform.platform())
print("Machine:", platform.machine())
print("Python bitness:", struct.calcsize("P") * 8, "bit")



#bron_map = Path(scipy_openblas64.get_lib_dir())



if platform.system() == "Darwin":
    sob_lib_dir_path = Path(sob.get_lib_dir()).parent
else:
    sob_lib_dir_path = Path(sob.get_lib_dir())

doel_map = Path("./superADMM")




files_moved = []
for bestand in sob_lib_dir_path.rglob("*.dylib" if platform.system() == "Darwin" else "*"):
    if not bestand.is_file():
        continue
    print("File moved:", bestand)
    files_moved.append(shutil.copy(bestand, doel_map / bestand.name))
    print(sob_lib_dir_path.parent)

if platform.system() == "Darwin":
    result = subprocess.run(
        ["install_name_tool", "-change", "@loader_path/../.dylibs/libgfortran.5.dylib", "@loader_path/libgfortran.5.dylib", str(doel_map)+"/"+sob.get_library(fullname=True)],
        #["echo", str(doel_map)+"/"+sob.get_library(fullname=True)],
        capture_output=True,
        text=True,
        check=True
    )

    print(result)

    result = subprocess.run(
        ["install_name_tool", "-change", "@loader_path/../.dylibs/libgcc_s.1.1.dylib", "@loader_path/libgcc_s.1.1.dylib", str(doel_map)+"/"+sob.get_library(fullname=True)],
        #["echo", str(doel_map)+"/"+sob.get_library(fullname=True)],
        capture_output=True,
        text=True,
        check=True
    )

    print(result)

    result = subprocess.run(
        ["install_name_tool", "-change", "@loader_path/../.dylibs/libquadmath.0.dylib", "@loader_path/libquadmath.0.dylib", str(doel_map)+"/"+sob.get_library(fullname=True)],
        #["echo", str(doel_map)+"/"+sob.get_library(fullname=True)],
        capture_output=True,
        text=True,
        check=True
    )

    print(result)

extra_link_args = []
extra_compile_args = []
if platform.system() == "Windows":
    pass
elif platform.system() == "Linux":
    extra_link_args.append("-Wl,-rpath,$ORIGIN")
    extra_compile_args.append("-g")
    extra_compile_args.append("-O0")
elif platform.system() == "Darwin":
    extra_link_args.append("-Wl,-rpath,@loader_path")

ext = Extension(
    "superADMM.superADMM",
    sources=["superADMM/pysuperADMM.pyx", "src/superADMM.c", "src/csparse.c", "src/ldl.c"],
    include_dirs=["src", np.get_include(), sob.get_include_dir()], # 
    library_dirs=["src", "./superADMM"],
    libraries=[sob.get_library()],
    extra_link_args=extra_link_args,
    define_macros=[("BUILDING_PYTHON", "1")],
    extra_compile_args=extra_compile_args
)

    # libraries=["mijnlib"],
    # library_dirs=["pad/naar/lib"],
    # 

setup(
    name="superADMM",
    version="0.7.0",
    packages=["superADMM"],
    package_data={"superADMM": ["libopenblas.dll",sob.get_library()]}, # bundle DLL 
    include_package_data=True,
    install_requires=[ "numpy>=2.0.2", "scipy>=1.13.1","scipy_openblas64"],
    ext_modules=cythonize(
        ext,
        compiler_directives={"language_level": "3"},
    )
)

for file in files_moved:
    file.unlink()