# SuperADMM
Quadratic Program Solver with dynamic weighting ADMM

This algorithm specifically solves problems of shape
```
minimize_x   0.5 x'Px + x'q
subject to   l <= Ax <= u
```
where `x in R^n` is the optimization variable, `P in R^(n x n)` and `q in R^n` describe the quadratic cost function (P is positive semi-definite), `A in R^(m x n)` is the linear constraint mapping and `l in R^m`, `u in R^m` denote the lower and upper bounds of the constraints. Note that by setting `l_i = u_i` for some index `i`, one can also include equality constraints in the problem. Furthermore, `l_i = -inf` or `u_i = inf` allows to eliminate a lower or upper bound, respectively.

Code will follow shortly

# Installation
Installation is currently only available in Windows, Linux and MACOS will follow later
## Python
- Download or clone this package
- install `gcc` (if you havent already)
- Download and unpack **OpenBLAS** binaries (OpenBLAS-x.x.xx_x86.zip): https://github.com/OpenMathLib/OpenBLAS/releases/tag/v0.3.29
- Open Command Prompt at the superADMM folder
- `gcc -shared -o superADMM.dll superADMM.c csparse.c ldl.c -I"C:\PATHtoOPENBLAS\include" -L"C:\PATHtoOPENBLAS\lib" -lopenblas -DBUILD_DLL`
- Enjoy your fast solver in Python with `import superADMM`

## MATLAB (Express installation)
- Download the latest release
- unzip the contents in some folder
- optionally - add this folder to the default paths in MATLAB
- Enjoy your fast solver in MATLAB with `superADMM.m`

## MATLAB (Manual installation)
- Download or clone this package
- Ensure that MATLAB `mex` add-on (code-generation) is installed, and `gcc` is installed as mex compiler
- Verify with `mex -setup` that `gcc` is the default compiler
- run `setup.m`
- Enjoy your fast solver in MATLAB with `superADMM.m`
# How to Cite
include bibtex
