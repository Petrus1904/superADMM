# SuperADMM
Quadratic Program Solver with dynamic weighting ADMM

This algorithm specifically solves problems of shape
```
minimize_x   0.5 x'Px + x'q
subject to   l <= Ax <= u
```
where `x in R^n` is the optimization variable, `P in R^(n x n)` and `q in R^n` describe the quadratic cost function (P is positive semi-definite), `A in R^(m x n)` is the linear constraint mapping and `l in R^m`, `u in R^m` denote the lower and upper bounds of the constraints. Note that by setting `l_i = u_i` for some index `i`, one can also include equality constraints in the problem. Furthermore, `l_i = -inf` or `u_i = inf` allows to eliminate a lower or upper bound, respectively.

# Installation
Installation is currently only available in Windows, Linux and MACOS will follow later
## Python
- Download or clone this package
- install `gcc` (if you havent already)
- Download and unpack **OpenBLAS** binaries (OpenBLAS-x.x.xx_x86.zip): https://github.com/OpenMathLib/OpenBLAS/releases
- Open Command Prompt at the superADMM python folder
- `gcc -shared -o superADMM.dll ../src/superADMM.c ../src/csparse.c ../src/ldl.c -I"C:\PATHtoOPENBLAS\include" -L"C:\PATHtoOPENBLAS\lib" -lopenblas -DBUILD_DLL`
- Enjoy your fast solver in Python with `import superADMM`

## MATLAB (Express installation)
- Download the latest release (https://github.com/Petrus1904/superADMM/releases)
- Unzip the contents in some folder
- Run `setup.m`
- Enjoy your fast solver in MATLAB with `superADMM.m`

## MATLAB (Manual installation)
- Download or clone this package
- Ensure that MATLAB `mex` add-on (code-generation) is installed, and `gcc` (MinGW64) is installed as mex compiler
- Verify with `mex -setup` that `MinGW64` is the default compiler
- Run `compile.m`
- Enjoy your fast solver in MATLAB with `superADMM.m`

# How to Cite

**APA:**

Verheijen, P.C.N., Goswami, D., and Lazar, M. (2025). *SuperADMM: Solving Quadratic Programs Faster with Dynamic Weighting ADMM.* arXiv: 2506.11608 [math.OC]. URL: https://arxiv.org/abs/2506.11608

**Bibtex:**
```
@misc{superADMM:Verheijen2025,
      title={{SuperADMM: Solving Quadratic Programs Faster with Dynamic Weighting ADMM}}, 
      author={P. C. N. Verheijen and D. Goswami and M. Lazar},
      year={2025},
      eprint={2506.11608},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2506.11608}, 
}
```

# References
The fast execution of various linear algebraic operations is provided by a set of third party libraries, which we list below:
- [OpenBLAS] Wang Qian, Zhang Xianyi, Zhang Yunquan, Qing Yi, *AUGEM: Automatically Generate High Performance Dense Linear Algebra Kernels on x86 CPUs*, In the International Conference for High Performance Computing, Networking, Storage and Analysis (SC'13), Denver CO, November 2013
- [LAPACK] E. Anderson, Z. Bai, C. Bischof, S. Blackford, J. Demmel, J. Dongarra, J. Du Croz, A. Greenbaum, S. Hammarling, A. McKenney, and D. Sorensen, *LAPACK Users' Guide,* 3rd e. Philadelphia, PA: Society for Industrial and Applied Mathematics, 1999.
- [CSPARSE] T.A. Davis, *Direct Methods for Sparse Linear Systems.* SIAM, 2006
- [LDL] T.A. Davis, *Algorithm 849: A concise sparse cholesky factorization package,* ACM Trans. Math. Softw. vol. 31, no. 4, p. 587-591, Dec. 2005
