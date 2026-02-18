# SuperADMM
Fast Quadratic Program Solver with dynamic weighting ADMM

This algorithm specifically solves problems of shape

$$
    \begin{equation}
        \begin{aligned}
            & \underset{x}{\textbf{minimize}}
            & & \tfrac12 x^\top P x + x^\top q\\
            & \textbf{subject to}
            & & l \le Ax \le u \\
        \end{aligned}
    \end{equation}
$$

where $x \in \mathrm{R}^n$ is the optimization variable, $P \in \mathrm{R}^{n \times n}$ and $q \in \mathrm{R}^n$ describe the quadratic cost function ($P$ is positive semi-definite), $A \in \mathrm{R}^{m \times n}$ is the linear constraint mapping and $l \in\mathrm{R}^m$, $u\in\mathrm{R}^m$ denote the lower and upper bounds of the constraints. Note that by setting $l_i = u_i$ for some index $i$, one can also include equality constraints in the problem. Furthermore, $l_i = -\infty$ or $u_i = \infty$ allows users to only consider lower or upper bounds.

# Installation
Installation is currently only available in Windows, Linux and MACOS will follow later
## Python
At this moment, we are working on changing this install routine to `pip install superADMM`, but as of now, that is yet unavailable. Instead, this should be the correct (and very cumbersome) install procedure.
- Download or clone this package
- Download and unpack **OpenBLAS** binaries (OpenBLAS-x.x.xx_x86.zip): https://github.com/OpenMathLib/OpenBLAS/releases
- In `setup.py`, change `openblas_include` and `openblas_lib` to the folder path where you unpacked OpenBLAS.
- Open Command Prompt at the superADMM folder (that is, the main, not superADMM/superADMM)
- `python -m build --wheel` (this will throw a bunch of warnings, which is fine. However, it can give some errors that some packages might need to be installed, install them if so)
- `pip install dist\superadmm-0.7.0-cpVERSION-cpVERSION-win_amd64.whl`, where VERSION is your python version (for me it says `cp313`). You can find the `.whl` file in the dist folder to check the correct name.
- Enjoy your fast solver in Python with `import superADMM`

## MATLAB (Express installation)
- Download the latest release (https://github.com/Petrus1904/superADMM/releases)
- Unzip the contents in some folder
- Run `superADMM_setup.m`
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

# License
SuperADMM is licensed under LGPL 2.1.
