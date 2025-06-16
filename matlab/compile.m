% compile.m -- Compiles the superADMM solver for first use
%   
%   see also: superADMM, getDefaultSettings

% (c) Peter Verheijen, 2025

mex -R2018a superADMM_mex.c superADMM.c csparse.c ldl.c ccBlas.c -output superADMM -lmwblas -lmwlapack -DMATLAB_COMP
