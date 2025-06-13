% setup.m -- Compiles the superADMM solver from the source code
%   it will also test that the code compiled properly.
%   
%   see also: superADMM, getDefaultSettings

% (c) Peter Verheijen, 2025

mex -R2018a ../src/superADMM_mex.c ../src/superADMM.c ../src/csparse.c ../src/ldl.c ../src/ccBlas.c -output superADMM -I../src -lmwblas -lmwlapack -DMATLAB_COMP
