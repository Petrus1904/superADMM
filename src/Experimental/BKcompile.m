% compile.m -- Compiles the BK LDL solver

% (c) Peter Verheijen, 2025

mex -R2018a ..\src\BKLDL_mex.c ..\src\BK_ldl.c ..\src\csparse.c ..\src\ldl.c -output BK_LDL -DMATLAB_COMP
