% compile.m -- Compiles the superADMM solver for first use
%   
%   see also: superADMM, superADMM_preCompute, superADMM_getOpts

% (c) Peter Verheijen, 2026

mex -R2018a ..\src\superADMM_mex.c ..\src\superADMM.c ..\src\csparse.c ..\src\ldl.c ..\src\ccBlas.c -output superADMM -lmwblas -lmwlapack -DMATLAB_COMP

thisFile = mfilename('fullpath');
thisFolder = fileparts(thisFile);
addpath(thisFolder);
status = savepath;

if status == 0
    fprintf('superADMM compile completed.\n');
else
    fprintf('superADMM compile failed, try running again as Administrator,\n or add this folder manually to the default paths (Home->Environments->Set Path). \n');
end
