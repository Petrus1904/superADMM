% compile.m -- Compiles the superADMM solver for first use
%   
%   see also: superADMM, getSuperADMMopts

% (c) Peter Verheijen, 2025

mex -R2018a superADMM_mex.c superADMM.c csparse.c ldl.c ccBlas.c -output superADMM -lmwblas -lmwlapack -DMATLAB_COMP

thisFile = mfilename('fullpath');
thisFolder = fileparts(thisFile);
addpath(thisFolder);
status = savepath;

if status == 0
    fprintf('superADMM compile completed.\n');
else
    fprintf('superADMM compile failed, try running again as Administrator,\n or add this folder manually to the default paths (Home->Environments->Set Path). \n');
end
