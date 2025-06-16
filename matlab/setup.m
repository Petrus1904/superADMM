% setup.m -- Sets up the use of the superADMM solver
%   Effectively, it just adds the path to the default paths :)
%   
%   see also: superADMM, getDefaultSettings

% (c) Peter Verheijen, 2025

thisFile = mfilename('fullpath');
thisFolder = fileparts(thisFile);
addpath(thisFolder);
status = savepath;

if status == 0
    fprintf('superADMM setup completed.\n');
else
    fprintf('superADMM setup failed, try running again as Administrator,\n or add this folder manually to the default paths (Home->Environments->Set Path). \n');
end
