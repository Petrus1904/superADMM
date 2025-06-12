function opts = getDefaultSettings()
%opts = GETDEFAULTSETTINGS() returns the default settings of the superADMM solver
%   see also: superADMM, setup 

% (c) Peter Verheijen, 2025

opts.verbose = 0;
opts.maxIter = 500;
opts.sigma = 1e-6;
opts.alpha = 500;
opts.tau = 0.5;
opts.rho_0 = 1;
opts.RBound = 1e8;
opts.eps_abs = 1e-8;
opts.eps_inf = 1e-8;
opts.repInterval = 10;
opts.timeLimit = 0;

end