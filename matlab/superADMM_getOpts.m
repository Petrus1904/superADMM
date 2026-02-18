function opts = superADMM_getOpts()
%opts = superADMM_getOpts() returns the default settings of the superADMM solver
% for help on the opts, see: superADMM

% (c) Peter Verheijen, 2026

opts.verbose = 0;
opts.maxIter = 500;
opts.sigma = 1e-6;
opts.alpha = 500;
opts.tau = 0.5;
opts.rho_0 = 1;
opts.RBound = 1e8;
opts.eps_abs = 1e-8;
opts.eps_rel = 1e-8;
opts.eps_inf = 1e-8;
opts.repInterval = 10;
opts.timeLimit = 0;
opts.lowRankPer = 0.05;
opts.Pamd = [];

end