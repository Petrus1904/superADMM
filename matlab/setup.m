
% setup.m -- Compiles the superADMM solver for first use
%   run this file before using superADMM for the first time.
%   it will also test that the code compiled properly.
%   
%   see also: superADMM, getDefaultSettings

% (c) Peter Verheijen, 2025

% mex -R2018a superADMM_mex.c superADMM.c csparse.c ldl.c BK_ldl.c -output superADMM -IC:\OpenBLAS\include -LC:\OpenBLAS\lib -lopenblas -DMATLAB_COMP
mex -R2018a superADMM_mex.c superADMM.c csparse.c ldl.c ccBlas.c -output superADMM -lmwblas -lmwlapack -DMATLAB_COMP

keyboard;
%% test function with test QP problem
P = 2*eye(5);
q = [-2; -6; -8; -4; -10];
Aeq = [1 1 1 1 1; ...
       1 -1 1 -1 1; ...
      -1 -1 -1 -1 -1; ...
       1 2 3 4 5];

beq = [10;3;-10;20];
A = [Aeq; eye(5)];
l = [-inf*ones(4,1); zeros(5,1)];
u = [beq; 10*ones(5,1)];

xt = [3.4; 3.6; 2.8; 0; 0.2];
yt = [0; 0; 8.4; 3.6; 0; 0; 0; -2; 0];


opts.verbose = 2;
opts.maxIter = 500;
opts.sigma = 1e-6;
opts.rho_0 = 1;
opts.tau = 0.5;
opts.alpha = 1;
opts.RBound = 1e8;
opts.eps_abs = 1e-8;
opts.repInterval = 25;
opts.lowRankPer = 0.05;
opts.timeLimit = -1;

% test dense
[x,y,eflag,info1] = superADMM(P, q, A, l, u, [], [], opts);
if(eflag == 1 && norm(xt-x, inf) < 1e-7 && norm(A'*(yt-y), inf) < 1e-7)
    %OK
    nbytes = fprintf('superADMM dense solver executed successfully!\n');
else
    %error("Setup failed, dense solver did not compute correctly");
end
% test sparse
% x = x + 1e-1*randn(5,1);
% % y = y + 1e-1*randn(9,1);
[x,y,eflag,info2] = superADMM(sparse(P), q, sparse(A), l, u, [], [], opts);
if(eflag == 1 && norm(xt-x, inf) < 1e-7 && norm(A'*(yt-y), inf) < 1e-7)
    %OK
    nbytes = fprintf('superADMM sparse solver executed successfully!\n');
else
    %error("Setup failed, sparse solver did not compute correctly");
end
