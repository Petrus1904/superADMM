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

% use opts = getSuperADMMopts() to obtain a list of the default settings. However
% this is completely optional. superADMM() can also take a struct with only the options
% you want changed. Do note that the names of the options have to match exactly.
%opts = getSuperADMMopts();
opts.verbose = 1; %turn on iterative display

% test dense
[x,y,eflag,info1] = superADMM(P, q, A, l, u, [], [], opts);
if(eflag == 1 && norm(xt-x, inf) < 1e-7 && norm(A'*(yt-y), inf) < 1e-7)
    %OK
    nbytes = fprintf('superADMM dense solver executed successfully!\n');
else
    error("Test failed, dense solver did not compute correctly");
end

% test sparse
[x,y,eflag,info2] = superADMM(sparse(P), q, sparse(A), l, u, [], [], opts);
if(eflag == 1 && norm(xt-x, inf) < 1e-7 && norm(A'*(yt-y), inf) < 1e-7)
    %OK
    nbytes = fprintf('superADMM sparse solver executed successfully!\n');
else
    error("Test failed, sparse solver did not compute correctly");
end
