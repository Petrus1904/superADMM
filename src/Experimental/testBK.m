% rng(1904);

nx = 15;
N = 50;

n = nx*N;
P = sparse(diag(abs((2*randn(n,1)).^2))); %+1e-6*speye(n);
Ax = sprand(nx,nx, 0.3);
Ax = (Ax - 0.5*(Ax~=0)).^5;
A = [speye(n)] - [sparse(nx,n); [kron(speye(N-1), Ax), sparse((N-1)*nx, nx)]];
A = [speye(n); A];
m = 2*n;
R = blkdiag(-1e8*speye(n), -1e-8*speye(n));
PARA = sparse([P, A'; A, R]);

[L1,D1,Pamd] = ldl(PARA, 'vector');
x0 = randn(n+m,1);
tic;
x1 = PARA\x0;
MA57t = toc;
fprintf("MA57 err: %.8e\n", norm(PARA*x1-x0, inf));

tic;
[x2, L2, D2, P2] = BK_LDL(PARA, x0);
BKt = toc;
fprintf("BK LDL err: %.8e\n", norm(PARA*x2-x0, inf));

tic;
[L,D,Pa] = ldl(PARA, 1e-12, 'vector');
LDLt = toc;
x4 = (L'\(D\(L\x0(Pa))));
Pi = zeros(size(Pa));
Pi(Pa) = 1:length(Pa);
x3 = x4(Pi);

fprintf("standard LDL err: %.8e\n", norm(PARA*x3-x0, inf));

fprintf("| %d | %d | %4.3fms | %4.3fms | %5.3fms | %.1e | %.1e | %.1e  |\n", size(PARA,1), nnz(PARA), MA57t*1000, LDLt*1000, BKt*1000, norm(PARA*x1-x0, inf), norm(PARA*x3-x0, inf), norm(PARA*x2-x0, inf));