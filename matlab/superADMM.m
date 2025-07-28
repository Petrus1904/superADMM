% [x, y, eflag, info] = superADMM(P, q, A, l, u, x0, y0, options)
%   QP solver for problems of structure
%
%   min_x 0.5x'*P*x + x'*q
%   s.t.  l<=A*x<=u
%
%   Uses an ADMM adaptation that considers super linear convergence.
%   If your problem has constraints in the M*x <= g; E*x = b format, you can
%   rewrite it as: [-inf*ones(length(gam),1); b] <= [M;E]*x <= [gam;b];
%
%   =Arguments===============================================================
%     P          -- (n x n) Positive semi-definite quadratic cost function matrix
%     q          -- (n x 1) linear cost vector
%     A          -- (m x n) matrix mapping the linear constraints
%     l          -- (m x 1) Vector of lower bounds
%     u          -- (m x 1) Vector of upper bounds
%     x0         -- (n x 1) Optional, the initial guess for the solution x, 
%                   can be empty
%     y0         -- (m x 1) Optional, the initial guess for the dual 
%                   variable y, can be empty
%     options    -- Struct Optional, set of options and parameters for the 
%                   solver (see below)
%   
%   =Returns=(All returns are optional)=====================================
%     x          -- the solution vector
%     y          -- the dual variable
%     eflag*     -- 1: Solved Succesfully (OK)
%                   2: Requested tolerance cannot be achieved (maybe OK)
%                   0: Maximum iterations reached (maybe OK)
%                  -1: Cholesky/LDL solver failed (check definiteness of P 
%                      and or increase sigma)
%                  -2: Problem infeasible or unbounded**
%                  -3: Time limit exceeded
%                  -4: Problem non-convex
%     info       -- a struct with properties:
%                   nIter:   number of iterations
%                   rPrim:   primal condition residual
%                   rDual:   dual condition residual
%                   runtime: solver time elapsed
%                   objVal:  the objective value
%                   status:  a string that shows the solvers exit status
%   
%   =options (struct) parameters============================================
%     verbose     -- 0 to silence
%                    1 to print intermediate output (default)
%                    2 to show elaborate timing
%     maxIter     -- (>0) default 500
%     sigma       -- (>=0) relaxation weight parameter on x, default 1e-6
%     rho_0       -- (>0) the initial rho value, default 1
%     tau         -- (0<tau<1) the exponential decrease rate on the bound, default 0.5
%     alpha       -- (>=1) exponential increase rate on R, default 500
%     RBound      -- (>>1) limits on r, default 1e8
%     eps_abs     -- (>=0) absolute convergence bound, default 1e-8
%     eps_inf     -- (>0) infeasibility bounds, default 1e-8
%     repInterval -- (only if verbose = 1), number of iterations between
%                    reporting the status. Default 10.
%     timeLimit   -- solver time limit, set to <= 0 for unlimited, default 0
%     lowRankPer  -- (0<=lowRankPer<=1) threshold for low rank updates, default 0.05
%   
%    one can create a struct manually, or use getSuperADMMopts() to
%    obtain a full struct with the default settings.
%
%    If both P and A are in a sparse format, a sparse solver is used, which 
%    is commonly much faster. Similarly, a dense implementation is used if 
%    both P and A are dense. If either P or A is sparse, an error is thrown. 
%    The header of the verbose text shows which solver is considered.
%
%    * Regardless of the exitflag, x and y will always contain values. Be
%    aware of any further effects of using incorrect solution can have in
%    the code. Always check the exitflag of the solver.
%
%    ** If the problem is infeasible, the outputs x and y will be the
%    certificate of infeasibility.
%    
%    see also: getSuperADMMopts

% (c) Peter Verheijen, 2025

