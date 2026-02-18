function Pamd = superADMM_preCompute(P,A)
    % Pamd = superADMM_preCompute(P,A)
    %
    % Advanced usage of the superADMM solver, allows for more efficient repetive 
    % calls of the superADMM solver if the same problem sparsity structure in 
    % P and A is guaranteed to remain the same. Otherwise, this function is not necessary.
    % This function only works for sparse matrices, not for any dense representation
    %
    % =Arguments===============================================================
    %  P    -- (n x n) sparse positive semi-definite quadratic cost function
    %  A    -- (m x n) sparse matrix mapping the linear constraints
    %
    % =Returns=================================================================
    %  Pamd -- Permutation vector corresponding to the sparsity in P and A.
    %
    % =Usage===================================================================
    %     %Assume P,q,A,l,u are already defined
    %     options.verbose = 0;
    %     options.Pamd = superADMM_preCompute(P,A);
    %     for i = 0:maxIter
    %         %update q,l,u, even the values in P and A. only the locations of the
    %         %non-zeros in P and A must remain static.
    %         [x, y, eflag, info] = superADMM.superADMM(P,q,A,l,u,[],[], options);
    %     end
    %
    % Please note that the same can be achieved with the following example, but 
    % can be considered "less neat":
    %     %Assume P,q,A,l,u are already defined
    %     options.verbose = 0;
    %     for i = 0:maxIter
    %         %update q,l,u, even the values in P and A. only the locations of the
    %         %non-zeros in P and A must remain static.
    %         [x, y, eflag, info] = superADMM.superADMM(P,q,A,l,u,[],[], options);
    %         options.Pamd = info.Pamd;
    %      end
    %
    % Manually defining Pamd or editing Pamd should be avoided unless you are familiar
    % with AMD or similar permutation algorithms for optimal LDL decomposition.
    % incorrectly defining/editing Pamd can severely slow down execution time or even break
    % the solver (or possibly unexpectedly *CRASH* MATLAB!).
    %
    % superADMM will throw a warning if the given Pamd could not be for this sparsity
    % pattern, as this can significantly slow down the solver.

    % (c) Peter Verheijen, 2026

    if ~issparse(P) || ~issparse(A)
        error("P and A must be both Sparse matrices");
    end
    nPrim = size(P,1);
    nDual = size(A,1);
    q = zeros(nPrim,1);
    l = zeros(nDual,1);
    u = zeros(nDual,1);
    opts.verbose = 0;
    opts.maxIter = 0;
    %opts.maxIter = 0 ensures only initialization/termination routine is executed
    [~,~,~,info] = superADMM(P, q, A, l, u, [], [], opts);
    Pamd = info.Pamd;
end