# distutils: language = c
# cython: boundscheck=False, wraparound=False, initializedcheck=False
# cython: cdivision=True

import numpy as np
import scipy.sparse as sp
cimport numpy as np

np.import_array()   #<-- Required!

# C declarations
cdef extern from "superADMM.h":

    ctypedef struct ADMMopts:
        int verbose
        int maxIter
        double sigma
        double rho_0
        double tau
        double alpha
        double RBound
        double eps_abs
        double eps_rel
        double eps_inf
        int repInterval
        double timeLimit
        double lowRankPer
        int* Pamd

    ctypedef struct ADMMinfo:
        int nIter
        double prim_res
        double dual_res
        double prim_tol
        double dual_tol
        double runtime
        double objVal
        char* status
        int* Pamd

    int superADMMsolverDense(
        const double* P, 
        const double* q,
        const double* A,
        const double* l,
        const double* u,
        double* x,
        double* y,
        int nPrim,
        int nDual,
        ADMMopts opts,
        ADMMinfo* info
    )

    int superADMMsolverSparse(
        double* Pdata, 
        int* Prowptr, 
        int* Pcolidx, 
        const double* q,
        double* Adata,
        int* Arowptr,
        int* Acolidx,
        const double* l,
        const double* u,
        double* x,
        double* y,
        int nPrim,
        int nDual,
        ADMMopts opts,
        ADMMinfo* info
    )

cdef _superADMM_dense(np.ndarray[np.float64_t, ndim=2, mode="fortran"] P,
                      np.ndarray[np.float64_t, ndim=1] q,
                      np.ndarray[np.float64_t, ndim=2, mode="fortran"] A,
                      np.ndarray[np.float64_t, ndim=1] l,
                      np.ndarray[np.float64_t, ndim=1] u,
                      np.ndarray[np.float64_t, ndim=1] x,
                      np.ndarray[np.float64_t, ndim=1] y,
                      ADMMopts c_opts):
    # all error checking performed by calling function
    # pointers
    cdef int nPrim = int(P.shape[0])
    cdef int nDual = int(A.shape[0])
    cdef const double* q_ptr = <double*> q.data
    cdef const double* l_ptr = <double*> l.data
    cdef const double* u_ptr = <double*> u.data
    cdef double* x_ptr = <double*>x.data
    cdef double* y_ptr = <double*>y.data

    cdef ADMMinfo out_info
    cdef int eflag
    cdef const double* P_ptr = <double*> P.data
    cdef const double* A_ptr = <double*> A.data
    eflag = superADMMsolverDense(P_ptr, q_ptr, A_ptr, l_ptr, u_ptr,
                                    x_ptr, y_ptr, nPrim, nDual, c_opts, &out_info)

    info = {
        "nIter": out_info.nIter,
        "prim_res": out_info.prim_res,
        "dual_res": out_info.dual_res,
        "prim_tol": out_info.prim_tol,
        "dual_tol": out_info.dual_tol,
        "runtime": out_info.runtime,
        "objVal": out_info.objVal,
        "status": out_info.status,
        "Pamd": None
    }
    return x, y, eflag, info

cdef _superADMM_sparse(object P,
                      np.ndarray[np.float64_t, ndim=1] q,
                      object A,
                      np.ndarray[np.float64_t, ndim=1] l,
                      np.ndarray[np.float64_t, ndim=1] u,
                      np.ndarray[np.float64_t, ndim=1] x,
                      np.ndarray[np.float64_t, ndim=1] y,
                      ADMMopts c_opts, np.ndarray[np.int32_t, ndim=1] pyPamd):
    # all error checking performed by calling function
    # pointers
    cdef int nPrim = int(P.shape[0])
    cdef int nDual = int(A.shape[0])
    cdef const double* q_ptr = <double*> q.data
    cdef const double* l_ptr = <double*> l.data
    cdef const double* u_ptr = <double*> u.data
    cdef double* x_ptr = <double*>x.data
    cdef double* y_ptr = <double*>y.data

    cdef ADMMinfo out_info
    cdef int eflag
    cdef np.ndarray[np.float64_t, ndim=1] Pdata = P.data
    cdef np.ndarray[np.int32_t, ndim=1] Pind = P.indices
    cdef np.ndarray[np.int32_t, ndim=1] Pptr = P.indptr

    cdef double* Pdata_ptr = <double*> Pdata.data
    cdef int* Pind_ptr = <int*> Pind.data
    cdef int* Pptr_ptr = <int*> Pptr.data

    cdef np.ndarray[np.float64_t, ndim=1] Adata = A.data
    cdef np.ndarray[np.int32_t, ndim=1] Aind = A.indices
    cdef np.ndarray[np.int32_t, ndim=1] Aptr = A.indptr

    cdef double* Adata_ptr = <double*> Adata.data
    cdef int* Aind_ptr = <int*> Aind.data
    cdef int* Aptr_ptr = <int*> Aptr.data
    eflag = superADMMsolverSparse(Pdata_ptr, Pptr_ptr, Pind_ptr, q_ptr,
                                    Adata_ptr, Aptr_ptr, Aind_ptr, l_ptr, u_ptr,
                                    x_ptr, y_ptr, nPrim, nDual, c_opts, &out_info)

    info = {
        "nIter": out_info.nIter,
        "prim_res": out_info.prim_res,
        "dual_res": out_info.dual_res,
        "prim_tol": out_info.prim_tol,
        "dual_tol": out_info.dual_tol,
        "runtime": out_info.runtime,
        "objVal": out_info.objVal,
        "status": out_info.status,
        "Pamd": pyPamd
    }
    return x, y, eflag, info

#Python interface
def _fixInputData(x, vname, ndims=1, dattype=np.float64):
    if hasattr(x, "detach"):
        x = x.detach().cpu().numpy() #extract numpy array
    
    if np.isscalar(x) and ndims==1:
        return np.array([x], dtype=dattype)
    if np.isscalar(x) and ndims==2:
        return np.array([[x]], dtype=dattype)

    x = np.asarray(x, dtype=dattype)
    if x.ndim < ndims:
        #try brute forcing it
        x = np.array([x])
    if x.ndim > ndims:
        #try squeezing it
        x = np.squeeze(x)
    
    if x.ndim != ndims:
        raise IndexError(vname + " is of incorrect dimension: expected " + str(ndims) + " got " + str(x.ndim))

    return x    
        

def superADMM(P, q, A, l, u, x0=None, y0=None, options=None):
    '''[x, y, eflag, info] = superADMM(P, q, A, l, u, x0, y0, options)
    
    QP solver for problems of structure
     min_x 0.5x'*P*x + x'*q
     s.t.  l<=A*x<=u
    
    Uses an ADMM adaptation that considers super linear convergence.
    If your problem has constraints in the M*x <= g; E*x = b format, you can
    rewrite it as:
    A = np.vstack((M,E))
    l = np.hstack((-np.inf*np.ones(len(gam)), b))
    u = np.hstack((gam,b))
    
    =Arguments===============================================================
     P (ndarray or csc_matrix): 
         (n x n) Positive semi-definite quadratic cost function matrix.
     q (ndarray):
         (n) or (n x 1) Vector with length equal to the rows and columns of P.
     A (ndarray or csc_matrix): 
         (m x n) Constraint matrix.
     l (ndarray):               
         (m) or (m x 1) Vector of lower bounds.
     u (ndarray):               
         (m) or (m x 1)  Vector of upper bounds.
     x0 (ndarray, optional):    
         (n) or (n x 1) The initial guess for the solution x
     y0 (ndarray, optional):    
         (m) or (m x 1)  The initial guess for the dual variable y
     options (dict, optional): 
         Hyperparameters, see getOptimOpts()
    
    =Returns=================================================================
      x (nparray): 
          The solution vector
      y (nparray): 
          The dual variable
      eflag (int): 
          * 1: Solved Succesfully
          * 2: Requested tolerance cannot be achieved
          * 0: Maximum iterations reached
          * -1: Cholesky solver failed
          * -2: Problem infeasible or unbounded
          * -3: Time limit exceeded
          * -4: problem non-convex
      info (dict): 
          * nIter: number of iterations
          * prim_res: primal solution residual
          * dual_res: dual solution residual
          * prim_tol: primal convergence tolerance
          * dual_tol: dual convergence tolerance 
          * runtime: time elapsed of the C-solver, python overhead not included
          * objVal: the objective value
          * status: the solver status on terminations
    
    Use getOptimOpts() to get the ADMMopts object with the default initialization.
    
    Note that both P and A have to be of the same type. If both are a numpy array,
    a Dense solver will be used. If both are csc_matrix, a sparse solver is used.
    '''

    #assert types
    if not (sp.isspmatrix_csc(P) and sp.isspmatrix_csc(A)):
        P = _fixInputData(P, "P", ndims=2)
        A = _fixInputData(A, "A", ndims=2)
    q = _fixInputData(q, "q")
    l = _fixInputData(l, "l")
    u = _fixInputData(u, "u")
    if x0 is not None:
        x0 = _fixInputData(x0, "x0")
    if y0 is not None:
        y0 = _fixInputData(y0, "y0")


    if not (sp.isspmatrix_csc(P) or isinstance(P, np.ndarray)):
        raise TypeError("P must be a scipy.sparse.csc_matrix or numpy.ndarray")
    if not (sp.isspmatrix_csc(A) or isinstance(A, np.ndarray)):
        raise TypeError("A must be a scipy.sparse.csc_matrix or numpy.ndarray")
    if not isinstance(q, np.ndarray):
        raise TypeError("q must be a numpy.ndarray")
    if not isinstance(l, np.ndarray):
        raise TypeError("l must be a numpy.ndarray")
    if not isinstance(u, np.ndarray):
        raise TypeError("u must be a numpy.ndarray")
    if x0 is not None and not isinstance(x0, np.ndarray):
        raise TypeError("x0 must be a numpy.ndarray")
    if y0 is not None and not isinstance(y0, np.ndarray):
        raise TypeError("y0 must be a numpy.ndarray")

    if A.ndim != 2:
        raise IndexError("A must be a 2 dimensional matrix")
    if P.ndim != 2:
        raise IndexError("P must be a 2 dimensional matrix")

    cdef int nPrim = int(P.shape[0])
    cdef int nDual = int(A.shape[0])

    q = np.asarray(q, dtype=np.float64)
    if q.ndim != 1:
        q = np.squeeze(q)

    l = np.asarray(l, dtype=np.float64)
    if l.ndim != 1:
        l = np.squeeze(l)

    u = np.asarray(u, dtype=np.float64)
    if u.ndim != 1:
        u = np.squeeze(u)

    if x0 is not None:
        x0 = np.asarray(x0, dtype=np.float64)
        if x0.ndim != 1:
            x0 = np.squeeze(x0)
            
    if y0 is not None:
        y0 = np.asarray(y0, dtype=np.float64)
        if y0.ndim != 1:
            y0 = np.squeeze(y0)

    if P.shape[1] != nPrim:
        raise IndexError("P must be a square matrix")
    if A.shape[1] != nPrim:
        raise IndexError("Number of columns in A must match the number of columns in P")
    if len(q) != nPrim:
        raise IndexError("q must have the same number of elements as rows in P")
    if len(u) != nDual:
        raise IndexError("u must have the same number of elements as rows in A")
    if len(l) != nDual:
        raise IndexError("l must have the same number of elements as rows in A")


    cdef np.ndarray[np.int32_t, ndim=1] pyPamd = np.zeros(nPrim+nDual+1, dtype=np.int32)
    pyPamd[0] = -1 #this flags superADMM to perform AMD
    #ADMM options
    cdef ADMMopts c_opts 
    c_opts.verbose      = 1
    c_opts.maxIter      = 500
    c_opts.sigma        = 1e-6
    c_opts.rho_0        = 1.0
    c_opts.tau          = 0.5
    c_opts.alpha        = 500.0
    c_opts.RBound       = 1.0e8
    c_opts.eps_abs      = 1e-8
    c_opts.eps_rel      = 1e-8
    c_opts.eps_inf      = 1e-8
    c_opts.repInterval  = 10
    c_opts.timeLimit    = 0.0
    c_opts.lowRankPer   = 0.05

    if options is not None:
        #verify for each option if it is included, overwrite, we check later
        if "verbose" in options:
            c_opts.verbose  = int(options["verbose"])
        if "maxIter" in options:
            c_opts.maxIter  = int(options["maxIter"])
        if "sigma" in options:
            c_opts.sigma    = <double>options["sigma"]
        if "rho_0" in options:
            c_opts.rho_0    = <double>options["rho_0"]
        if "tau" in options:
            c_opts.tau      = <double>options["tau"]
        if "alpha" in options:
            c_opts.alpha    = <double>options["alpha"]
        if "RBound" in options:
            c_opts.RBound   = <double>options["RBound"]
        if "eps_abs" in options:
            c_opts.eps_abs  = <double>options["eps_abs"]
        if "eps_rel" in options:
            c_opts.eps_rel  = <double>options["eps_rel"]
        if "eps_inf" in options:
            c_opts.eps_inf  = <double>options["eps_inf"]
        if "repInterval" in options:
            c_opts.repInterval = int(options["repInterval"])
        if "timeLimit" in options:
            c_opts.timeLimit = <double>options["timeLimit"]
        if "lowRankPer" in options:
            c_opts.lowRankPer = <double>options["lowRankPer"]
        if "Pamd" in options:
            Pin = options["Pamd"]
            if len(Pin) != nDual+nPrim+1:
                raise IndexError("Permutation array Pamd is of incorrect length")
            #copy
            pyPamd = np.copy(Pin)

    c_opts.Pamd = <int*> pyPamd.data

    #check options
    if c_opts.verbose < 0 or c_opts.verbose > 2:
        raise ValueError("options.verbose must be either 0, 1, or 2")
    if c_opts.alpha < 1:
        raise ValueError("options.alpha must be larger or equal to 1")
    if c_opts.rho_0 <=0:
        raise ValueError("options.rho_0 must be positive non-zero")
    if c_opts.tau > 1 or c_opts.tau <= 0:
        raise ValueError("options.tau must be (0, 1]")
    if c_opts.RBound < 1:
        raise ValueError("options.RBound must be greater than 1")
    if c_opts.eps_abs < 0 or c_opts.eps_rel < 0:
        raise ValueError("options.eps_abs and options.eps_rel must be greater or equal to zero")
    if c_opts.eps_inf <= 0:
        raise ValueError("options.eps_inf must be greater than zero")
    if c_opts.repInterval < 1:
        raise ValueError("options.repInterval must be greater than 1")
    if c_opts.lowRankPer > 1 or c_opts.lowRankPer < 0:
        raise ValueError("options.lowRankPer must be [0, 1]")

    cdef np.ndarray[np.float64_t, ndim=1] x = np.zeros(nPrim, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] y = np.zeros(nDual, dtype=np.float64)

    if x0 is not None:
        if len(x0) != nPrim:
            raise IndexError("x0 must have the same number of elements as rows in P")
        x = np.copy(x0)
    if y0 is not None:
        if len(y0) != nDual:
            raise IndexError("y0 must have the same number of elements as rows in A")
        y = np.copy(y0)
        
    cdef ADMMinfo out_info
    cdef int eflag
    if sp.isspmatrix_csc(P) and sp.isspmatrix_csc(A):
        return _superADMM_sparse(P, q, A, l, u, x, y, c_opts, pyPamd)
    elif isinstance(P, np.ndarray) and isinstance(A, np.ndarray):
        P = np.asarray(P, dtype=np.float64, order='F')
        A = np.asarray(A, dtype=np.float64, order='F')
        return _superADMM_dense(P, q, A, l, u, x, y, c_opts)
    else:
        raise TypeError("Both A and P must be of the same datatype")


def getOptimOpts():
    '''options = superADMM.getOptimOpts()
    
    returns the default settings for the superADMM solver as a dict.
    
    =options (dict) fields=========================================
     - verbose (int):    
         0 to silence (default), 1 to print intermediate output, 2 to show elaborate timing
     - maxIter (int):     
         default 500
     - sigma (float):     
         relaxation weight parameter on x, default 1e-6
     - rho_0 (float):     
         the initial rho value > 0, default 1
     - tau (float):       
         the exponential decrease rate on the bound < 1, default 0.5
     - alpha (float):     
         exponential increase rate on R, alpha => 1, default 500
     - RBound (float):    
         limits on r, default 1e8
     - eps_abs (float):   
         absolute convergence tolerance default 1e-8
     - eps_rel (float):   
         relative convergence tolerance default 1e-8
     - eps_inf (float):   
         infeasibility bounds, > 0
     - repInterval (int): 
         (only if verbose = 1), number of iterations between reporting the status. Default 10.
     - timeLimit (float): 
         Time limit in seconds, use timeLimit <= 0 for unlimited, default 0.
     - lowRankPer (float);
         Threshold for low rank updates, 0<=lowRankPer<=1, default 0.05.
    '''
    
    options = {"verbose"    : 1,
               "maxIter"    : 500,
               "sigma"      : 1e-6,
               "rho_0"      : 1.0,
               "tau"        : 0.5,
               "alpha"      : 500,
               "RBound"     : 1e8,
               "eps_abs"    : 1e-8,
               "eps_rel"    : 1e-8,
               "eps_inf"    : 1e-8,
               "repInterval": 10,
               "timeLimit"  : 0,
               "lowRankPer" : 0.05}
    return options

def preCompute(P, A):
    ''' Pamd = superADMM.preCompute(P,A)

    Advanced usage of the superADMM solver, allows for more efficient repetive 
    calls of the superADMM solver if the same problem sparsity structure in 
    P and A is guaranteed to remain the same. Otherwise, this function is not necessary.
    This function only works for CSC_matrix, not for any dense representation

    =Arguments===============================================================
     P (csc_matrix): 
         (n x n) Positive semi-definite quadratic cost function matrix.
     A (csc_matrix): 
         (m x n) Constraint matrix.

    =Returns=================================================================
     Pamd (ndarray, int32):
         Permutation vector corresponding to the sparsity in P and A.

    =Usage===================================================================
        #Assume P,q,A,l,u are already defined
        options = {"verbose":0}
        options["Pamd"] = superADMM.preCompute(P,A)
        for i in range(0,maxIter):
            #update q,l,u, even the values in P and A. only the locations of the
            #non-zeros in P and A must remain static.
            [x, y, eflag, info] = superADMM.superADMM(P,q,A,l,u,options=options)

    Please note that the same can be achieved with the following example, but 
    can be considered "less neat":
        #Assume P,q,A,l,u are already defined
        options = {"verbose":0}
        for i in range(0,maxIter):
            #update q,l,u, even the values in P and A. only the locations of the
            #non-zeros in P and A must remain static.
            [x, y, eflag, info] = superADMM.superADMM(P,q,A,l,u,options=options)
            options["Pamd"] = info["Pamd"]

    Manually defining Pamd or editing Pamd should be avoided unless you are familiar
    with the AMD or similar permutation algorithms for sparse LDL decomposition.
    incorrectly defining/editing Pamd can severely slow down execution time or even break
    the solver.
    
    superADMM will throw a warning if the given Pamd could not be for this sparsity
    pattern, as this can significantly slow down the solver.
    '''

    #assert types
    if not sp.isspmatrix_csc(P):
        raise TypeError("P must be a scipy.sparse.csc_matrix")
    if not sp.isspmatrix_csc(A):
        raise TypeError("A must be a scipy.sparse.csc_matrix")

    options = {"verbose":0,"maxIter":0} 
    cdef int nPrim = P.shape[0]
    cdef int nDual = A.shape[0]
    cdef np.ndarray[np.float64_t, ndim=1] q = np.zeros(nPrim, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] l = np.zeros(nDual, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] u = np.zeros(nDual, dtype=np.float64)
    #setting maxIter to zero ensures we only run the initialization/termination routine
    info = superADMM(P,q,A,l,u,options=options)[3] #only return the info object
    return info["Pamd"]
    
