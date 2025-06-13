# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 10:24:37 2024

@author: Peter Verheijen
"""
import ctypes
import numpy as np
from scipy.sparse import csc_matrix

lib = ctypes.CDLL('superADMM.dll')

class ADMMopts(ctypes.Structure):
    _fields_ = [("verbose", ctypes.c_int32),
                ("maxIter", ctypes.c_int32),
                ("sigma", ctypes.c_double),
                ("rho_0", ctypes.c_double),
                ("tau", ctypes.c_double),
                ("alpha", ctypes.c_double),
                ("RBound", ctypes.c_double),
                ("eps_abs", ctypes.c_double),
                ("eps_inf", ctypes.c_double),
                ("repInterval", ctypes.c_int32),
                ("timeLimit", ctypes.c_double),
                ("lowRankPer", ctypes.c_double)]
    
class ADMMinfo(ctypes.Structure):
    _fields_ = [("nIter", ctypes.c_int32),
                ("rPrim", ctypes.c_double),
                ("rDual", ctypes.c_double),
                ("runtime", ctypes.c_double),
                ("objVal", ctypes.c_double),
                ("status", ctypes.c_char_p)]
    
    
lib.superADMMsolverDense.argtypes=[ctypes.POINTER(ctypes.c_double), #P
                                   ctypes.POINTER(ctypes.c_double), #q
                                   ctypes.POINTER(ctypes.c_double), #A
                                   ctypes.POINTER(ctypes.c_double), #l
                                   ctypes.POINTER(ctypes.c_double), #u
                                   ctypes.POINTER(ctypes.c_double), #x
                                   ctypes.POINTER(ctypes.c_double), #y
                                   ctypes.c_int32, ctypes.c_int32, ADMMopts, ctypes.POINTER(ADMMinfo)] #nPrim, nDual, opts
lib.superADMMsolverDense.restype = ctypes.c_int32

lib.superADMMsolverSparse.argtypes=[ctypes.POINTER(ctypes.c_double), #Pdata
                                    ctypes.POINTER(ctypes.c_int32), #Pcol
                                    ctypes.POINTER(ctypes.c_int32), #Pidx
                                    ctypes.c_int32, #Pnnz
                                    ctypes.POINTER(ctypes.c_double), #q
                                    ctypes.POINTER(ctypes.c_double), #Adata
                                    ctypes.POINTER(ctypes.c_int32), #Acol
                                    ctypes.POINTER(ctypes.c_int32), #Aidx
                                    ctypes.c_int32, #Annz
                                    ctypes.POINTER(ctypes.c_double), #l
                                    ctypes.POINTER(ctypes.c_double), #u
                                    ctypes.POINTER(ctypes.c_double), #x
                                    ctypes.POINTER(ctypes.c_double), #y
                                    ctypes.c_int32, ctypes.c_int32, #nPrim, #nDual
                                    ADMMopts, ctypes.POINTER(ADMMinfo)] #opts, #info
lib.superADMMsolverSparse.restype = ctypes.c_int32

def getOptimOpts():
    '''ADMMoptions = getOptimOpts()
    
    returns the default settings for the superADMM solver.
    
    =ADMMopts (struct) parameters=========================================
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
         absolute convergence bound
     - eps_inf (float):   
         infeasibility bounds, > 0
     - repInterval (int): 
         (only if verbose = 1), number of iterations between reporting the status. Default 10.
     - timeLimit (float): 
         Time limit in seconds, use timeLimit <= 0 for unlimited, default 0.
    '''
    
    options = ADMMopts(verbose=0,    maxIter=500,
                       sigma=1e-6,   rho_0=1.0,
                       tau=0.5,      alpha=500,
                       RBound=1e8,   eps_abs=1e-8,
                       eps_inf=1e-8, repInterval=10, timeLimit=0, lowRankPer=0.0)
    return options

def superADMM(P, q, A, l, u, x0 = None, y0 = None, ADMMoptions = None):
    '''[x, y, eflag, info] = superADMM(P, q, A, l, u, x0, y0, ADMMoptions)
    
    QP solver for problems of structure
     min_x 0.5x'*P*x + x'*q
     s.t.  l<=A*x<=u
    
    Uses an ADMM adaptation that considers super linear convergence.
    If your problem has constraints in the M*x <= g; E*x = b format, you can
    rewrite it as:
    [-inf*ones(length(gam),1); b] <= [M;E]*x <= [gam;b];
    
    =Arguments===============================================================
     P (ndarray or csc_matrix): 
         (n x n) Positive semi-definite quadratic cost function matrix.
     q (ndarray):
         (n x 1) Vector with length equal to the rows and columns of P.
     A (ndarray or csc_matrix): 
         (m x n) Constraint matrix.
     l (ndarray):               
         (m x 1) Vector of lower bounds.
     u (ndarray):               
         (m x 1) Vector of upper bounds.
     x0 (ndarray, optional):    
         (n x 1) The initial guess for the solution x
     y0 (ndarray, optional):    
         (m x 1) The initial guess for the dual variable y
     ADMMoptions (ADMMopts, optional): 
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
          * -2: Problem infeasible
          * -3: Time limit exceeded
          * -4: Problem non-convex
      info (obj): 
          * nIter: number of iterations
          * rPrim: primal convergence error
          * rDual: dual convergence error
          * runtime: time elapsed of the C-solver, python overhead not included
          * objVal: the objective value
          * status: the solver status on terminations
    
    Use getOptimOpts() to get the ADMMopts object with the default initialization.
    
    Note that both P and A have to be of the same type. If both are a numpy array,
    a Dense solver will be used. If both are csc_matrix, a sparse solver is used.
    '''

    #assert things
    assert type(P) == np.ndarray or type(P) == csc_matrix, "P must be a numpy array or csc_matrix"
    assert type(A) == np.ndarray or type(A) == csc_matrix, "A must be a numpy array or csc_matrix"
    assert type(q) == np.ndarray, "q must be a numpy array"
    assert type(l) == np.ndarray, "l must be a numpy array"
    assert type(u) == np.ndarray, "u must be a numpy array"
    
    assert np.shape(P)[0] == np.shape(P)[1], "P must be square"
    assert np.shape(A)[1] == np.shape(P)[0], "A must have the same number of columns as rows in P"
    
    nPrim = np.shape(P)[0]
    nDual = np.shape(A)[0]
    assert len(q) == nPrim, "q must have the same number of elements as rows in P"
    assert len(l) == nDual, "l must have the same number of elements as rows in A"
    assert len(u) == nDual, "u must have the same number of elements as rows in A"
    
    if(ADMMoptions is None):
        ADMMoptions = getOptimOpts()
    else:
        if(not type(ADMMoptions) == ADMMopts):
            raise TypeError("ADMMoptions must be of type ADMMopts. Use getOptimOpts() to obtain the default settings")
        assert ADMMoptions.tau > 0 and ADMMoptions.tau < 1, "tau out of range, choose 0<tau<1"
        assert ADMMoptions.alpha >= 1, "alpha must be greater or equal to 1"
        assert ADMMoptions.sigma >= 0, "sigma must be nonnegative"
        assert ADMMoptions.RBound > 1, "RBound must be bigger than 1"
        assert ADMMoptions.rho_0 > 0, "rho_0 must be positive"
        
    if(x0 is None):
        x = np.zeros(nPrim, dtype=np.double)
    else:
        assert type(x0) == np.ndarray, "x0 must be numpy array"
        assert len(x0) == nPrim, "x0 must have the same number of elements as rows in P"
        x = np.array(x0, dtype=np.double) #copy to prevent confusion
    if(y0 is None):
        y = np.zeros(nDual, dtype=np.double)
    else:
        assert type(y0) == np.ndarray, "y0 must be numpy array"
        assert len(y0) == nDual, "x0 must have the same number of elements as rows in A"
        y = np.array(y0, dtype = np.double) #copy to prevent confusion
        
    #get pointers
    qt = np.array(q, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    lt = np.array(l, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    ut = np.array(u, dtype=np.double).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    xt = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    yt = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    
    info = ADMMinfo(nIter=0, rPrim=0, rDual=0, runtime=0, objVal=0)
    infopt = ctypes.byref(info)
        
    if(type(P) == np.ndarray and type(A) == np.ndarray):
        P = np.array(P, dtype = np.double, order='F') #convert to column major
        A = np.array(A, dtype = np.double, order='F')
        Pt = P.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        At = A.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        eflag = lib.superADMMsolverDense(Pt, qt, At, lt, ut, xt, yt, nPrim, nDual, ADMMoptions, infopt)
        
    elif(type(P) == csc_matrix and type(A) == csc_matrix):
        Ptdata = P.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        Ptcol = P.indptr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        Ptind = P.indices.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        Ptnnz = P.nnz
        Atdata = A.data.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        Atcol = A.indptr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        Atind = A.indices.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
        Atnnz = A.nnz
        eflag = lib.superADMMsolverSparse(Ptdata, Ptcol, Ptind, Ptnnz, qt,
                                         Atdata, Atcol, Atind, Atnnz, lt,
                                         ut, xt, yt, nPrim, nDual, ADMMoptions, infopt)
    else:
        raise Exception("Both P and A must be either dense or sparse CSC matrices")
        
    
    return x, y, eflag, info


if __name__ == '__main__':
    #sample code
    Q = 2 * np.eye(5)
    c = np.array([-2, -6, -8, -4, -10])
    A = np.array([[1, 1, 1, 1, 1],
                  [1, -1, 1, -1, 1],
                  [-1, -1, -1, -1, -1],
                  [1, 2, 3, 4, 5]])
    b = np.array([10, 3, -10, 20])
    A = np.vstack((A, np.eye(5)))
    l = np.hstack((-np.inf*np.ones(4), np.array([0, 0, 0, 0, 0])))
    u = np.hstack((b, np.array([10, 10, 10, 10, 10])))
    
    opts = getOptimOpts()
    opts.sigma = 1e-6
    # mind that printing to python can be really slow
    opts.verbose = 1
    opts.alpha = 500
    #solve dense
    [x, y, eflag, info] = superADMM(Q, c, A, l, u, ADMMoptions = opts)
    print(x)
    print(y)
    print(info.runtime)
    #solve sparse
    A = csc_matrix(A)
    Q = csc_matrix(Q)
    [x1, y1, eflag1, info1] = superADMM(Q, c, A, l, u, ADMMoptions = opts)
    print(x)
    print(y)
    print(info1.runtime)
    
    #x = [3.4, 3.6, 2.8, 0, 0.2]
                  
    
    




