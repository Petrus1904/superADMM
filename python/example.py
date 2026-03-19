"""
@author: Peter Verheijen

Example python script on superADMM usage.
"""
import numpy as np
from scipy.sparse import csc_matrix
import os, sys
import superADMM

if __name__ == '__main__':
    #sample code
    Q = 2 * np.eye(5)
    c = np.array([-2, -6, -8, -4, -10])
    A = np.array([[1, -1, 1, -1, 1],
                  [-1, -1, -1, -1, -1],
                  [1, 2, 3, 4, 5]])
    b = np.array([3, -10, 20])
    A = np.vstack((A, np.eye(5)))
    l = np.hstack((-np.inf*np.ones(3), np.array([0, 0, 0, 0, 0])))
    u = np.hstack((b, np.array([10, 10, 10, 10, 10])))

    x_true = [3.4, 3.6, 2.8, 0, 0.2]
    y_true = [0, 8.4, 3.6, 0, 0, 0, -2, 0]
    # opts = superADMM.getOptimOpts()
    opts = {"verbose": 0}
    #solve dense
    [x, y, eflag, info] = superADMM.superADMM(Q, c, A, l, u, options = opts)
    if(eflag == 1):
        print("superADMM_dense success!, took: %.3e seconds" % info["runtime"])
        print("  iterations: %d" % info["nIter"])
        print("  x_err:      %.3e" % np.linalg.norm(x_true-x))
        print("  y_err:      %.3e" % np.linalg.norm(y_true-y))
    else:
        print("superADMM_dense failed, with message: " + info["status"])
    #solve 
    A = csc_matrix(A)
    Q = csc_matrix(Q)
    [x1, y1, eflag1, info1] = superADMM.superADMM(Q, c, A, l, u, options = opts)
    if(eflag1 == 1):
        print("superADMM_sparse success!, took: %.3e seconds" % info1["runtime"])
        print("  iterations: %d" % info1["nIter"])
        print("  x_err:      %.3e" % np.linalg.norm(x_true-x1))
        print("  y_err:      %.3e" % np.linalg.norm(y_true-y1))
    else:
        print("superADMM_dense failed, with message: " + info1["status"])