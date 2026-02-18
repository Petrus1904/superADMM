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
    A = np.array([[1, 1, 1, 1, 1],
                  [1, -1, 1, -1, 1],
                  [-1, -1, -1, -1, -1],
                  [1, 2, 3, 4, 5]])
    b = np.array([10, 3, -10, 20])
    A = np.vstack((A, np.eye(5)))
    l = np.hstack((-np.inf*np.ones(4), np.array([0, 0, 0, 0, 0])))
    u = np.hstack((b, np.array([10, 10, 10, 10, 10])))

    # opts = superADMM.getOptimOpts()
    opts = {"verbose": 2}
    #solve dense
    [x, y, eflag, info] = superADMM.superADMM(Q, c, A, l, u, options = opts)
    print(x)
    print(y)
    print(info["runtime"])
    #solve 
    A = csc_matrix(A)
    Q = csc_matrix(Q)
    [x1, y1, eflag1, info1] = superADMM.superADMM(Q, c, A, l, u, options = opts)
    print(x1)
    print(y1)
    print(info1["runtime"])