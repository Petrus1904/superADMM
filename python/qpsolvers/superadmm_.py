"""Implementation of the superADMM solver in the QPbenchmark repository
https://github.com/qpsolvers/qpsolvers
"""

"""Solver interface for `superADMM <https://superadmm.org>, <https://github.com/Petrus1904/superADMM> `__.

SuperADMM is an Augmented Lagrangian, operator splitting ADMM solver.
Unlike standard ADMM methods, superADMM assigns an individual penalty to each constraint
and updates each penalty based on it being currently in the active or inactive sets.
This results in similar behaviour to the Method of Multipliers (see Bertsekas, 1976), where for
increasing penalty weights, the convergence becomes superlinear. Although superADMM is implemented
both for sparse and dense problems, we exclusively test the sparse implementation, as it is most often
faster. For full reference, see: 
``Verheijen, P.C.N., Goswami, D., and Lazar, M. (2025). SuperADMM: Solving Quadratic Programs Faster 
with Dynamic Weighting ADMM. arXiv: 2506.11608 [math.OC]. URL: https://arxiv.org/abs/2506.11608``

**Warm-start:** this solver interface supports warm starting
"""

import warnings
from typing import Optional, Union

import numpy as np
import scipy.sparse as spa
from scipy.sparse import csc_matrix

from ..conversions import ensure_sparse_matrices
from ..problem import Problem
from ..solution import Solution
import superADMM

def superadmm_solve_problem(
    problem: Problem,
    initvals: Optional[np.ndarray] = None,
    verbose: bool = False,
    **kwargs,
) -> Solution:
    """Solve a quadratic program using superADMM.

    Parameters
    ----------
    problem :
        Quadratic program to solve.
    initvals :
        Warm-start guess vector for the primal solution.
    verbose :
        Set to `True` to print out extra information.

    Returns
    -------
    :
        Solution returned by the solver.

    Raises
    ------

    Note
    ----

    Notes
    -----
    Keyword arguments are forwarded to superADMM 
    see help(superADMM.superADMM) for detailed list of kwargs and returns

    Lower values for absolute or relative tolerances yield more precise
    solutions at the cost of computation time. See *e.g.* [Caron2022]_ for an
    overview of solver tolerances.
    """
    P, q, G, h, A, b, lb, ub = problem.unpack()
    P, G, A = ensure_sparse_matrices("osqp", P, G, A)

    #Pack the constraints to the superADMM format, i.e.,
    # l_super <= A_super * x <= u_super
    A_super = None
    l_super = None
    u_super = None
    if G is not None and h is not None:
        A_super = G
        l_super = np.full(h.shape, -np.inf)
        u_super = h
    if A is not None and b is not None:
        A_super = A if A_super is None else spa.vstack([A_super, A], format="csc")
        l_super = b if l_super is None else  np.hstack([l_super, b])
        u_super = b if u_super is None else  np.hstack([u_super, b])
    if lb is not None or ub is not None:
        lb = lb if lb is not None else np.full(q.shape, -np.inf)
        ub = ub if ub is not None else np.full(q.shape, +np.inf)
        E = spa.eye(q.shape[0], format="csc")
        A_super = E  if A_super is None else spa.vstack([A_super, E], format="csc")
        l_super = lb if l_super is None else np.hstack([l_super, lb])
        u_super = ub if u_super is None else np.hstack([u_super, ub])

    opts = kwargs
    opts["verbose"] = int(verbose)
    [x,y,eflag,info] = superADMM.superADMM(P, q, A_super, l_super, u_super, x0=initvals, options=opts)

    solution = Solution(problem)
    solution.extras = {
        "info": info,
    }
    solution.found = (eflag == 1)
    if not solution.found:
        warnings.warn(f"SuperADMM exited with status '{info["status"]}'")
    solution.x = x
    m = G.shape[0] if G is not None else 0
    meq = A.shape[0] if A is not None else 0
    solution.z = y[:m] if G is not None else np.empty((0,))
    solution.y = y[m : m + meq] if A is not None else np.empty((0,))
    solution.z_box = (
        y[m + meq :]
        if lb is not None or ub is not None
        else np.empty((0,))
    )
    return solution


def superadmm_solve_qp(
    P: Union[np.ndarray, csc_matrix],
    q: np.ndarray,
    G: Optional[Union[np.ndarray, csc_matrix]] = None,
    h: Optional[np.ndarray] = None,
    A: Optional[Union[np.ndarray, csc_matrix]] = None,
    b: Optional[np.ndarray] = None,
    lb: Optional[np.ndarray] = None,
    ub: Optional[np.ndarray] = None,
    initvals: Optional[np.ndarray] = None,
    verbose: bool = False,
    **kwargs,
) -> Optional[np.ndarray]:
    problem = Problem(P, q, G, h, A, b, lb, ub)
    solution = superadmm_solve_problem(problem, initvals, verbose, **kwargs)
    return solution.x if solution.found else None
