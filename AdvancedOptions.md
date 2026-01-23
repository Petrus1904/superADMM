# Solver Settings
In this document, we elaborate on the solver settings and the advanced usage of superADMM
First, for accessing the quick explanation of the solvers' arguments and options etc, one can refer to the help files in the functions themselves, accessed as 
MATLAB: `help superADMM`
Python: `help(pysuperADMM.superADMM)`

| option | dtype | range | Default | usage |
| --- | --- | --- | --- | --- |
| verbose | Int | {0,1,2} | 1 |controls the solver print behaviour. 0 for silence, 1 for intermediate output, 2 shows data usage in KB and elaborate timing. Mind that the timing can be offset since measuring time can be time consuming. |
| maxIter | Int | > 0 | 500 |Maximum number of iterations before the solver terminates with exitflag 0. |
| sigma | Double |>= 0 | 1e-6| The weight corresponding to a relaxation weight on x. This robustifies the solver, but can slow down convergence. One can increase sigma if the solver returns exitflag 2 to see if it can improve the solution accuracy. |
| rho_0 | Double | > 0 | 1.0 |Initial penalty of the Augmented Lagrangian term. There is ADMM theory that picking the best rho_0 can accelerate convergence, however superADMM is less sensitive to this. |
| tau | Double | (0, 1) | 0.5 |every iteration, superADMM verifies how accurate the solution was of the previous step. Due to increasing weights, this can suffer from numerical instabilities. If this harms convergence, the large weights are reduced with factor tau. |
| alpha | Double | >= 1 | 500.0 |Alpha controls how fast the penalties in the augmented lagrangian increase or decrease. Small alphas ensure a controlled convergence, but is slower. 
| RBound | Double | >> 1 | 1.0e8 | The penalties in R (see paper) only increase up to RBound (or 1/RBound). To properly assert superlinear convergence, RBound >> max(P,A). |
| eps_abs & eps_rel | Double | >= 0 | 1e-8 (both) | eps_abs and eps_rel define the solution accuracy. The solver terminates if: $\|Ax^k-z^k\|\leq \epsilon_{abs}+\epsilon_{rel}\max(\|Ax^k\|, \|z^k\|)$ and $\|Px^k+Ay^k+q\|\leq \epsilon_{abs}+\epsilon_{rel}\max(\|Px^k\|, \|Ay^k\|, \|q\|)$. |
| eps_inf | Double | > 0 | 1e-8 | Controls the infeasibility detection. One could increase this if the solver does not terminate correctly to detect if infeasibility is an issue. |
| repInterval | Int | >= 0 | 10 | If verbose = 1, repInterval defines how many iterations should be between a printed update. |
| timeLimit | Double | - | 0 | In seconds, if timeLimit > 0, the solver terminates at the iteration the timelimit is exceeded. Use for time-critical situations |
| lowRankPer | Double | [0, 1] | 0.05 | Instead of recomputing a large LDL decomposition, if the number of updates in R are less than lowRankPer*nDual, superADMM performs a low-rank update only, which is more efficient. |
| Pamd | Int* | - | - | See below |

# Advanced usage of superADMM
Every superADMM iteration, the solver computes the solution of the following linear system of equations: 
```math
\begin{bmatrix}P+\sigma I & A^T \\ A & -(R^k){-1}\end{bmatrix}\begin{bmatrix}x \\ \nu \end{bmatrix} = \begin{bmatrix}\sigma x^k-q \\  z^k-(R^k)^{-1}y^k\end{bmatrix}
```
For sparse $P$ and $A$, this is solved by computing the sparse LDL decomposition and then two sparse triangular solves. As the LDL decomposition of sparse matrices can become dense (and thus very slow and memory expensive), superADMM permutes the rows and the columns such that the decomposition is as sparse as possible. This permutation algorithm is known as Approximate Minimum Degree ordering or AMD.
Depending on the structure of $P$ and $A$, AMD can be quite expensive to compute. Fortunately, since the sparse structure remains static throughout the superADMM iterations, it is only performed once.

There are many situations where one has to solve the same problem over and over again, only changing the values of the problem, not the sparse structure (think of MPC, SQP, etc). Then, superADMM allows to precompute the permutation vector Pamd only once, and then use that for every subsequent superADMM call, thus saving some computation time. We include two ways to do this, which we illustrate next.
## Matlab:
```
solveropts.verbose = 0 %prevent flooding the command window
solveropts.Pamd = superADMM_preCompute(P,A);
for i = 1:maxIter
  %modify q, l, u and even the values in P and A here.
  [x,y,eflag,info] = superADMM(P,q,A,l,u, [], [], solveropts);
end
%OR:
solveropts.verbose = 0 %prevent flooding the command window
for i = 1:maxIter
  %modify q, l, u and even the values in P and A here.
  [x,y,eflag,info] = superADMM(P,q,A,l,u, [], [], solveropts);
  solveropts.Pamd = info.Pamd;
end
```
## Python:
```
import pysuperADMM
solveropts = { "verbose": 0}
solveropts["Pamd"] = pysuperADMM.preCompute(P,A)
for i in range(0, maxIter-1)
  #modify q, l, u and even the values in P and A here.
  [x,y,eflag,info] = pysuperADMM.superADMM(P, q, A, l, u, options=solveropts)
#OR:
solveropts = { "verbose": 0}
for i in range(0, maxIter-1)
  #modify q, l, u and even the values in P and A here.
  [x,y,eflag,info] = pysuperADMM.superADMM(P, q, A, l, u, options=solveropts)
  solveropts["Pamd"] = info["Pamd"]
```
Please note that since Pamd is openly accessible, possible augmentations or alternative AMD algorithms can compute (perhaps) better permutation vectors. We thus allow to change this, however, superADMM will throw a warning. 
Additionally, incorrectly changing Pamd can result in a significant increase in execution time, or even break the solver entirely. Therefore, we advice to keep Pamd as is.
Since you have now read this warning, if you still seek to change Pamd, please note that the length is `nPrim+nDual+1`. The last entry defines the number of non-zeros in the permuted L, which verifies if a correct permutation is used. As such, if you put in the correct number of nonzeros there, superADMM will not throw a warning.

