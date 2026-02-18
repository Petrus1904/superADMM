# SPARSE Bunch-Kaufman LDL' Decomposition
The files in this folder contain a (mostly) homemade LDL decomposition based on the Bunch-Kaufman (BK) algorithm.
Instead of the common LDL' decomposition, where D is strictly diagonal, the BK decomposition includes possible 2x2 diagonal blocks in D. As a result, the decomposition is often numerically more stable, for an (ideally) marginal computational increase, as the 2x2 blocks can be efficiently inverted using Cramers' rule.
In practice, a version BK-LDL decomp is employed in the MA57 algorithm, which is used by intel-MKL and MATLAB (MATLAB only for solving indefinite non-singular sparse matrices).

The execution of the decomposition and symbolic factoring is build/heavily inspired from the LDL-implementation made by Timothy A. Davis, which is also the LDL method currently employed in superADMM. However, as of this point, the code is substantially different.

The BK-decomposition made here is completely handwritten by me (Peter). Performance wise, it does show that MA57 remains more accurate, and definitely faster for large scale matrices. For smaller matrices (say less than 1000x1000) this BK-LDL decomp could execute faster. Additionally, BK-LDL produces more accurate solutions that naive LDL.

Despite the fact that this BK-sparse LDL decomposition does work*, it is currently not employed in the superADMM solver, as the resulting decomposition is substantially more dense, which slows down the over al solver time. In fact, I made a couple observations:
 - BK needs to recompute permutation & symbolic every iteration, as the permutation depends now on the value scale too.
 - BK is in general slower than LDL because of 2x2 pivots, even if the number of nonzeros are comparable.
 - However, using BK instead of LDL reduces the iteration count by 40%!
 - BK can solve problems that are numerically terribly scaled (very nice!)
 - In general --> BK scales terrible with problem size, and is as of now never faster than standard LDL.

Some comparison results, note that the sparse matrices are generated in a similar structure as commonly seen in superADMM
| size n | nnz    | tMA57    | tLDL      | tBKLDL     | epsMA57 | epsLDL  | epsBKLDL |
| ------ | ------ | -------- | --------- | ---------- | ------- | ------- | -------- |
| 30     | 84     | 0.053ms  | 0.028ms   | 0.035ms    | 2.2e-15 | 3.4e-08 | 8.9e-16  |
| 150    | 476    | 0.309ms  | 0.200ms   | 0.168ms    | 4.3e-15 | 5.7e-08 | 5.3e-15  |
| 750    | 2338   | 0.924ms  | 1.046ms   | 0.473ms    | 2.9e-15 | 7.0e-08 | 4.2e-14  |
| 1200   | 4672   | 2.809ms  | 1.204ms   | 3.581ms    | 5.8e-15 | 1.4e-07 | 5.7e-14  |
| 1500   | 6146   | 1.871ms  | 2.136ms   | 19.566ms   | 3.8e-15 | 1.1e-07 | 1.3e-13  |
| 1800   | 8880   | 4.046ms  | 2.221ms   | 35.592ms   | 6.2e-15 | 1.1e-07 | 1.2e-13  |
| 3150   | 25330  | 9.632ms  | 8.201ms   | 289.726ms  | 1.5e-14 | 2.2e-07 | 7.4e-12  |
| 5250   | 42630  | 19.309ms | 17.343ms  | 2387.984ms | 1.0e-14 | 2.8e-07 | 2.3e-11  |
| 6750   | 33928  | 10.283ms | 8.854ms   | 9226.626ms | 6.2e-15 | 1.5e-07 | 4.3e-11  |

From inspecting the decomposed L matrix, I noticed that this BK algorithm leaves a significant chuck of non-zero (but *very* small) elements, which MA57 does not have. This scales easily to over 10x the number of nonzeros and thus significantly slow down decomposition and inversion.
    

*I tested the code against expected inputs, where it works as intended (see table above). I cannot guarantee it also works under unexpected inputs.
