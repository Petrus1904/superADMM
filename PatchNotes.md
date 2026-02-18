# Patch Notes
## V0.7.0
- Relative termination criteria are now included. The solver now terminates if $\|Ax^k-z^k\|\leq \epsilon_{abs}+\epsilon_{rel}\max(\|Ax^k\|, \|z^k\|)$ and $\|Px^k+Ay^k+q\|\leq \epsilon_{abs}+\epsilon_{rel}\max(\|Px^k\|, \|Ay^k\|, \|q\|)$. `opts.eps_rel` is added as a solver options to tweak this bound. Set to zero to disable relative termination, default `1e-8`.
- As a result of the change above, the info package at the end of the solver also includes the considered termination criteria, to avoid confusion.
- The permutation vector obtained by performing AMD ordering on the KKT matrix can now be re-used between runs. Please read `AdvancedOptions.md` for usage.
- Adding on the previous note, the permutation vector can be pre-obtained using `superADMM_preCompute.m` and `superADMM.preCompute()` in MATLAB and Python, respectively.
- The Python interface is rewritten in Cython and comes with a pip-installable module. The idea is to make the install routine for Python also `pip install superADMM` for near future versions.
- For Python users, the info and options structs are now dictionaries. This allows for flexible usage, where the user only has to define the options that it wants to change. This does however alter the syntax.

## V0.6.0
- Sparse solver can now update the LDL matrices with low-rank updates. This is enabled if the updates in R are less than 5% of the number of constraints (which can be adjusted).
- Dense solver updates P+A'RA with low rank updates, this is enabled if the updates are less than 50% of the number of constraints.
- New option: `opts.lowRankPer`, which thresholds when low rank updates are favored over full decomposition (number of updates is less than `opts.lowRankPer*m`, where `m` is the number of constraints). Mind that this differs per problem type and can result in substantial computational improvements for problems that take more than ~8 iterations to complete.
- A^T is no longer explicitly stored in the sparse solver, reducing memory usage.
- Function headers of `superADMMsolverSparse` and `superADMMsolverDense` now include argument documentation, nice for `c` users.
- Both solvers now display their memory usage of both the problem size and workspace memory in KB. Useful for future embedded implementations with limited memory sizes. This is enabled by setting `opts.verbose = 2`.

## V0.5.0
Initial Release
