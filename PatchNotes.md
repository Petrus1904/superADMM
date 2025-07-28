# Patch Notes
## V0.6.0
- Sparse solver can now update the LDL matrices with low-rank updates. This is enabled if the updates in R are less than 5% of the number of constraints (which can be adjusted).
- Dense solver updates P+A'RA with low rank updates, this is enabled if the updates are less than 50%.
- New option: 0<=`opts.lowRankPer`<=1, which controls when low rank updates are favored over full decomposition. Mind that this differs per problem type and can result in substantial computational improvements for problems that take more than ~8 iterations to complete.
- A' is no longer explicitly stored in the sparse solver, reducing memory usage.
- Function headers of `superADMMsolverSparse` and `superADMMsolverDense` now include argument documentation, nice for `c` users.
- Both solvers now display their memory usage of both the problem size and workspace memory in KB. Useful for future embedded implementations with limited memory sizes. This is enabled by setting `opts.verbose = 2`.

## V0.5.0
Initial Release
