/*
    Header file for the Bunch Kaufman sparse LDL decomposition and solver
    Created by Peter Verheijen, 2025.
*/
#define ADMMint long long int
#define ADMMfloat double
/*==============================================================================
bk_permute computes the Bunch-Kaufman row-column permutation, which changes
The existing permutation given on the input (as a result of symbolic decomp for 
instance), and returns a vector indicating when a 1x1 pivot must be done, and when a
2x2 pivot is expected.
================================================================================*/
ADMMint bk_permute
(
    const ADMMint n,            /*col/rows of A*/
    const ADMMfloat *Ax,    /*const nonzero values of A*/
    const ADMMint *Ai,          /*const row indices of A*/
    const ADMMint *Ap,          /*const column pointers of A*/
    ADMMint *P,                 /*in/out permutation vector*/
    ADMMint *Pinv,              /*in/out inverse of P*/
    ADMMint *b2                 /*output vector of size n+1 defining when 2x2 pivots are needed*/
);

/* ========================================================================== */
/* === bk_dsolve:  solve Dx=b ============================================== */
/* ========================================================================== */

void bk_dsolve
(
    ADMMint n,		    /* length of x */
    ADMMfloat *x,	/* b on the input, x on the output */
    ADMMfloat *D,	/* D vector, block diagonal */
    ADMMint *b2         /*when to do 1x1 of when 2x2 inverses */
);

/* ========================================================================== */
/* === bk_numeric: decompose A=LDL^T ======================================== */
/* ========================================================================== */

ADMMint bk_numeric	/* returns n if successful, k if D (k,k) is zero */
(
    ADMMint n,		/* A and L are n-by-n, where n >= 0 */
    ADMMint *Ap,	/* input of size n+1, not modified */
    ADMMint *Ai,	/* input of size nz=Ap[n], not modified */
    ADMMfloat *Ax,	/* input of size nz=Ap[n], not modified */
    ADMMint *Lp,	/* input of size n+1, not modified */
    ADMMint *Parent,	/* input of size n, not modified */
    ADMMint *Lnz,	/* output of size n, not defn. on input */
    ADMMint *Li,	/* output of size lnz=Lp[n], not defined on input */
    ADMMfloat *Lx,	/* output of size lnz=Lp[n], not defined on input */
    ADMMfloat *D,	/* output of size 2n, not defined on input */
    ADMMfloat *Y,	/* workspace of size 2n, not defn. on input or output */
    ADMMint *Pattern,/* workspace of size n, not defn. on input or output */
    ADMMint *Flag,	/* workspace of size n, not defn. on input or output */
    ADMMint *b2,   /* Vector of size n, indicates if a 1x1 or 2x2 pivot is needed*/
    ADMMint *P,	/* optional input of size n */
    ADMMint *Pinv	/* optional input of size n */
);

/* ========================================================================== */
/* === bk_symbolic: Determine column sizes of decomposed L=================== */
/* ========================================================================== */

void bk_symbolic
(
    ADMMint n,		/* A and L are n-by-n, where n >= 0 */
    ADMMint *Ap,	/* input of size n+1, not modified */
    ADMMint *Ai,	/* input of size nz=Ap[n], not modified */
    ADMMint *Lp,	/* output of size n+1, not defined on input */
    ADMMint *Parent,	/* output of size n, not defined on input */
    ADMMint *Lnz,	/* output of size n, not defined on input */
    ADMMint *Flag,	/* workspace of size n, not defn. on input or output */
    ADMMint *b2,     /*vector of size n+1, indicates when 2x2 pivots are used*/
    ADMMint *P,	/* optional input of size n */
    ADMMint *Pinv	/* optional output of size n (used if P is not NULL) */
);

/* ========================================================================== */
/* === bk_ltsolve: solves L^Tx = b ========================================== */
/* ========================================================================== */

void bk_ltsolve
(
    ADMMint n,		/* L is n-by-n, where n >= 0 */
    ADMMfloat *X,	/* size n.  right-hand-side on input, soln. on output */
    ADMMint *Lp,	/* input of size n+1, not modified */
    ADMMint *Li,	/* input of size lnz=Lp[n], not modified */
    ADMMfloat *Lx,	/* input of size lnz=Lp[n], not modified */
    ADMMint *Lnz /*nnz per column*/
);

/* ========================================================================== */
/* === bk_lsolve: solves Lx = b ============================================= */
/* ========================================================================== */

void bk_lsolve
(
    ADMMint n,		/* L is n-by-n, where n >= 0 */
    ADMMfloat *X,	/* size n.  right-hand-side on input, soln. on output */
    ADMMint   *Lp,	/* input of size n+1, not modified */
    ADMMint   *Li,	/* input of size lnz=Lp[n], not modified */
    ADMMfloat *Lx,	/* input of size lnz=Lp[n], not modified */
    ADMMint   *Lnz /*nnz per column*/
);