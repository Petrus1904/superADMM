/*
    Source file for the Bunch Kaufman sparse LDL decomposition and solver
    Created by Peter Verheijen, 2025.

    The execution of the decomposition and symbolic factoring is build/heavily inspired from the LDL-
    implementation made by Timothy A. Davis, which is also the LDL method currently
    employed in superADMM. However, as of this point, the code is substantially different.

    In short: the Bunch-Kaufman LDL decomposition considers that instead of just a diagonal D,
    some diagonal parts are 2x2 blocks instead. As a result, the method is much more robust against
    numerically poorly scaled matrices. Additionally, 2x2 matrices can be easily inverted using Cramers rule
    and have a low computational overhead. In practice, a version BK-LDL decomp is employed in the MA57 algorithm, 
    which is used by intel-MKL and MATLAB (MATLAB only for solving indefinite non-singular sparse matrices).

    The BK-decomposition made here is completely handwritten by me (Peter). Performance wise, it does show that MA57
    remains more accurate, and definitely faster for large scale matrices. For smaller matrices (say less than 1000x1000)
    this BK-LDL decomp could execute faster. Additionally, BK-LDL produces more accurate solutions that naive LDL.

    Despite the fact that this BK-sparse LDL decomposition does work*, it is currently not
    employed in the superADMM solver, as the resulting decomposition is substantially
    more dense, which slows down the over al solver time. In fact, I made a couple observations:
     - BK needs to recompute permutation & symbolic every iteration, as the permutation depends now on the value scale too.
     - Due to increased number of fill-ins (as the permutation is now less efficient), 
       BK needs more memory and can run to overflows for large problems.
     - BK is in general slower than LDL because of 2x2 pivots
     - However, using BK instead of LDL reduces the iteration count by 40%!
     - BK can solve problems that are numerically terribly scaled (very nice!)
     - In general --> BK scales terrible with problem size, and is as of now never faster than standard LDL.

    Some comparison results, note that the sparse matrices are generated in a similar structure as commonly seen in superADMM
    | size n | nnz    | tMA57    | tLDL      | tBKLDL     | epsMA57 | epsLDL  | epsBKLDL |
    +--------+--------+----------+-----------+------------+---------+---------+----------+
    | 30     | 84     | 0.053ms  | 0.028ms   | 0.035ms    | 2.2e-15 | 3.4e-08 | 8.9e-16  |
    | 150    | 476    | 0.309ms  | 0.200ms   | 0.168ms    | 4.3e-15 | 5.7e-08 | 5.3e-15  |
    | 750    | 2338   | 0.924ms  | 1.046ms   | 0.473ms    | 2.9e-15 | 7.0e-08 | 4.2e-14  |
    | 1200   | 4672   | 2.809ms  | 1.204ms   | 3.581ms    | 5.8e-15 | 1.4e-07 | 5.7e-14  |
    | 1500   | 6146   | 1.871ms  | 2.136ms   | 19.566ms   | 3.8e-15 | 1.1e-07 | 1.3e-13  |
    | 1800   | 8880   | 4.046ms  | 2.221ms   | 35.592ms   | 6.2e-15 | 1.1e-07 | 1.2e-13  |
    | 3150   | 25330  | 9.632ms  | 8.201ms   | 289.726ms  | 1.5e-14 | 2.2e-07 | 7.4e-12  |
    | 5250   | 42630  | 19.309ms | 17.343ms  | 2387.984ms | 1.0e-14 | 2.8e-07 | 2.3e-11  |
    | 6750   | 33928  | 10.283ms | 8.854ms   | 9226.626ms | 6.2e-15 | 1.5e-07 | 4.3e-11  |
    +--------+--------+----------+-----------+------------+---------+---------+----------+
    Which shows clearly that this BKLDL algorithm scales badly... I recon there are some improvements possible on the pivoting rule.
    

    * I say "work", but in no way should this code be taken for granted, verify inputs, check data, and verify results.
      Even now, I cannot ensure the code is safe to use, there can still be memory issues hiding.
*/

// #include <math.h>
#include <string.h>
#include "BK_ldl.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define abs(a) (a < 0 ? -a : a)

static ADMMint countMatchingCols(ADMMint *A, ADMMint n, ADMMint *B, ADMMint m){
    ADMMint cnt = 0;
    ADMMint i = 0;
    ADMMint j = 0;
    //quick tests
    if(A[0] > B[m-1]){
        //nothing in range
        return 0;
    }
    if(B[0] > A[n-1]){
        return 0;
    }
    while(i<n && j<m){
        if(A[i] == B[j]){
            cnt++;
            i++;
            j++;
        } else if(A[i]<B[j]){
            i++;
        } else {
            j++;
        }
    }
    return cnt;
}

/*==============================================================================
bk_permute computes the Bunch-Kaufman row-column permutation, which changes
The existing permutation given on the input (as a result of symbolic decomp for 
instance), and returns a vector indicating when a 1x1 pivot must be done, and when a
2x2 pivot is expected. This contains the elemental B-K pivoting rules, with an
inclusion to prevent dense fill-ins (mind that we also need to keep this minimal)

There are two sorting orders to perform the permutation, forwards and backwards 
(forwards is commented currently). Backwards seems to work better with the AMD
procedure, which favors dense rows at the back. This can indeed reduce the chance
of swapping a very dense row further down with a sparse row earlier up, which makes
the LDL decomposition very dense.
================================================================================*/
ADMMint bk_permute
(
    const ADMMint n,            /*col/rows of A*/
    const ADMMfloat *Ax,        /*const nonzero values of A*/
    const ADMMint *Ai,          /*const row indices of A*/
    const ADMMint *Ap,          /*const column pointers of A*/
    ADMMint *P,                 /*in/out permutation vector*/
    ADMMint *Pinv,              /*in/out inverse of P*/
    ADMMint *b2                 /*output vector of size n+1 defining when 2x2 pivots are needed*/
)
{
    ADMMint kk, r, t, cnt, p2;
    ADMMfloat infmax, sigma, d;
    const ADMMfloat alpha =  0.525; // 0.6404;
    cnt = 0;
    // Forwards implementation.
    b2[0] = 0;
    b2[1] = 1; //skip first index -- following MA57 standard
    for(ADMMint k = 0; k < n; k++){
        if(k == 0 || b2[k-1] != b2[k]){
            kk = P[k];
            //inf norm on column
            infmax = 0;
            d = 0;
            p2 = Ap[kk+1];
            for(ADMMint p = Ap[kk]; p < p2; p++){
                ADMMint ac = Pinv[Ai[p]];
                if(ac == k){
                    d = abs(Ax[p]);
                }
                if(ac > k){
                    if(infmax < abs(Ax[p])){
                        infmax = abs(Ax[p]);
                        r = ac;
                    }
                }
            }
            if(infmax*alpha > d){
                sigma = 0;
                p2 = Ap[P[r]+1];
                for(ADMMint p = Ap[P[r]]; p < p2; p++){
                    if(Pinv[Ai[p]] > k+1){
                        sigma = max(sigma,abs(Ax[p]));
                    }   
                }
                if(alpha*infmax*infmax < sigma*d){
                    //1x1 pivot
                    b2[k+1] = b2[k]+1;
                } else {
                    //swap
                    t       = P[k+1];
                    P[k+1]  = P[r];
                    P[r]    = t;
                    Pinv[P[k+1]] = k+1;
                    Pinv[P[r]] = r;
                    //2x2 pivot
                    b2[k+1] = b2[k];
                    cnt++;
                }
            } else {
                //1x1 pivot
                b2[k+1] = b2[k]+1;
            }
        } else {
            //size for 2x2 pivot
            b2[k+1] = b2[k]+3;
        }
    }
    // b2[n] = 1;
    // // backwards implementation
    // for(ADMMint k = n-1; k >= 0; k--){
    //     if(b2[k+1] != 3){
    //         kk = P[k];
    //         //inf norm on column
    //         infmax = 0;
    //         d = 0;
    //         p2 = Ap[kk+1];
    //         for(ADMMint p = Ap[kk]; p < p2; p++){
    //             ADMMint ac = Pinv[Ai[p]];
    //             if(ac == k){
    //                 d = abs(Ax[p]);
    //             } else if(ac < k){
    //                 if(infmax < abs(Ax[p])){
    //                     infmax = abs(Ax[p]);
    //                     r = ac;
    //                 }
    //             }
    //         }
    //         if(infmax*alpha > d){
    //             sigma = 0;
    //             p2 = Ap[P[r]+1];
    //             for(ADMMint p = Ap[P[r]]; p < p2; p++){
    //                 if(Pinv[Ai[p]] < k-1){
    //                     sigma = max(sigma,abs(Ax[p]));
    //                 }   
    //             }
    //             if(alpha*infmax*infmax < sigma*d){
    //                 //1x1 pivot
    //                 b2[k] = 1;
    //             } else {
    //                 //swap
    //                 t       = P[k-1];
    //                 P[k-1]  = P[r];
    //                 P[r]    = t;
    //                 Pinv[P[k-1]] = k-1;
    //                 Pinv[P[r]] = r;
    //                 //2x2 pivot
    //                 b2[k] = 3;
    //                 cnt++;
    //             }
    //         } else {
    //             //1x1 pivot
    //             b2[k] = 1;
    //         }
    //     } else {
    //         //size for 2x2 pivot
    //         b2[k] = 0;
    //     }
    // }
    // //cumsum
    // ADMMint ii = 0;
    // ADMMint tt;
    // for(ADMMint k = 0; k <= n; k++){
    //     tt = b2[k];
    //     b2[k] = ii;
    //     ii += tt;
    // }
    return cnt;
}

/* ========================================================================== */
/* === bk_dsolve:  solve Dx=b ============================================== */
/* ========================================================================== */

void bk_dsolve
(
    ADMMint n,		    /* length of x */
    ADMMfloat *x,	    /* b on the input, x on the output */
    ADMMfloat *D,	    /* D vector of size b2[n], block diagonal */
    ADMMint *b2         /*when to do 1x1 of when 2x2 inverses */
)
{
    ADMMint k, cc;
    ADMMfloat c1, c2, inv;
    cc = 0;
    for (ADMMint j = 0; j < n; j++)
    {
        k = b2[j];
        if(cc == 0){ //this is much cheaper to compute
            if(k == b2[j+1]){
                c1 = x[j];
                c2 = x[j+1];
                x[j] = (c1*D[k]+D[k+1]*c2);
                x[j+1] = (c2*D[k+2]+D[k+1]*c1);
                cc = 1;
            } else {
                x[j] /= D[k]; //inverse only for diagonals
            }
        } else {
            cc = 0;
        }
    }
}

/* ========================================================================== */
/* === ldl_numeric ========================================================== */
/* ========================================================================== */

/* Given a sparse matrix A (the arguments n, Ap, Ai, and Ax) and its symbolic
 * analysis (Lp and Parent, and optionally P and Pinv), compute the numeric LDL'
 * factorization of A or PAP'.  The outputs of this routine are arguments Li,
 * Lx, and D.  It also requires three size-n workspaces (Y, Pattern, and Flag).
 * Original from TA Davis, Modified by me (Peter) to do also 2x2 pivots
 */

 ADMMint bk_numeric	/* returns n if successful, k if D (k,k) is zero */
 (
    ADMMint n,		    /* A and L are n-by-n, where n >= 0 */
    ADMMint *Ap,	    /* input of size n+1, not modified */
    ADMMint *Ai,	    /* input of size nz=Ap[n], not modified */
    ADMMfloat *Ax,	    /* input of size nz=Ap[n], not modified */
    ADMMint *Lp,	    /* input of size n+1, not modified */
    ADMMint *Parent,	/* input of size n, not modified */
    ADMMint *Lnz,	    /* output of size n, not defn. on input */
    ADMMint *Li,	    /* output of size lnz=Lp[n], not defined on input */
    ADMMfloat *Lx,	    /* output of size lnz=Lp[n], not defined on input */
    ADMMfloat *D,	    /* output of size dnz = b2[n], not defined on input */
    ADMMfloat *Y,	    /* workspace of size 2n, not defn. on input or output */
    ADMMint *Pattern,   /* workspace of size n, not defn. on input or output */
    ADMMint *Flag,	    /* workspace of size n, not defn. on input or output */
    ADMMint *b2,        /* Vector of size n+1, not modified, indicates if a 1x1 or 2x2 pivot is needed*/
    ADMMint *P,	        /* permutation vector of size n */
    ADMMint *Pinv	    /* inv. permutation vector of size n */
 )
 {
    ADMMfloat yi, l_ki, l_ko, inv, d1, d2, d3;
    ADMMfloat *Y2 = Y + n;
    ADMMint i, k, p, kk, p2, len, top, j, jj, cc;
    cc = 0;
    if(!Y) return -2;
    if(!Lnz) return -3;
    if(!P) return -4;
    if(!b2) return -5;
    if(!Lp) return -6;
    if(!Flag) return -7;
    for (k = 0; k < n; k++)
    {
        /* compute nonzero Pattern of kth row of L, in topological order */
        Y[k] = 0.0;		    /* Y(0:k) is now all zero */
        Y2[k] = 0.0;       /* Y(n+1:k) is also zero */
        top = n;		    /* stack for pattern is empty */
        Flag[k] = k;		    /* mark node k as visited */
        Lnz[k] = Lp[k];		    /* count of nonzeros in column k of L */
        kk = P[k];  /* kth original, or permuted, column */
        p2 = Ap[kk+1];
        jj = b2[k];
        if(cc == 1){ //if previous one was the top corner of the block, the next is the bottom corner!
            cc = 2;
            Flag[k-1] = k; //mark as visited
        } else if(jj == b2[k+1]){
            cc = 1;
        } else {
            cc = 0;
        }
        for (p = Ap[kk]; p < p2; p++)
        {
            i = Pinv[Ai[p]];	/* get A(i,k) */
            if (i <= k)
            {
                Y[i] = Ax[p];  /* scatter A(i,k) into Y (do not sum duplicates) */
                for (len = 0; Flag[i] != k; i = Parent[i])
                {
                    Pattern[len++] = i;   /* L(k,i) is nonzero */
                    Flag[i] = k;	    /* mark i as visited */
                }
                
                while (len > 0) Pattern[--top] = Pattern[--len];
            }
        }
        /* compute numerical values kth row of L (a sparse triangular solve) */
        if(cc == 2){
            D[jj+1] = Y[k-1]; //this works because of symmetry :)
            D[jj+2] = Y[k];
            Y[k-1] = 0.0; //also make this zero as the block ends at k-1
        } else {
            D[jj] = Y[k];
        }
        Y[k] = 0.0;
        
        for (; top < n; top++)
        {
            //i < k -- always
            i = Pattern[top];	    /* Pattern[top:n-1] is pattern of L(k,:) */
            yi = Y[i];	            /* get and clear Y(i) */
            Y[i] = 0.0;
            p2 = Lnz[i]; //Lp[i] + 
            if(p2 > 0 && Li[p2-1] == k){
                //we already set this entry last iteration, just make sure we stop at k-1
                p2--;
            }
            for (p = Lp[i]; p < p2; p++)
            {
                Y[Li[p]] -= Lx[p] * yi;
            }
            j = b2[i];
            if(b2[i+1] == j){
                d1 = D[j];
                d2 = D[j+1];
                if(Li[p] == k){
                    
                    l_ki = 2*Lx[p] + d1*yi;
                    p2 = Lnz[i+1]-1; //check last entry of previous column

                    Lx[p]  += d1*yi;
                    //update Lx[p2]
                    Lx[p2] += d2*yi;
                } else {
                    //and append for the next column
                    p2 = Lnz[i+1];

                    //update L[p] and L[p+1]
                    l_ki    = yi*d1;
                    Lx[p]   = l_ki;
                    Lx[p2]  = d2*yi;
                    Li[p]   = k;
                    Li[p2]  = k;
                    Lnz[i]++;
                    Lnz[i+1]++;
                }
                
                if(cc == 2) l_ko = d1*Y2[i] + d2*Y2[i+1];                
            } else if(i > 0 && j == b2[i-1]){
                d2 = D[j+1];
                d3 = D[j+2];
                if(Li[p] == k){
                    
                    l_ki = 2*Lx[p] + d3*yi; //2*Lx[p] fixes the cross terms for D.
                    Lx[p] += d3*yi;
                    //update Lx[p2]
                    p2 = Lnz[i-1]-1; //check last entry of previous column
                    Lx[p2] += d2*yi;
                } else {
                    //only update if the previous y was zero
                    //otherwise everything is already done there
                    p2 = Lnz[i-1];   //check last entry of previous column
                    
                    l_ki    = yi*d3;
                    Lx[p]   = l_ki;
                    //also update the previous one
                    Lx[p2]  = d2*yi;
                    Li[p2]  = k;
                    Li[p]   = k;
                    Lnz[i]++;
                    Lnz[i-1]++;
                    
                }
                if(cc == 2) l_ko = d3*Y2[i] + d2*Y2[i-1];
            } else {
                l_ki = yi / D[j];	    /* the nonzero entry L(k,i) */
                Li[p] = k;	            /* store L(k,i) in column form of L */
                Lx[p] = l_ki;
                Lnz[i]++;
                if(cc == 2) l_ko = Y2[i]/D[j];
            }
            //update D
            if(cc == 1){
                D[jj] -= l_ki * yi;
                //store y
                Y2[i] = yi;
            } else if(cc == 2){
                D[jj+2] -= l_ki * yi;
                D[jj+1] -= l_ko * yi; //Y[n+i];
            } else {
                D[jj] -= l_ki * yi;
            }
        }
        if(cc==2){
            //new D block computed, invert it!
            d1 = D[jj];
            d2 = D[jj+1];
            d3 = D[jj+2];
            inv = d1*d3-d2*d2;
            if(inv == 0.0) return(k); /*failure, det(D(k-1:k,k-1:k)) is zero*/
            D[jj] = d3/inv;
            D[jj+1] = -d2/inv;
            D[jj+2] = d1/inv;

            memset(Y2, 0.0, (k+1) * sizeof(ADMMfloat));
        }
        if(k>0) Y[k-1] = 0.0; //security check
        if (cc == 0.0 && D[k] == 0.0) return (k);/* failure, D(k,k) is zero */
    }
    return (n);	/* success, diagonal of D is all nonzero */
 }

 /* ========================================================================== */
/* === ldl_symbolic ========================================================= */
/* ========================================================================== */

/* The input to this routine is a sparse matrix A, stored in column form, and
 * an optional permutation P.  The output is the elimination tree
 * and the number of nonzeros in each column of L.  Parent[i] = k if k is the
 * parent of i in the tree.  The Parent array is required by ldl_numeric.
 * Lnz[k] gives the number of nonzeros in the kth column of L, excluding the
 * diagonal.
 *
 * One workspace vector (Flag) of size n is required.
 *
 * If P is NULL, then it is ignored.  The factorization will be LDL' = A.
 * Pinv is not computed.  In this case, neither P nor Pinv are required by
 * ldl_numeric.
 *
 * If P is not NULL, then it is assumed to be a valid permutation.  If
 * row and column j of A is the kth pivot, the P[k] = j.  The factorization
 * will be LDL' = PAP', or A (p,p) in MATLAB notation.  The inverse permutation
 * Pinv is computed, where Pinv[j] = k if P[k] = j.  In this case, both P
 * and Pinv are required as inputs to ldl_numeric.
 *
 * The floating-point operation count of the subsequent call to ldl_numeric
 * is not returned, but could be computed after ldl_symbolic is done.  It is
 * the sum of (Lnz[k]) * (Lnz[k] + 2) for k = 0 to n-1.
 */

void bk_symbolic
(
    ADMMint n,		    /* A and L are n-by-n, where n >= 0 */
    ADMMint *Ap,	    /* input of size n+1, not modified */
    ADMMint *Ai,	    /* input of size nz=Ap[n], not modified */
    ADMMint *Lp,	    /* output of size n+1, not defined on input */
    ADMMint *Parent,	/* output of size n, not defined on input */
    ADMMint *Lnz,	    /* output of size n, not defined on input */
    ADMMint *Flag,	    /* workspace of size n, not defn. on input or output */
    ADMMint *b2,        /* vector of size n+1, indicates when 2x2 pivots are used*/
    ADMMint *P,	        /* permutation vector of size n */
    ADMMint *Pinv	    /* inv. permutation vector of size n */
)
{
    ADMMint i, k, p, kk, p2;
    for (k = 0; k < n; k++)
    {
        /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
        Parent[k] = -1;	    /* parent of k is not yet known */
        Flag[k] = k;		    /* mark node k as visited */
        Lnz[k] = 0;		    /* count of nonzeros in column k of L */
        kk = P[k];  /* kth original, or permuted, column */
        p2 = Ap[kk+1];
        for (p = Ap[kk]; p < p2; p++)
        {
            /* A (i,k) is nonzero (original or permuted A) */
            i = Pinv[Ai[p]];
            if (i < k)
            {
                /* follow path from i to root of etree, stop at flagged node */
                for (; Flag[i] != k; i = Parent[i])
                {
                    /* find parent of i if not yet determined */
                    if (Parent[i] == -1) Parent[i] = k;
                    Lnz[i]++;				/* L (k,i) is nonzero */
                    Flag[i] = k;			/* mark i as visited */
                }
            }
        }
    }
    /* construct Lp index array from Lnz column counts */
    Lp[0] = 0;
    ADMMint cn = 0;
    ADMMint cc = 0;
    for (k = 0; k < n; k++)
    {
        cn = 0;
        if(b2[k] == b2[k+1]){
            cn = Lnz[k+1];
            cc = 1;
        }
        else if(cc == 1){
            cn = Lnz[k-1];
            cc = 0;
        }
	    Lp[k+1] = Lp[k] + Lnz[k]+cn;
    }
}

/* ========================================================================== */
/* === ldl_lsolve:  solve Lx=b ============================================== */
/* ========================================================================== */

void bk_lsolve
(
    ADMMint n,		/* L is n-by-n, where n >= 0 */
    ADMMfloat *X,	/* size n.  right-hand-side on input, soln. on output */
    ADMMint *Lp,	/* input of size n+1, not modified */
    ADMMint *Li,	/* input of size lnz=Lp[n], not modified */
    ADMMfloat *Lx,	/* input of size lnz=Lp[n], not modified */
    ADMMint *Lnz    /* vector of size n, nnz per column*/
)
{
    ADMMint j, p, p2;
    ADMMfloat xj;
    for (j = 0; j < n; j++)
    {
        p = Lp[j];
        p2 = Lnz[j]; //p+
        xj = X[j];
        for (; p < p2; p++)
        {
            X[Li[p]] -= Lx[p] * xj;
        }
    }
}

/* ========================================================================== */
/* === ldl_ltsolve: solve L'x=b  ============================================ */
/* ========================================================================== */

void bk_ltsolve
(
    ADMMint n,		/* L is n-by-n, where n >= 0 */
    ADMMfloat *X,	/* size n.  right-hand-side on input, soln. on output */
    ADMMint *Lp,	/* input of size n+1, not modified */
    ADMMint *Li,	/* input of size lnz=Lp[n], not modified */
    ADMMfloat *Lx,	/* input of size lnz=Lp[n], not modified */
    ADMMint *Lnz    /* vector of size n, nnz per column*/
)
{
    ADMMint j, p, p2;
    ADMMfloat xj;
    for (j = n-1; j >= 0; j--)
    {
        p = Lp[j];
        p2 = Lnz[j]; //p+
        xj = X[j];
        for (; p < p2; p++)
        {
            xj -= Lx[p] * X[Li[p]];
        }
        X[j] = xj;
    }
}