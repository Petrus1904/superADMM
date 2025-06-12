#ifndef _CS_H
#define _CS_H

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 1		    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"	    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

/*-------Custom overwrite to include compiling in single precision float format----------*/
#ifndef USE_SINGLE_PRECISION
    #define CSfloat double
#else
    #define CSfloat float
#endif

#ifdef MATLAB_COMP
    #define CSint long long int
#else
    #define CSint int
#endif

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    CSint nzmax ;	    /* maximum number of entries */
    CSint m ;	    /* number of rows */
    CSint n ;	    /* number of columns */
    CSint *p ;	    /* column pointers (size n+1) or col indices (size nzmax) */
    CSint *i ;	    /* row indices, size nzmax */
    CSfloat *x ;	    /* numerical values, size nzmax */
    CSint nz ;	    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add (const cs *A, const cs *B, CSfloat alpha, CSfloat beta) ;
CSint cs_cholsol (const cs *A, CSfloat *b, CSint order) ;
CSint cs_dupl (cs *A) ;
CSint cs_entry (cs *T, CSint i, CSint j, CSfloat x) ;
CSint cs_lusol (const cs *A, CSfloat *b, CSint order, CSfloat tol) ;
CSint cs_gaxpy (const cs *A, const CSfloat *x, CSfloat *y) ;
cs *cs_multiply (const cs *A, const cs *B) ;
CSint cs_qrsol (const cs *A, CSfloat *b, CSint order) ;
cs *cs_transpose (const cs *A, CSint values) ;
cs *cs_triplet (const cs *T) ;
CSfloat cs_norm (const cs *A) ;
CSint cs_print (const cs *A, CSint brief) ;
cs *cs_load (FILE *f) ;
/* utilities */
void *cs_calloc (CSint n, size_t size) ;
void *cs_free (void *p) ;
void *cs_realloc (void *p, CSint n, size_t size, CSint *ok) ;
cs *cs_spalloc (CSint m, CSint n, CSint nzmax, CSint values, CSint triplet) ;
cs *cs_spfree (cs *A) ;
CSint cs_sprealloc (cs *A, CSint nzmax) ;
void *cs_malloc (CSint n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    CSint *Pinv ;	    /* inverse row perm. for QR, fill red. perm for Chol */
    CSint *Q ;	    /* fill-reducing column permutation for LU and QR */
    CSint *parent ;   /* elimination tree for Cholesky and QR */
    CSint *cp ;	    /* column pointers for Cholesky, row counts for QR */
    CSint *b2 ;
    CSint m2 ;	    /* # of rows for QR, after adding fictitious rows */
    CSint lnz ;	    /* # entries in L for LU or Cholesky; in V for QR */
    CSint unz ;	    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;	    /* L for LU and Cholesky, V for QR */
    cs *U ;	    /* U for LU, R for QR, not used for Cholesky */
    CSint *Pinv ;	    /* partial pivoting for LU */
    CSfloat *B ;	    /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    CSint *P ;	    /* size m, row permutation */
    CSint *Q ;	    /* size n, column permutation */
    CSint *R ;	    /* size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q) */
    CSint *S ;	    /* size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q) */
    CSint nb ;	    /* # of blocks in fine dmperm decomposition */
    CSint rr [5] ;    /* coarse row decomposition */
    CSint cc [5] ;    /* coarse column decomposition */
} csd ;

CSint *cs_amd (const cs *A, CSint order) ;
csn *cs_chol (const cs *A, const css *S) ;
csd *cs_dmperm (const cs *A) ;
CSint cs_droptol (cs *A, CSfloat tol) ;
CSint cs_dropzeros (cs *A) ;
CSint cs_happly (const cs *V, CSint i, CSfloat beta, CSfloat *x) ;
CSint cs_ipvec (CSint n, const CSint *P, const CSfloat *b, CSfloat *x) ;
CSint cs_lsolve (const cs *L, CSfloat *x) ;
CSint cs_ltsolve (const cs *L, CSfloat *x) ;
csn *cs_lu (const cs *A, const css *S, CSfloat tol) ;
cs *cs_permute (const cs *A, const CSint *P, const CSint *Q, CSint values) ;
CSint *cs_pinv (const CSint *P, CSint n) ;
CSint cs_pvec (CSint n, const CSint *P, const CSfloat *b, CSfloat *x) ;
csn *cs_qr (const cs *A, const css *S) ;
css *cs_schol (const cs *A, CSint order) ;
css *cs_sqr (const cs *A, CSint order, CSint qr) ;
cs *cs_symperm (const cs *A, const CSint *Pinv, CSint values) ;
CSint cs_usolve (const cs *U, CSfloat *x) ;
CSint cs_utsolve (const cs *U, CSfloat *x) ;
CSint cs_updown (cs *L, CSint sigma, const cs *C, const CSint *parent) ;
/* utilities */
css *cs_sfree (css *S) ;
csn *cs_nfree (csn *N) ;
csd *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
CSint *cs_counts (const cs *A, const CSint *parent, const CSint *post, CSint ata) ;
CSint cs_cumsum (CSint *p, CSint *c, CSint n) ;
CSint cs_dfs (CSint j, cs *L, CSint top, CSint *xi, CSint *pstack, const CSint *Pinv) ;
CSint *cs_etree (const cs *A, CSint ata) ;
CSint cs_fkeep (cs *A, CSint (*fkeep) (CSint, CSint, CSfloat, void *), void *other) ;
CSfloat cs_house (CSfloat *x, CSfloat *beta, CSint n) ;
CSint *cs_maxtrans (const cs *A) ;
CSint *cs_post (CSint n, const CSint *parent) ;
CSint cs_reach (cs *L, const cs *B, CSint k, CSint *xi, const CSint *Pinv) ;
csd *cs_scc (cs *A) ;
CSint cs_scatter (const cs *A, CSint j, CSfloat beta, CSint *w, CSfloat *x, CSint mark,
    cs *C, CSint nz) ;
CSint cs_splsolve (cs *L, const cs *B, CSint k, CSint *xi, CSfloat *x,
    const CSint *Pinv) ;
CSint cs_tdfs (CSint j, CSint k, CSint *head, const CSint *next, CSint *post,
    CSint *stack) ;
/* utilities */
csd *cs_dalloc (CSint m, CSint n) ;
cs *cs_done (cs *C, void *w, void *x, CSint ok) ;
CSint *cs_idone (CSint *p, cs *C, void *w, CSint ok) ;
csn *cs_ndone (csn *N, cs *C, void *w, void *x, CSint ok) ;
csd *cs_ddone (csd *D, cs *C, void *w, CSint ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (int) size)
#endif
