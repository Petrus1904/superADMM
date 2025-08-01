/*
    Header file for the superADMM solver
    Created by Peter Verheijen, 2024.
*/

#ifndef SUPERADMM_H
#define SUPERADMM_H

#ifdef BUILD_DLL
#define BUILDTAG  __declspec(dllexport)
#else
#define BUILDTAG 
#endif

#ifdef MATLAB_COMP
    #define ADMMint long long int
#else
    #define ADMMint int
#endif

#ifndef USE_SINGLE_PRECISION
    #define ADMMfloat double
    #define cblas_ger cblas_dger
    #define cblas_syr cblas_dsyr
    #define cblas_copy cblas_dcopy
    #define cblas_scal cblas_dscal
    #define cblas_syrk cblas_dsyrk
    #define cblas_gemv cblas_dgemv
    #define cblas_symv cblas_dsymv
    #define cblas_trsv cblas_dtrsv
    #define cblas_amax cblas_damax
    #define cblas_axpy cblas_daxpy
    #define cblas_dot cblas_ddot

    // extern void dgesv_(const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, ADMMint *ipiv, double *b, const ADMMint *ldb, ADMMint *info);
    extern void dposv(const char *uplo, const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, double *b, const ADMMint *ldb, ADMMint *info);
    // extern void dlaset_(const char *uplo, const ADMMint *m, const ADMMint *n, const double *alpha, const double *beta, double *a, const ADMMint *lda);
    extern void dlacpy(const char* uplo, const ADMMint *m, const ADMMint *n, const double* a, const ADMMint *lda, double* b, const ADMMint *ldb);
    // extern void dsysv_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, ADMMint *ipiv, double *b, const ADMMint *ldb, double *work, const ADMMint *lwork, ADMMint *info);

    // #define lapack_gesv dgesv_
    #define lapack_posv dposv
    // #define lapack_laset dlaset_
    #define lapack_lacpy dlacpy
    // #define lapack_sysv dsysv_
#else
    #define ADMMfloat float
    #define cblas_ger cblas_sger
    #define cblas_copy cblas_scopy
    #define cblas_scal cblas_sscal
    #define cblas_syrk cblas_ssyrk
    #define cblas_gemv cblas_sgemv
    #define cblas_symv cblas_ssymv
    #define cblas_amax cblas_samax
    #define cblas_axpy cblas_saxpy
    #define cblas_dot cblas_sdot

    extern void sgesv_(const ADMMint *n, const ADMMint *nrhs, float *a, const ADMMint *lda, ADMMint *ipiv, float *b, const ADMMint *ldb, ADMMint *info);
    extern void sposv_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, float *a, const ADMMint *lda, float *b, const ADMMint *ldb, ADMMint *info);
    extern void slaset_(const char *uplo, const ADMMint *m, const ADMMint *n, const float *alpha, const float *beta, float *a, const ADMMint *lda);
    extern void slacpy_(const char* uplo, const int* m, const int* n, const float* a, const int* lda, float* b, const int* ldb);
    extern void ssysv_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, float *a, const ADMMint *lda, ADMMint *ipiv, float *b, const ADMMint *ldb, float *work, const ADMMint *lwork, ADMMint *info);

    #define lapack_gesv sgesv_
    #define lapack_posv sposv_
    #define lapack_laset slaset_
    #define lapack_lacpy slacpy_
    #define lapack_sysv ssysv_
#endif

typedef struct ADMMopts{
    ADMMint verbose;
    ADMMint maxIter;
    ADMMfloat sigma;    //slack on x -- ensures posDEF in P+ARA
    ADMMfloat rho_0;    //starting coeff
    ADMMfloat tau;      //exponential Bound decrease
    ADMMfloat alpha;    //exponential increase
    ADMMfloat RBound;   //limit on R
    ADMMfloat eps_abs;  //absolute tolerance
    ADMMfloat eps_inf;  //relative tolerance
    ADMMint repInterval; //reporting interval
    ADMMfloat timeLimit;//execution time limit, use 0 for unlimited
    ADMMfloat lowRankPer; //percentage required to execute low rank updates instead of full LDL, default 0.05 (5%)
} ADMMopts;

//default options
#define DEFVERBOSE 1
#define DEFMAXITER 500
#define DEFSIGMA   1e-6
#define DEFRHO0    1.0
#define DEFTAU     0.5
#define DEFALPHA   500.0
#define DEFRBOUND  1.0e8
#define DEFEPSABS  1e-8
#define DEFEPSINF  1e-8
#define DEFREPIVAL 10
#define DEFTIMELIM 0.0
#define DEFLRPER   0.05

typedef struct ADMMinfo{
    ADMMint nIter;          //number of iterations it took
    ADMMfloat rPrim;    //primal convergence score
    ADMMfloat rDual;    //dual convergence score
    ADMMfloat runtime;  //run time in seconds
    ADMMfloat objVal;   //The objective value
    char* status;       //the solver status
} ADMMinfo;

typedef struct csLDL{
    ADMMint n;          /* length of L, D*/
    ADMMint nnz;        /* number of nonzeros in Lx */
    ADMMfloat *D;   /* n   D, size b2[n] for BK */
    ADMMint *Lp;        /* n+1 column pointers in L*/
    ADMMint *Li;        /* nnz row indexes in L */
    ADMMfloat *Lx;  /* nnz L data */
    ADMMint *P;         /* n   Permutation from AMD */
    ADMMint *Pinv;      /* n   inverse permutation, Pinv[P[i]] = i */
    ADMMint *Parent;    /* n   etree from symbolic */
    ADMMint *Lnz;       /* n   nnz per column in L */
    ADMMint *b2;        /* n+1 2x2 blocks indexes only for BK, unused in native */
    ADMMint *Flag;      /* n   workspace data for LDL */
    ADMMint *Pattern;   /* n   workspace data for LDL */
    ADMMfloat *Y;   /* n   workspace data for LDL, 2n in BK */
} csldl;

BUILDTAG ADMMint superADMMsolverDense(const ADMMfloat* P, const ADMMfloat* q, const ADMMfloat* A, const ADMMfloat* l, const ADMMfloat* u,
                                        ADMMfloat* x, ADMMfloat* y, ADMMint nPrim, ADMMint nDual, ADMMopts opts, ADMMinfo* info);

BUILDTAG ADMMint superADMMsolverSparse(ADMMfloat* Pdata, ADMMint *Prowptr, ADMMint *Pcolidx, const ADMMfloat* q, 
                                        ADMMfloat* Adata, ADMMint *Arowptr, ADMMint *Acolidx, const ADMMfloat* l, const ADMMfloat* u,
                                        ADMMfloat* x, ADMMfloat* y, ADMMint nPrim, ADMMint nDual, ADMMopts opts, ADMMinfo* info);

#endif