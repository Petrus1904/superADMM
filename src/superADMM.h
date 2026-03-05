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

#define VERSION "0.7.1"

#ifdef MATLAB_COMP
    /*
    To clarify these actions:
    Since MATLAB is specific in all its types, and we compile against MATLABs internal blas and lapack libraries,
    we define everything independend from the Python and other languages.
    As such, 64-bit ints and 64-bit floats (i.e., doubles). See ccBlas.h for the translations from MATLABs internal blas
    to the standard version (which uses int* instead of ptrdiff_t).
    */
    #define ADMMint long long int
    #define ADMMfloat double
    #define FMT_INT "d"

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

    extern void dposv(const char *uplo, const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, double *b, const ADMMint *ldb, ADMMint *info);
    extern void dlacpy(const char* uplo, const ADMMint *m, const ADMMint *n, const double* a, const ADMMint *lda, double* b, const ADMMint *ldb);

    #define lapack_posv dposv
    #define lapack_lacpy dlacpy

#else
    #include <cblas.h>
    #include <inttypes.h>
    #define ADMMint blasint

    #ifdef OPENBLAS_USE64BITINT
        #ifdef OPENBLAS_OS_DARWIN
            #define FMT_INT "ld"
        #else
            #define FMT_INT PRId64
        #endif
        
    #else
        #define FMT_INT PRId32
    #endif

    #ifndef USE_SINGLE_PRECISION
        #define ADMMfloat double
        #if INTPTR_MAX == INT64_MAX
            #define cblas_ger scipy_cblas_dger64_
            #define cblas_syr scipy_cblas_dsyr64_
            #define cblas_copy scipy_cblas_dcopy64_
            #define cblas_scal scipy_cblas_dscal64_
            #define cblas_syrk scipy_cblas_dsyrk64_
            #define cblas_gemv scipy_cblas_dgemv64_
            #define cblas_symv scipy_cblas_dsymv64_
            #define cblas_trsv scipy_cblas_dtrsv64_
            #define cblas_amax scipy_cblas_damax64_
            #define cblas_axpy scipy_cblas_daxpy64_
            #define cblas_dot scipy_cblas_ddot64_

            extern void scipy_dposv_64_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, double *b, const ADMMint *ldb, ADMMint *info);
            extern void scipy_dlacpy_64_(const char* uplo, const ADMMint *m, const ADMMint *n, const double* a, const ADMMint *lda, double* b, const ADMMint *ldb);

            #define lapack_posv scipy_dposv_64_
            #define lapack_lacpy scipy_dlacpy_64_
        #elif INTPTR_MAX == INT32_MAX
            #define cblas_ger scipy_cblas_dger
            #define cblas_syr scipy_cblas_dsyr
            #define cblas_copy scipy_cblas_dcopy
            #define cblas_scal scipy_cblas_dscal
            #define cblas_syrk scipy_cblas_dsyrk
            #define cblas_gemv scipy_cblas_dgemv
            #define cblas_symv scipy_cblas_dsymv
            #define cblas_trsv scipy_cblas_dtrsv
            #define cblas_amax scipy_cblas_damax
            #define cblas_axpy scipy_cblas_daxpy
            #define cblas_dot scipy_cblas_ddot

            extern void scipy_dposv_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, double *a, const ADMMint *lda, double *b, const ADMMint *ldb, ADMMint *info);
            extern void scipy_dlacpy_(const char* uplo, const ADMMint *m, const ADMMint *n, const double* a, const ADMMint *lda, double* b, const ADMMint *ldb);

            #define lapack_posv scipy_dposv_
            #define lapack_lacpy scipy_dlacpy_
        #endif
    #else
        //Compile with the flag USE_SINGLE_PRECISION to use 32-bit floats in the solver
        //WARNING: numerical stability with 32-bit floats is significantly worse than using 64-bit floats.
        //         this severely limits the execution speed and precision of the solution. I highly recommend
        //         using 64-bit floats in every scenario possible.
        #define ADMMfloat float
        #if INTPTR_MAX == INT64_MAX
            #define cblas_ger scipy_cblas_sger64_
            #define cblas_copy scipy_cblas_scopy64_
            #define cblas_scal scipy_cblas_sscal64_
            #define cblas_syrk scipy_cblas_ssyrk64_
            #define cblas_gemv scipy_cblas_sgemv64_
            #define cblas_symv scipy_cblas_ssymv64_
            #define cblas_amax scipy_cblas_samax64_
            #define cblas_axpy scipy_cblas_saxpy64_
            #define cblas_dot scipy_cblas_sdot64_

            extern void scipy_sposv_64_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, float *a, const ADMMint *lda, float *b, const ADMMint *ldb, ADMMint *info);
            extern void scipy_slacpy_64_(const char* uplo, const int* m, const int* n, const float* a, const int* lda, float* b, const int* ldb);

            #define lapack_posv scipy_sposv_64_
            #define lapack_lacpy scipy_slacpy_64_
        #elif INTPTR_MAX == INT32_MAX
            #define cblas_ger scipy_cblas_sger
            #define cblas_copy scipy_cblas_scopy
            #define cblas_scal scipy_cblas_sscal
            #define cblas_syrk scipy_cblas_ssyrk
            #define cblas_gemv scipy_cblas_sgemv
            #define cblas_symv scipy_cblas_ssymv
            #define cblas_amax scipy_cblas_samax
            #define cblas_axpy scipy_cblas_saxpy
            #define cblas_dot scipy_cblas_sdot

            extern void scipy_sposv_(const char *uplo, const ADMMint *n, const ADMMint *nrhs, float *a, const ADMMint *lda, float *b, const ADMMint *ldb, ADMMint *info);
            extern void scipy_slacpy_(const char* uplo, const int* m, const int* n, const float* a, const int* lda, float* b, const int* ldb);

            #define lapack_posv scipy_sposv_
            #define lapack_lacpy scipy_slacpy_
        #endif
    #endif
#endif

typedef struct ADMMopts{
    ADMMint verbose;      //silent = 0; iter = 1; timing + memory = 2;
    ADMMint maxIter;      //maximum iterations for the solver
    ADMMfloat sigma;      //slack on x -- ensures posDEF in P+ARA
    ADMMfloat rho_0;      //starting coeff
    ADMMfloat tau;        //exponential Bound decrease
    ADMMfloat alpha;      //exponential increase
    ADMMfloat RBound;     //limit on R
    ADMMfloat eps_abs;    //absolute tolerance
    ADMMfloat eps_rel;    //relative tolerance
    ADMMfloat eps_inf;    //infeasibility tolerance
    ADMMint repInterval;  //reporting interval
    ADMMfloat timeLimit;  //execution time limit, use 0 for unlimited
    ADMMfloat lowRankPer; //percentage required to execute low rank updates instead of full LDL, default 0.05 (5%)
    ADMMint *Pamd;      //preinitialized AMD sorting vector of length nPrim+nDual.
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
#define DEFEPSREL  1e-8
#define DEFEPSINF  1e-8
#define DEFREPIVAL 10
#define DEFTIMELIM 0.0
#define DEFLRPER   0.05

//exitflags
#define EFLAG_MAXITER 0
#define EFLAG_SOLVED 1
#define EFLAG_INACCURATE 2
#define EFLAG_FAILED -1
#define EFLAG_INFEASIBLE -2
#define EFLAG_TIMELIMIT -3
#define EFLAG_NONCONVEX -4

typedef struct ADMMinfo{
    ADMMint nIter;      //number of iterations it took
    ADMMfloat prim_res; //primal convergence score
    ADMMfloat dual_res; //dual convergence score
    ADMMfloat prim_tol; //primal tolerance
    ADMMfloat dual_tol; //dual tolerance
    ADMMfloat runtime;  //run time in seconds
    ADMMfloat objVal;   //The objective value
    char* status;       //the solver status
    ADMMint* Pamd;      //the resulting permutation vector
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