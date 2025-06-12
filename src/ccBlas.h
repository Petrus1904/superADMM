/*
Custom CBLAS header for MATLAB compile. This allows the superADMM solver
to be compiled using MATLABs own BLAS library. Saves the user a download :)
Created by Peter Verheijen, 2025.
*/
// #define CblasUpper 'U'
// #define CblasLower 'L'
// #define CblasTrans 'T'
// #define CblasNoTrans 'N'
// #define CblasColMajor 'C'

#include "superADMM.h"

#define CblasUpper 1
#define CblasLower 2
#define CblasTrans 3
#define CblasNoTrans 4
#define CblasColMajor 5 //its always col major

//BLAS LEVEL 1 Functions in double
void cblas_dcopy(const ADMMint n, const double *x, const ADMMint incx, double *y, const ADMMint incy);
void cblas_dscal(const ADMMint n, const double alpha, double *x, const ADMMint incx);
double cblas_ddot(const ADMMint n, const double *x, const ADMMint incx, const double *y, const ADMMint incy);
void cblas_daxpy(const ADMMint n, const double alpha, const double *x, const ADMMint incx, double *y, const ADMMint incy);
double cblas_damax(const ADMMint n, const double *x, const ADMMint incx);
//BLAS LEVEL 2 Functions in double
void cblas_dgemv(const ADMMint order, const ADMMint trans, const ADMMint m, const ADMMint n, 
                 const double alpha, const double *A, const ADMMint lda, 
                 const double *x, const ADMMint incx, const double beta, 
                 double *y, const ADMMint incy);
void cblas_dsymv(const ADMMint order, const ADMMint uplo, const ADMMint n, 
                 const double alpha, const double *A, const ADMMint lda, 
                 const double *x, const ADMMint incx, const double beta, 
                 double *y, const ADMMint incy);
void cblas_dger(const ADMMint order, const ADMMint m, const ADMMint n, 
                const double alpha, const double *x, const ADMMint incx,
                const double *y, const ADMMint incy, double *A, const ADMMint lda);
//BLAS LEVEL 3 Functions in double
void cblas_dsyrk(const ADMMint order, const ADMMint uplo, const ADMMint trans, const ADMMint n, const ADMMint k,
                 const double alpha, const double *A, const ADMMint lda,
                 const double beta, double *C, const ADMMint ldc);
