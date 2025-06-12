/*
Custom CBLAS wrapper for MATLAB compile. This allows the superADMM solver
to be compiled using MATLABs own BLAS library. Saves the user a download :)
Created by Peter Verheijen, 2025.
*/
#include "ccBlas.h"
#include "blas.h" //matlab knows

#define abs(a) ((a) < 0 ? (-a) : (a)) // a < 0 ensures that zero is always seen as "positive"

//SIDENOTE we will always cast to ptrdiff_t, to be sure. even if the presented int type is the same

//BLAS LEVEL 1 Functions in double
void cblas_dcopy(const ADMMint n, const double *x, const ADMMint incx, double *y, const ADMMint incy){
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    dcopy(&pn, x, &pincx, y, &pincy);
}
void cblas_dscal(const ADMMint n, const double alpha, double *x, const ADMMint incx){
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    dscal(&pn, &alpha, x, &pincx);
}
double cblas_ddot(const ADMMint n, const double *x, const ADMMint incx, const double *y, const ADMMint incy){
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    return ddot(&pn, x, &pincx, y, &pincy);
}
void cblas_daxpy(const ADMMint n, const double alpha, const double *x, const ADMMint incx, double *y, const ADMMint incy){
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    daxpy(&pn, &alpha, x, &pincx, y, &pincy);
}
double cblas_damax(const ADMMint n, const double *x, const ADMMint incx){
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t idx = idamax(&pn, x, &pincx);
    return idx ? abs(x[idx-1]) : -1; //fortran is 1 based indexing, lets pick -1 as error
}
//BLAS LEVEL 2 Functions in double
void cblas_dgemv(const ADMMint order, const ADMMint trans, const ADMMint m, const ADMMint n, 
                 const double alpha, const double *A, const ADMMint lda, 
                 const double *x, const ADMMint incx, const double beta, 
                 double *y, const ADMMint incy){
    char t = 'N';
    if(trans == CblasTrans){
        t = 'T';
    }
    ptrdiff_t pm = (ptrdiff_t)m;
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t plda = (ptrdiff_t)lda;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    dgemv(&t, &pm, &pn, &alpha, A, &plda, x, &pincx, &beta, y, &pincy);                
}
void cblas_dsymv(const ADMMint order, const ADMMint uplo, const ADMMint n, 
                 const double alpha, const double *A, const ADMMint lda, 
                 const double *x, const ADMMint incx, const double beta, 
                 double *y, const ADMMint incy){
    char u = 'U';
    if(uplo == CblasLower){
        u = 'L';
    }
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t plda = (ptrdiff_t)lda;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    dsymv(&u, &pn, &alpha, A, &plda, x, &pincx, &beta, y, &pincy);                
}
void cblas_dger(const ADMMint order, const ADMMint m, const ADMMint n, 
                const double alpha, const double *x, const ADMMint incx,
                const double *y, const ADMMint incy, double *A, const ADMMint lda){
    ptrdiff_t pm = (ptrdiff_t)m;
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t plda = (ptrdiff_t)lda;
    ptrdiff_t pincx = (ptrdiff_t)incx;
    ptrdiff_t pincy = (ptrdiff_t)incy;
    dger(&pm, &pn, &alpha, x, &pincx, y, &pincy, A, &plda);                
}

//BLAS LEVEL 3 Functions in double
void cblas_dsyrk(const ADMMint order, const ADMMint uplo, const ADMMint trans, const ADMMint n, const ADMMint k,
                 const double alpha, const double *A, const ADMMint lda,
                 const double beta, double *C, const ADMMint ldc){
    char t = 'N';
    char u = 'U';
    if(trans == CblasTrans){
        t = 'T';
    }
    if(uplo == CblasLower){
        u = 'L';
    }
    ptrdiff_t pn = (ptrdiff_t)n;
    ptrdiff_t pk = (ptrdiff_t)k;
    ptrdiff_t plda = (ptrdiff_t)lda;
    ptrdiff_t pldc = (ptrdiff_t)ldc;
    dsyrk(&u, &t, &pn, &pk, &alpha, A, &plda, &beta, C, &pldc);                
}