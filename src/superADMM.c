/*
    Source file for the superADMM solver
    Created by Peter Verheijen, 2024.
    compile string for python: gcc -shared -o superADMM.dll ../src/superADMM.c ../src/csparse.c ../src/ldl.c -I"C:\OpenBLAS\include" -L"C:\OpenBLAS\lib" -lopenblas -DBUILD_DLL
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <sys\time.h>
#include "superADMM.h"
#include "csparse.h"
#include "ldl.h"
// #include "BK_ldl.h"

//overwrite
#ifdef MATLAB_COMP
    #include "ccBlas.h"
    #include "mex.h" //MATLAB codegen knows what to do :)
    #define print mexPrintf
    #define max(a,b) ((a) > (b) ? (a) : (b))
    #define min(a,b) ((a) < (b) ? (a) : (b))
#else
    #include <cblas.h>
    #include <time.h>
    #include <stdint.h>
    #include <windows.h>
    #define print printf
    #define clock_gettime __pthread_clock_gettime
    // char buf[10000];
    // setvbuf(stdout, buf, _IOFBF, sizeof(buf));
#endif

// Define clock IDs for compatibility
#define CLOCK_REALTIME  0
#define CLOCK_MONOTONIC 1

int setenv(const char *name, const char *value, ADMMint overwrite)
{
    ADMMint errcode = 0;
    return _putenv_s(name, value);
}

ADMMfloat dTime(struct timespec* start){
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC, &end);
    ADMMfloat timeElapsed;
    timeElapsed = (ADMMfloat)(end.tv_sec-start->tv_sec)*1e9;
    timeElapsed = (timeElapsed + (end.tv_nsec - start->tv_nsec))*1e-9;
    return timeElapsed;
}

void printVector(const char* str, ADMMint* x, ADMMint n){
    print("%s = [", str);
    for(ADMMint i = 0; i < n; i ++){
        print("%d ", x[i]);
        if((i+1)%50 == 0){
            print("...\n");
        }
    }
    print("];\n");
}

void printVectorf(const char* str, ADMMfloat* x, ADMMint n){
    print("%s = [", str);
    for(ADMMint i = 0; i < n; i ++){
        if(isnan(x[i])){
            print("nan ");
        } else {
            print("%.7e ", x[i]);
        }
        if((i+1)%50 == 0){
            print("...\n");
        }
    }
    print("];\n");
}


void updatePARA(const ADMMfloat* A, const ADMMfloat* Rup, ADMMfloat* PARA, ADMMint nDual, ADMMint nPrim){
    //Update PARA => PARA + A_i'*Rup_i*A_i, for all nonzero i in Rup.
    //More efficient that full matrix inner product for lower number of nonzeros in Rup.
    for(ADMMint i = 0; i < nDual; i++){
        if(Rup[i] != 0){
            cblas_syr(CblasColMajor, CblasLower, nPrim, Rup[i], &A[i], nDual, PARA, nPrim);
        }
    }
}

void computePARA(const ADMMfloat* P, const ADMMfloat* A, const ADMMfloat* R, const ADMMfloat sigma, ADMMfloat* tmp, ADMMfloat* out, ADMMint nDual, ADMMint nPrim){
    //computes P+A^T*R*A -- we assume column major data --ONLY UPPER TRIANGLE
    //copy
    char uplo = 'L';
    cblas_copy(nDual*nPrim, A, 1, tmp, 1);
    lapack_lacpy(&uplo, &nPrim, &nPrim, P, &nPrim, out, &nPrim);
    // cblas_copy(nPrim*nPrim, P, 1, out, 1);
    //scale tmp = R*A;
    for(ADMMint i = 0; i < nDual; i++){
        cblas_scal(nPrim, sqrt(R[i]), &tmp[i],nDual);
    }
    //cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nPrim, nPrim, nDual, 1.0, tmp, nDual, A, nDual, 1.0, out, nPrim);
    //Symmetric rank k update of full size is cheaper than dgemm
    cblas_syrk(CblasColMajor, CblasLower, CblasTrans, nPrim, nDual, 1.0, tmp, nDual, 1.0, out, nPrim);
    
    if(sigma > 0){ //add sigma on diagonal only if nonzero
        for(ADMMint i = 0; i < nPrim; i++){
            out[i*nPrim + i] += sigma;
        }
    }
}

int chol_lowRankUpdate(ADMMfloat* L, const ADMMfloat* A, const ADMMfloat* Rup, ADMMfloat* tmp, ADMMint nDual, ADMMint nPrim){
    //updates cholesky decomposition A = L*L^T with set of column vectors in A and gains in Rup.
    ADMMint i, j, k;
    ADMMfloat r, c, s, a, lk, xk, rk, rks;
    for(i = 0; i < nDual; i++){
        if(Rup[i] != 0){
            //update
            rk = Rup[i];
            rks = 1;
            cblas_copy(nPrim, &A[i], nDual, tmp, 1);
            for(k = 0; k < nPrim; k++){
                if(tmp[k] != 0){
                    lk = L[k*nPrim+k];
                    xk = tmp[k];
                    r = lk*lk + (rk/rks)*xk*xk;
                    if(r<0){
                        return -k-1; //not posdef
                    }
                    r = sqrt(r);
                    c = r/lk;
                    s = xk/lk;
                    a = xk*rk/(r*rks);
                    L[k*nPrim+k] = r;
                    //SIMD
                    for(j = k+1; j < nPrim-3; j+=4){
                        tmp[j]   -= L[j+k*nPrim]*s;
                        tmp[j+1] -= L[j+k*nPrim+1]*s;
                        tmp[j+2] -= L[j+k*nPrim+2]*s;
                        tmp[j+3] -= L[j+k*nPrim+3]*s;
                        L[j+k*nPrim]   = c*L[j+k*nPrim]   + a*tmp[j];
                        L[j+k*nPrim+1] = c*L[j+k*nPrim+1] + a*tmp[j+1];
                        L[j+k*nPrim+2] = c*L[j+k*nPrim+2] + a*tmp[j+2];
                        L[j+k*nPrim+3] = c*L[j+k*nPrim+3] + a*tmp[j+3];
                        
                    }
                    for(; j < nPrim; j++){
                        tmp[j] -= L[j+k*nPrim]*s;
                        L[j+k*nPrim] = c*L[j+k*nPrim] + a*tmp[j];
                    }
                    rks = rks + rk*(xk*xk)/(lk*lk);
                }
            }
        }
    }
    return 0;
}

cs* cs_copy(const cs* A){
    // creates a copy of A
    cs* out = cs_spalloc(A->m, A->n, A->nzmax, 1, 0);
    for(ADMMint i = 0; i < A->n+1; i++){ //copy column indices
        out->p[i] = A->p[i];
    }
    for(ADMMint j = 0; j<A->nzmax; j++){ //copy the rest
        out->x[j] = A->x[j];
        out->i[j] = A->i[j];
    }
    return out;
}

/* y = AT*x+y */
ADMMint cs_gatxpy (const cs *A, const ADMMfloat *x, ADMMfloat *y){
    // mat-vec multiplication transposed.
    ADMMint p, j, n, *Ap, *Ai;
    ADMMfloat *Ax;
    ADMMfloat val;
    if (!A || !x || !y) return (0);	    /* check inputs */
    n = A->n; 
    Ap = A->p; 
    Ai = A->i; 
    Ax = A->x;
    for (j = 0 ; j < n ; j++)
    {
        val = y[j];
        for (p = Ap[j]; p < Ap[j+1]; p++)
        {
            val += Ax[p] * x[Ai[p]];
        }
        y[j] = val;
    }
    return (1) ;
}

ADMMint cs_syaxpy (const cs *A, const CSfloat *x, CSfloat *y) 
{
    //symmetric y = Ax+y, only upper triangle of A is considered
    ADMMint p, pp, j, n, i, *Ap, *Ai;
    CSfloat *Ax;
    CSfloat val, yval;
    if (!A || !x || !y) return (0) ;	    /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        val = x[j];
        yval = y[j];
        pp = Ap[j+1];
        for (p = Ap [j]; p < pp; p++)
        {
            i = Ai[p];
            if(i< j){
                y[i] += Ax[p] * val;
                yval += Ax[p] * x[i];
            }
            else if(i==j){
                yval += Ax[p] * val;
                //my columns are sorted, I should be able to do this :)
                break;
            } else {
                break; //only if the column has no diagonal entry
            }
        }
        y[j] = yval;
    }
    return (1);
}

ADMMfloat cs_dotCol(const cs* A, ADMMint col, ADMMfloat* x){
    //Computes A_i^Tx with row i.
    ADMMfloat out = 0.0;
    ADMMint pp = A->p[col+1];
    for(ADMMint i = A->p[col]; i < pp; i++){
        out += A->x[i]*x[A->i[i]];
    }
    return out;
}

ADMMint cs_isdiag(const cs* A){
    //checks if sparse matrix is diagonal
    if(A->n != A->nzmax){
        return 0;
    }
    ADMMint n = A->n;
    for(ADMMint i = 0; i < n; i++){
        if(A->p[i] != i || A->i[i] != i){
            return 0;
        }
    }
    return 1;
}

/* free a symbolic factorization */
csldl *cs_ldlfree (csldl *S)
{
    if (!S) return (NULL);	/* do nothing if S already NULL */
    cs_free(S->D);
    cs_free(S->Lp);
    cs_free(S->Li);
    cs_free(S->Lx);
    cs_free(S->P);
    cs_free(S->Pinv);
    cs_free(S->Parent);
    cs_free(S->Lnz);
    if(S->b2) cs_free(S->b2); /* protect if b2 is not defined */
    cs_free(S->Flag);
    cs_free(S->Pattern);
    cs_free(S->Y);
    return (cs_free(S));	/* free the css struct and return NULL */
}
void vec_dset(const ADMMint n, const ADMMfloat alpha, ADMMfloat *x){
    //sets a vector of length n with alpha
    ADMMint i;
    for(i = 0; i < n-3; i+=4){
        x[i] = alpha;
        x[i+1] = alpha;
        x[i+2] = alpha;
        x[i+3] = alpha;
    }
    for(; i<n; i++){
        x[i] = alpha;
    }
}
/*
LDL_rank1update update the sparse LDL decomposition with a rank 1 update, i.e.:
LDL^T + aa*v*v';
but specifically, v only contains 1 nonzero entry. if the nonzero entry of v = 1, aa can be the gain of the update
*/
ADMMint LDL_rank1update(ADMMfloat *Lx,      /* L data */
                        ADMMint *Lp,        /* L column pointer */
                        ADMMint *Li,        /* L row indexes */
                        ADMMfloat *D,       /* D data (assumed diagonal) */
                        ADMMfloat *v,       /* the update vector, contains strictly 1 nonzero */
                        ADMMint starti,     /* the index where v is nonzero */
                        ADMMfloat aa,       /* the gain of the update, aa*v*v' */
                        ADMMint *Flag,      /* Workspace vector, must not have values equal to 1 */
                        ADMMint *queue      /* Workspace vector, undef. on input and output */
                        ){

    ADMMint j, jj, r, r2, top, cc, ccr, ii;
    ADMMfloat p, b, a, Db;
    if(!queue) return 0;
    if(!Flag) return 0;
    queue[0] = starti;
    top = 1;
    a = aa;
    for(jj = 0; jj < top; jj++){
        j = queue[jj];
        p = v[j];
        Flag[j] = 0; //continuation plan
        v[j] = 0; //set to zero to clear array
        Db = D[j] + a*p*p;
        if(Db == 0.0) return -j; //exit -- problem became singular
        b = p*a/Db;
        a = D[j]*a/Db;
        D[j] = Db;
        r2 = Lp[j+1];
        for(r = Lp[j]; r<r2; r++){
            cc = Li[r];
            if(Flag[cc] != 1){
                Flag[cc] = 1;
                queue[top] = cc;
                ii = top;
                //sort -- this can become expensive 
                //     -- but def. cheaper than finding the next one every iteration
                while(ii > 0 && queue[ii-1] > queue[ii]){
                    ccr = queue[ii-1];
                    queue[ii-1] = queue[ii];
                    queue[ii] = ccr;
                    ii--;
                }
                top++;
            }
            v[cc] -= p*Lx[r];
            Lx[r] += b*v[cc];
        }
    }
    return 1;
}

/*
LDL_rankn_solve updates the existing LDL matrices such that 
R = R + sum(Rup_i*Rup_i^T), updates the diagonal entries of R with the new value.
Then also solves the problem LDL^Tx = b
Effective if Rup is only nonzero in a few entries.
*/
ADMMint LDL_rankn_solve(csldl* S,           /* The LDL data package */
                        ADMMfloat *b,       /* the rhs vector, i.e. (LDLT)xout = b*/
                        ADMMfloat *xout,    /* the solution vector, i.e. (LDLT)xout = b*/
                        ADMMfloat *Rup,     /* input, zero on output. A vector containing all is nonzero for all R updates */
                        ADMMfloat *Rtmp,    /* Workspace vector. must be zero on input, zero on output */
                        ADMMint nDual,      /* the length of Rup, Rtmp */
                        ADMMint nPrim       /* nDual + nPrim is the length of b and all stuff in S */
                    ){
    
    ADMMint n = S->n;
    ADMMint eflag = 1;
    ADMMfloat r;
    for(ADMMint i = 0; i<nDual; i++){
        // do rank 1 updates
        if(Rup[i] != 0){
            r = Rup[i];
            Rup[i] = 0; //clear;
            Rtmp[S->Pinv[i+nPrim]] = 1; // set to one, all others are zero, we define the scalar multiplier as the value to update
            eflag = LDL_rank1update(S->Lx, S->Lp, S->Li, S->D, Rtmp, S->Pinv[i+nPrim], r, S->Flag, S->Pattern);
            if(eflag != 1){
                return eflag; //something wrong
            }
        }
    }

    
    ADMMfloat* x = S->Y; //re-use Y here
    ldl_perm(n, x, b, S->P);                /* x = P*b  */
    ldl_lsolve(n, x, S->Lp, S->Li, S->Lx);  /* x = L\x  */
    ldl_dsolve(n, x, S->D);                 /* x = D\x  */
    ldl_ltsolve(n, x, S->Lp, S->Li, S->Lx); /* x = L'\x */
    ldl_permt(n, xout, x, S->P);            /* xout = P'*x */
    
    return nPrim;
}

csldl* LDL_symb(const cs *A, ADMMint order){
    /*
    Wrapper function for the LDL_symbolic from ldl.c. Included inverse permutation from CSparse.
    */
    if(!A) return NULL;
    if(!A->p) return NULL;

    csldl* S = cs_calloc(1, sizeof (csldl));
    if(!S) return NULL; /*failure*/
    
    ADMMint *P;
    ADMMint n = A->n;
    //AMD fails on diagonal matrices. Fortunately, those are the easiest to invert.
    if(!cs_isdiag(A)){
        P = cs_amd(A, order);               /* P = amd(A+A'), or natural */
    } else {
        P = cs_malloc(A->n, sizeof(ADMMint));
        for(ADMMint i = 0; i < A->n; i++){
            P[i] = i; //no permutation
        }
    }
    
    if(!P){
        cs_free(S);
        return NULL;
    } 

    ADMMint* Lp, *Parent, *Flag, *Pinv, *Lnz;
    Lp     = cs_malloc(n+1, sizeof(ADMMint));
    Parent = cs_malloc(n, sizeof(ADMMint));
    Lnz    = cs_malloc(n, sizeof(ADMMint));
    Flag   = cs_malloc(n, sizeof(ADMMint));
    Pinv   = cs_malloc(n, sizeof(ADMMint));
    Lp[n]  = -1; //safeguard check, this should not be 0 after symbolic

    ldl_symbolic(n, A->p, A->i, Lp, Parent, Lnz, Flag, P, Pinv);
    //print("LDL data: Anz=%d, lnz=%d, blkcnt=%d \n", A->p[n], Lp[n], 0);
    ADMMint lnz = Lp[n];
    if(lnz == -1){
        cs_free(Lp);
        cs_free(Parent);
        cs_free(Lnz);
        cs_free(Flag);
        cs_free(Pinv);
        cs_free(P);
        cs_free(S);
        return NULL;
    }
    S->Pinv = Pinv;
    S->P = P;
    S->nnz = lnz;
    S->n = n;
    S->Lp = Lp;
    S->Parent = Parent;
    S->Lnz = Lnz;
    S->Flag = Flag;
    
    S->Lx = cs_malloc(max(1,lnz), sizeof(ADMMfloat));
    S->D = cs_malloc(n, sizeof(ADMMfloat));
    S->Y = cs_malloc(n, sizeof(ADMMfloat));
    S->Li = cs_malloc(max(1,lnz), sizeof(ADMMint));
    S->Pattern = cs_malloc(n, sizeof(ADMMint));

    return S;
}

ADMMint LDLsolve(const cs* A, ADMMfloat *b, ADMMfloat *xout, csldl* S){
    /*
    Wrapper for the ldl.c functions around solving the LDL decomp stuff
    */
    if(!S) return 0; //not OK
    if(!S->Lp) return 0; //not OK

    ADMMint n = A->n;
    ADMMint res;
    res = ldl_numeric(n, A->p, A->i, A->x, 
                      S->Lp, S->Parent, S->Lnz,
                      S->Li, S->Lx, S->D, S->Y, S->Pattern, S->Flag, S->P, S->Pinv);

    if(res >= 0){
        //succes
        ADMMfloat *x = S->Y; // re-use Y here
        ldl_perm(n, x, b, S->P);                /* x = P*b */
        ldl_lsolve(n, x, S->Lp, S->Li, S->Lx);  /* x = L\x */
        ldl_dsolve(n, x, S->D);                 /* x = D\x*/
        ldl_ltsolve(n, x, S->Lp, S->Li, S->Lx); /* x = L'\x*/
        ldl_permt(n, xout, x, S->P);            /* xout = P'*x */
    }

    return res;

}

cs* assemblePARA_cs(const cs* P, const cs* A, const ADMMfloat* R, const ADMMfloat sigma, ADMMint nDual, ADMMint nPrim){
    //constructs PARA = [P+sigma*I, A'; A, -1/R]
    cs* PARA = cs_spalloc(nDual+nPrim, nDual+nPrim, P->nzmax+2*A->nzmax+nDual+nPrim, 1, 0);
    ADMMint k = 0;
    ADMMint i,j,r;
    ADMMfloat *PAx = PARA->x;
    ADMMint *PAp = PARA->p;
    ADMMint *PAi = PARA->i;
    PARA->p[0] = 0;
    for(i = 0; i < nPrim; i++){
        //iterate per column
        ADMMint hasDiag = 0;
        for(j = P->p[i]; j < P->p[i+1]; j++){
            if(P->i[j] > i && hasDiag == 0){
                //col P_i zero on diagonal, include sigma
                PAx[k] = sigma;
                PAi[k] = i;
                k++;
                hasDiag = 1;
            }
            PAx[k] = P->x[j];
            PAi[k] = P->i[j];
            if(P->i[j] == i){
                PAx[k] +=sigma;
                hasDiag = 1;
            }
            k++;
        }
        if(hasDiag == 0){
            //col P_i was strict upper triangular
            PAx[k] = sigma;
            PAi[k] = i;
            k++;
            hasDiag = 1;
        }
        for(j = A->p[i]; j < A->p[i+1]; j++){
            PAx[k] = A->x[j];
            PAi[k] = A->i[j] + nPrim;
            k++;
        }
        PAp[i+1] = k;
    }
    ADMMint *PARAnz = cs_calloc(nDual, sizeof(ADMMint)); //calloc because it must be zero
    for(i = 0; i < A->p[nPrim]; i++){
        PARAnz[A->i[i]]++; //count nonzeros per row
    }
    for(i = 0; i < nDual; i++){
        //cumsum
        r = PAp[nPrim+i]+PARAnz[i];
        PAp[nPrim+i+1] = r+1; //+1 because of the R entry
        PAx[r] = -1/R[i]; //fill in R value
        PAi[r] = nPrim + i;
        PARAnz[i] = 0; //reset counter
    }
    for(i = 0; i < nPrim; i++){
        for(j = A->p[i]; j < A->p[i+1]; j++){
            r = PAp[nPrim+A->i[j]]+PARAnz[A->i[j]];
            PAx[r] = A->x[j];
            PAi[r] = i;
            PARAnz[A->i[j]]++;
        }
    }
    // for(i = 0; i < nDual; i++){
    //     //iterate per column
    //     for(ADMMint j = AT->p[i]; j < AT->p[i+1]; j++){
    //         PARA->x[k] = AT->x[j];
    //         PARA->i[k] = AT->i[j];
    //         k++;
    //     }
    //     PARA->x[k] = -1/R[i];
    //     PARA->i[k] = nPrim + i; //diagonal
    //     k++;
    //     PARA->p[nPrim+i+1] = k;
    // }
    //cs_sprealloc allocates PARA->p with the wrong size...
    cs_free(PARAnz);
    PARA->nzmax = PAp[nPrim+nDual];
    return PARA;
}


void printTime(struct timespec* start){
    ADMMfloat timeElapsed;

    timeElapsed = dTime(start);
    print("%2.6f s\n", timeElapsed);
    clock_gettime(CLOCK_MONOTONIC, start);
}

/* superADMMsolverDense solves quadratic program 
   min_x 0.5x'*P*x + x'*q
   s.t. l<=A*x<=u 
   Considering Dense matrices P and A.
*/
BUILDTAG ADMMint superADMMsolverDense(const ADMMfloat* P, /* (nPrim x nPrim) quadratic cost function matrix*/
                                      const ADMMfloat* q, /* (nPrim) cost function vector */
                                      const ADMMfloat* A, /* (nDual x nPrim) Constraint matrix */
                                      const ADMMfloat* l, /* (nDual) lower bounds vector */
                                      const ADMMfloat* u, /* (nDual) upper bounds vector */
                                      ADMMfloat* x,       /* (nPrim) initial solution guess on input, primal solution on output */
                                      ADMMfloat* y,       /* (nDual) initial dual guess on input, dual solution on output */
                                      ADMMint nPrim,      /* the primal size (i.e., size of cost function) */
                                      ADMMint nDual,      /* the dual size (i.e., size of constraints) */
                                      ADMMopts opts,      /* input struct with solver settings */
                                      ADMMinfo* info      /* output struct with solver information */
                                      ){

    //superADMMsolverDense solves a dense QP problem using the superADMM method
    //returns the exitflag, with meaning:
    //  1: OK
    //  2: requested tolerance not achievable (almost OK)
    //  0: Maximum iterations reached
    // -1: LSmin failed (x update) -- solver failed
    // -2: problem infeasible
    // -3: time limit exceeded
    // WARNING: superADMMsolverDense does not perform any safety checks regarding the sizes of the inputs.
    //          the user should be responsible to ensure they are correct. Interfaces to python and MATLAB
    //          check this for you.
    
    setenv("OPENBLAS_NUM_THREADS", "12", 1);

    if(opts.verbose > 0){
        print("===================== superADMM solver DENSE ====================\n");
    }
    if(opts.verbose == 1){
        print("| Iter | rPrim     | rDual     | rCond     | Bound     | MESSAGE \n");
        print("+------+-----------+-----------+-----------+-----------+---------\n");
    }

    //Timing stuff -- remove before release
    struct timespec tstart, ttotal, tLoopTask;
    ADMMfloat timeLoop1 = 0, timeLoop2 = 0, timeLoop3 = 0, timeLoop4 = 0, timeLoop5 = 0;

    clock_gettime(CLOCK_MONOTONIC, &tstart);
    clock_gettime(CLOCK_MONOTONIC, &ttotal);
    
    const ADMMint infeasInterval = 10;

    //Allocate R
    ADMMfloat* R;
    R = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* z;
    z = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    //Allocate all intermediate helper matrices
    ADMMfloat* PARA;
    PARA = (ADMMfloat*) malloc(nPrim*nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmpA; //temporary matrix of size A
    tmpA = (ADMMfloat*) malloc(nDual*nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmpP; //temporary matrix of size P
    tmpP = (ADMMfloat*) malloc(nPrim*nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmp_q; //temporary vector of size q
    tmp_q = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmp_x; //temporary vector of size x
    tmp_x = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmp_y; //temporary vector of size y
    tmp_y = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* tmp_z; //temporary vector of size y
    tmp_z = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* xp; //previous x.
    xp = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat* Rup;
    Rup = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    //TODO: include prescaling (superADMM often doesnt need it anyway)

    //stuff used for LAPACK FORTRAN calls
    ADMMint one = 1;
    char uplo = 'A'; //all elements
    char uplopos = 'L';
    ADMMfloat zero = 0.0;
    ADMMint LSinfo;
    ADMMint nrhs = 1;
    ADMMint eflag = 1; //OK
    ADMMint repHeadercnt = 0;
    ADMMfloat cond = 0.0;

    //nonconstant values
    ADMMfloat bound = opts.RBound;
    ADMMfloat rPrim = 2*opts.eps_abs;
    ADMMfloat rDual = 2*opts.eps_abs;
    ADMMint Rupdates = 2*nDual;

    //set R and Z
    vec_dset(nDual, opts.rho_0, R);
    //define z0 = A*x0 -- this allows for quick termination if x0 is close to optimal
    cblas_gemv(CblasColMajor, CblasNoTrans, nDual, nPrim, 1.0, A, nDual, x, 1, 0.0, z, 1);

    //PARA = P + A'*R*A;
    computePARA(P, A, R, opts.sigma, tmpA, PARA, nDual, nPrim);

    if(opts.verbose == 2){
        ADMMint probSize = (nPrim*nPrim + nDual*nPrim + 2*nPrim + 3*nDual)*sizeof(ADMMfloat);
        ADMMint workSize = (5*nDual + 3*nPrim + 2*nPrim*nPrim + nDual*nPrim)*sizeof(ADMMfloat);
        print("Problem size in memory: %d (KB)\nSolver workspace memory: %d (KB) \n", probSize/1000, workSize/1000);
        print("Initialization:  ");
        printTime(&tstart);
    }
    //========== END OF SETUP =================
    ADMMint nIter;
    for(nIter = 0; nIter < opts.maxIter; nIter++){

        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
        //compute A'*(R*z-y)-q
        cblas_copy(nPrim, q, 1, tmp_q, 1);
        if(opts.sigma > 0) cblas_axpy(nPrim, -opts.sigma, x, 1, tmp_q, 1); //-sigma because I later do -tmp_q, so it becomes plus then :)
        //tmp_y = R*z-y
        for(ADMMint i = 0; i < nDual; i++){
            tmp_y[i] = R[i]*z[i] - y[i];
        }
        cblas_gemv(CblasColMajor, CblasTrans, nDual, nPrim, 1.0, A, nDual, tmp_y, 1, -1.0, tmp_q, 1);

        timeLoop1 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
        //copy for precision check
        if((nIter+1) % infeasInterval == 0) cblas_copy(nPrim, x, 1, xp, 1); //copy instead
        cblas_copy(nPrim, tmp_q, 1, x, 1);
        

        //x^k+1 = solve using LAPACK -- I use the direct FORTRAN call here since this works with large matrices
        // Since PARA is posdef by definition, we use the dposv function. This uses Cholesky decomposition and more importantly, 
        // only considers the upper triangular part of PARA. Conversevely, we also only compute that part to save more time :).
        // if(Rupdates < opts.lowRankPer*nDual){ //low rank cholesky is slower than plain cholesky... practically regardless of the amount
        //     //tmpP is the cholesky decomp.
        //     LSinfo = chol_lowRankUpdate(tmpP, A, Rup, tmp_x, nDual, nPrim);
        //     cblas_trsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, nPrim, tmpP, nPrim, x, 1);
        //     cblas_trsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, nPrim, tmpP, nPrim, x, 1);
        // } else {
            lapack_lacpy(&uplopos, &nPrim, &nPrim, PARA, &nPrim, tmpP, &nPrim); //copy upper triangular part of PARA for condition check
            lapack_posv(&uplopos, &nPrim, &nrhs, tmpP, &nPrim, x, &nPrim, &LSinfo);
        // }
        
        
        if(LSinfo != 0){ //matrix is not positive definite. Either singular or contains negative eigenvalues (non-convex)
            if(opts.verbose == 1) print("Cholesky failed with info: %d exiting solver\n", LSinfo);
            eflag = -1;
            break; //something went wrong
        }
        timeLoop2 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
        
        //z^k+1 = Ax^k+1
        cblas_gemv(CblasColMajor, CblasNoTrans, nDual, nPrim, 1.0, A, nDual, x, 1, 0.0, z, 1);
        //also update y^k+1
        ADMMfloat zpp;
        for(ADMMint i = 0; i < nDual; i++){
            zpp = z[i] + y[i]/R[i];
            if(zpp > u[i]) zpp = u[i];
            if(zpp < l[i]) zpp = l[i];
            tmp_y[i] = z[i]-zpp; //keep for rPrim check
            y[i] = y[i] + R[i]*tmp_y[i];
            z[i] = zpp;
        }
        timeLoop3 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);

        //do Condition check before convergence -- allows me to re-use temps
        //tmpq = PARA*x-tmpq
        cblas_symv(CblasColMajor, CblasLower, nPrim, 1.0, PARA, nPrim, x, 1, -1.0, tmp_q, 1);
        cond = cblas_amax(nPrim, tmp_q, 1);

        //check convergence
        rPrim = cblas_amax(nDual, tmp_y, 1); //easy :)
        //dual convergence P*x + q + A'y;
        cblas_copy(nPrim, q, 1, tmp_q, 1); //tmpq = q
        cblas_symv(CblasColMajor, CblasUpper, nPrim, 1.0, P, nPrim, x, 1, 1.0, tmp_q, 1); //tmpq = q + P*x
        cblas_gemv(CblasColMajor, CblasTrans, nDual, nPrim, 1.0, A, nDual, y, 1, 1.0, tmp_q, 1); //tmpq = P*x + q+ A'*y;
        rDual = cblas_amax(nPrim, tmp_q, 1);

        if(rPrim < opts.eps_abs && rDual < opts.eps_abs){ //no relative check yet
            //converged
            break;
        }
        if(isnan(rPrim) || isnan(rDual)){
            //protection
            eflag = -1;
            break;
        }
        //check if time limit is exceeded
        if(opts.timeLimit > 0){
            if(dTime(&ttotal) > opts.timeLimit){
                eflag = -3;
                break;
            }
        }

        if((nIter+1) % infeasInterval == 0){ //Look -- if the problem was feasible, I want it to solve fast. if it isnt, I just want to know why
            //primal infeas checks
            //tmp_y = Ax-z, so technically R*tmp_y must be dy.
            for(ADMMint i = 0; i < nDual; i++){
                tmp_y[i] *= R[i];
            }
            cblas_axpy(nPrim, -1.0, x, 1, xp, 1); //xp = -x+xp;
            cblas_scal(nPrim, -1.0, xp, 1); //xp = x-xp;
            ADMMfloat dy_inf = cblas_amax(nDual, tmp_y, 1);
            cblas_gemv(CblasColMajor, CblasTrans, nDual, nPrim, 1.0, A, nDual, tmp_y, 1, 0.0, tmp_q, 1); //tmp_q = A'*dy
            if(cblas_amax(nPrim, tmp_q, 1) < opts.eps_inf*dy_inf){ 
                ADMMfloat primInfeas = 0;
                for(ADMMint i = 0; i < nDual; i++){
                    if(tmp_y[i] > 0){
                        primInfeas += u[i]*tmp_y[i];
                    }
                    else if(tmp_y[i] < 0){
                        primInfeas += l[i]*tmp_y[i];
                    }
                }
                if(primInfeas < 0){
                    //problem is primal infeasible
                    cblas_copy(nPrim, xp, 1, x, 1);
                    cblas_copy(nDual, tmp_y, 1, y,1);
                    eflag = -2;
                    break;
                }
            }

            //dual infeas
            bool isDualInfeas = true;
            ADMMfloat dx_inf = cblas_amax(nPrim, xp, 1);
            if(cblas_dot(nPrim, q, 1, xp, 1) < 0){
                cblas_symv(CblasColMajor, CblasUpper, nPrim, 1.0, P, nPrim, xp, 1, 0.0, tmp_x, 1); //tmpx = P*dx
                if(cblas_amax(nPrim, tmp_x, 1) < opts.eps_inf*dx_inf){
                    //only Adx left
                    for(ADMMint i = 0; i < nDual && isDualInfeas; i++){
                        if(!isinf(l[i]) && !isinf(u[i])){
                            if(fabs(cblas_dot(nPrim, &A[i], nDual, xp, 1)) > opts.eps_inf*dx_inf){
                                isDualInfeas = false;
                                break;
                            }
                        } else if(!isinf(l[i]) && isinf(u[i])) {
                            if(cblas_dot(nPrim, &A[i], nDual, xp, 1) < -opts.eps_inf*dx_inf) {
                                isDualInfeas = false;
                                break;
                            }
                        } else if(isinf(l[i]) && !isinf(u[i])) {
                            if(cblas_dot(nPrim, &A[i], nDual, xp, 1) >  opts.eps_inf*dx_inf) {
                                isDualInfeas = false;
                                break;
                            }
                        }
                    }
                    if(isDualInfeas){
                        cblas_copy(nPrim, xp, 1, x, 1);
                        cblas_copy(nDual, tmp_y, 1, y,1);
                        eflag = -2;
                        break;
                    }
                }
            }
        }//in feas check interval

        timeLoop4 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
        //update R
        Rupdates = 0;
        if(cond - rPrim > -1e-2*cond || cond - rDual > -1e-2*cond){
            bound *= opts.tau;
            Rupdates = nDual;
            if(bound < 1){
                eflag = 2;
                break;
            }
        }
        
        

        ADMMfloat Rt;
        for(ADMMint i = 0; i < nDual; i++){
            Rup[i] = 0;
            Rt = R[i];
            if(z[i] == u[i] || z[i] == l[i]){
                
                if(R[i] < bound){
                    R[i] = min(bound, R[i]*opts.alpha);
                    Rupdates++;
                    Rup[i] = (R[i]) - (Rt);
                } else {
                    R[i] = bound;
                }
                
                // R[i] = min(bound, R[i]*opts.alpha);
            } else {
                if(R[i] > 1/bound){
                    R[i] = max(1/bound, R[i]/opts.alpha); 
                    Rupdates++;
                    Rup[i] = (R[i]) - (Rt);
                } else {
                    R[i] = 1/bound;
                }
                
                // R[i] = max(1/bound, R[i]/opts.alpha); 
            }
            // Rup[i] = (R[i]) - (Rt);
        }
        if(Rupdates < opts.lowRankPer*nDual){
            updatePARA(A, Rup, PARA, nDual, nPrim);
        } else {
            computePARA(P, A, R, opts.sigma, tmpA, PARA, nDual, nPrim);
        }

        if(opts.verbose == 1 && nIter % opts.repInterval == 0){
            print("|  %3d | %.3e | %.3e | %.3e | %.3e | %d\n", nIter, rPrim, rDual, cond, bound, Rupdates);
            repHeadercnt++;
            if(repHeadercnt >= 25){
                print("| Iter | rPrim     | rDual     | rCond     | Bound     | MESSAGE \n");
                repHeadercnt = 0;
            }
        }
        
        timeLoop5 += dTime(&tLoopTask);

    }
    if(nIter == opts.maxIter){
        eflag = 0;
    }
    ADMMfloat objVal = 0;
    cblas_symv(CblasColMajor, CblasUpper, nPrim, 1.0, P, nPrim, x, 1, 0.0, tmp_x, 1);
    objVal = 0.5*cblas_dot(nPrim, tmp_x, 1, x, 1); //0.5*x^T*P*x
    objVal += cblas_dot(nPrim, x, 1, q, 1);
    switch(eflag){
        case 0:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | MAXIMUM ITERATIONS REACHED\n", nIter, rPrim, rDual, cond, bound);
            info->status = "max iteration reached";
            break;
        case 1:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | COMPLETED\n", nIter, rPrim, rDual, cond, bound);
            info->status = "solved";
            break;
        case 2:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | REQUESTED TOLERANCE NOT ACHIEVABLE \n", nIter, rPrim, rDual, cond, bound);
            info->status = "solved inaccurate";
            break;
        case -1:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | SOLVER FAILED \n", nIter, rPrim, rDual, cond, bound);
            info->status = "failed";
            break;
        case -2:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | PROBLEM INFEASIBLE OR UNBOUNDED \n", nIter, rPrim, rDual, cond, bound);
            info->status = "infeasible or unbounded";
            break;
        case -3:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | TIME LIMIT EXCEEDED \n", nIter, rPrim, rDual, cond, bound);
            info->status = "time limit exceeded";
            break;
    }
    if(opts.verbose == 1){
        print("+------+-----------+-----------+-----------+-----------+---------\n");
        
        //print the final stuff
        print("\n");
        print("Iterations:     %d\n", nIter);
        print("Run Time:       ");
        printTime(&tstart);
        print("Objective:      %.4f\n", objVal);
    } 
    if(opts.verbose == 2){
        print("Loop---------------------\n");
        print("--Time KKT:      %2.6f s\n", timeLoop1);
        print("--Time Solve:    %2.6f s\n", timeLoop2);
        print("--Time updatezy: %2.6f s\n", timeLoop3);
        print("--Time cond:     %2.6f s\n", timeLoop4);
        print("--Time PARA:     %2.6f s\n", timeLoop5);
        print("Loop Total:      ");
        printTime(&tstart);
    }
    
    //clean-up
    free(R);
    free(PARA);
    free(tmpA);
    free(tmpP);
    free(z);
    free(tmp_x);
    free(tmp_q);
    free(tmp_y);
    free(tmp_z);
    free(xp);
    free(Rup);

    info->runtime = dTime(&ttotal);
    info->nIter = nIter;
    info->rPrim = rPrim;
    info->rDual = rDual;
    info->objVal = objVal;

    if(opts.verbose == 2){
        print("Total time:      ");
        printTime(&ttotal);
    }

    return eflag;
}


/* superADMMsolverSparse solves quadratic program 
   min_x 0.5x'*P*x + x'*q
   s.t. l<=A*x<=u 
   Considering Sparse (CSC) matrices P and A.
*/
BUILDTAG ADMMint superADMMsolverSparse(ADMMfloat* Pdata,    /* (Pnnz) non-zero values in P */
                                       ADMMint *Pcolptr,    /* (nPrim+1) column pointers of P */
                                       ADMMint *Prowidx,    /* (Pnnz) row indices of P */
                                       const ADMMfloat* q,  /* (nPrim) quadratic cost function vector */
                                       ADMMfloat* Adata,    /* (Annz) non-zero values in A */
                                       ADMMint *Acolptr,    /* (nPrim+1) colum pointers of A */
                                       ADMMint *Arowidx,    /* (Annz) row indices of A */
                                       const ADMMfloat* l,  /* (nDual) vector of lower bounds */
                                       const ADMMfloat* u,  /* (nDual) vector of upper bounds */
                                       ADMMfloat* x,        /* (nPrim) initial primal guess on input, primal solution on output */
                                       ADMMfloat* y,        /* (nDual) initial dual guess on input, dual solution on output */
                                       ADMMint nPrim,       /* number of primal variables (i.e., size of cost function)*/
                                       ADMMint nDual,       /* number of dual variables (i.e., size of constraints) */
                                       ADMMopts opts,       /* input struct with solver settings */
                                       ADMMinfo* info       /* output struct with solver information */
                                    ){
    //superADMMsolverSparse solves a sparse QP problem using the superADMM method
    //returns the exitflag, with meaning:
    //  1: OK
    //  2: requested tolerance not achievable (almost OK)
    //  0: Maximum iterations reached
    // -1: LSmin failed (x update) -- solver failed
    // -2: problem infeasible
    // -3: time limit exceeded
    // -4: problem non-convex
    // WARNING: superADMMsolverSparse does not perform any safety checks regarding the sizes of the inputs.
    //          the user should be responsible to ensure they are correct. Interfaces to python and MATLAB
    //          check this for you.

    setenv("OPENBLAS_NUM_THREADS", "12", 1);

    if(opts.verbose > 0){
        print("==================== superADMM solver SPARSE ====================\n");
    }
    if(opts.verbose == 1){
        print("| Iter | rPrim     | rDual     | rCond     | Bound     | MESSAGE \n");
        print("+------+-----------+-----------+-----------+-----------+---------\n");
    }

    const cs P = {Pcolptr[nPrim], nPrim, nPrim, Pcolptr, Prowidx, Pdata, -1}; //-1 for CSC format
    const cs A = {Acolptr[nPrim], nDual, nPrim, Acolptr, Arowidx, Adata, -1}; //-1 for CSC format

    //Timing stuff -- remove before release
    struct timespec tstart, ttotal, tLoopTask;
    ADMMfloat timeLoop1 = 0, timeLoop2 = 0, timeLoop3 = 0, timeLoop4 = 0, timeLoop5 = 0;

    clock_gettime(CLOCK_MONOTONIC, &tstart);
    clock_gettime(CLOCK_MONOTONIC, &ttotal);

    const ADMMint infeasInterval = 10;
    //Allocate R
    ADMMfloat* R      = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* Rup    = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* Rtmp   = (ADMMfloat*) malloc((nPrim+nDual) * sizeof(ADMMfloat));
    ADMMfloat* z      = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));

    ADMMfloat* tmp_q  = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat* rhs    = (ADMMfloat*) malloc((nPrim + nDual)* sizeof(ADMMfloat));
    ADMMfloat* sol    = (ADMMfloat*) malloc((nPrim + nDual)* sizeof(ADMMfloat));
    ADMMfloat* tmp_x  = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat* tmp_y  = (ADMMfloat*) malloc(nDual * sizeof(ADMMfloat));
    ADMMfloat* xp     = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    ADMMfloat *tmp_yy = (ADMMfloat*) malloc(nDual*sizeof(ADMMfloat));

    //Allocate all intermediate helper matrices
    cs* PARA;
    csldl* S = NULL;

    //TODO: include prescaling (superADMM often doesnt need it anyway)

    //stuff used for LAPACK FORTRAN calls
    ADMMint LSinfo;
    ADMMint eflag = 1; //OK

    ADMMint repHeadercnt = 0; //if 25, we reprint the header.

    //nonconstant values
    ADMMfloat bound = opts.RBound;
    ADMMfloat rPrim = 1.0;
    ADMMfloat rDual = 1.0;
    ADMMfloat cond;

    ADMMint nKKT = nDual+nPrim;

    //set R and Z
    vec_dset(nKKT, 0.0, Rtmp);
    vec_dset(nDual, opts.rho_0, R);
     //Rup is set to zero every iteration by default if alpha isnt 1
    if(opts.alpha == 1) vec_dset(nDual, 0.0, Rup);
    vec_dset(nDual, 0.0, z);
    cs_gaxpy(&A, x, z);
    
    //PARA = [P+sigmaI, A'; A, -R^{-1}];
    PARA = assemblePARA_cs(&P, &A, R, opts.sigma, nDual, nPrim);
    S = LDL_symb(PARA, 0);
    if(!S){
        print("symbolic failed\n");
    }

    if(opts.verbose == 2){
        if(S){
            ADMMint probSize = (P.nzmax + 2*nPrim + A.nzmax + 3*nDual)*sizeof(ADMMfloat) + (A.nzmax + P.nzmax + 2*nPrim + 2)*sizeof(ADMMint);
            ADMMint workSize = (5*nDual + 3*nPrim + 3*(nPrim+nDual)+ S->nnz+2*nKKT+ PARA->nzmax)*sizeof(ADMMfloat) + (8*nKKT+2+ S->nnz+PARA->nzmax)*sizeof(ADMMint);
            print("Problem size in memory: %d (KB)\nSolver workspace memory: %d (KB) \n", probSize/1000, workSize/1000);
        }
        print("Initialization:  ");
        printTime(&tstart);
    }
    //========== END OF SETUP =================
    ADMMint nIter;
    ADMMint Rupdates = 2*nDual;
    for(nIter = 0; nIter < opts.maxIter; nIter++){

        if(!S){
            //exit loop immediately
            eflag = -1;
            break;
        }
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
        //compute [sigma*x -q; z-R^-1*y ]
        //For those raising eyebrows on this weird 4-step for loop
        // somewhere I read that this makes memory loading more efficient as the compiler
        // now loads a whole batch instead of every individual variable.
        // and it enables possible SIMD
        ADMMint jj;
        for(jj = 0; jj < nPrim-3; jj+=4){
            rhs[jj]   = opts.sigma*x[jj]   - q[jj];
            rhs[jj+1] = opts.sigma*x[jj+1] - q[jj+1];
            rhs[jj+2] = opts.sigma*x[jj+2] - q[jj+2];
            rhs[jj+3] = opts.sigma*x[jj+3] - q[jj+3];
        }
        for(; jj < nPrim; jj++){
            rhs[jj] = opts.sigma*x[jj] - q[jj];
        }
        //tmp_y = z-y/R
        for(jj = 0; jj < nDual-3; jj+=4){
            rhs[nPrim+jj]   = z[jj]   - y[jj]/R[jj];
            rhs[nPrim+jj+1] = z[jj+1] - y[jj+1]/R[jj+1];
            rhs[nPrim+jj+2] = z[jj+2] - y[jj+2]/R[jj+2];
            rhs[nPrim+jj+3] = z[jj+3] - y[jj+3]/R[jj+3];
        }
        for(; jj < nDual; jj++){
            rhs[nPrim+jj] = z[jj] - y[jj]/R[jj];
        }
        timeLoop1 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);

        if((nIter+1) % infeasInterval == 0) cblas_copy(nPrim, x, 1, xp, 1); //copy instead

        //x^k+1 = solve
        if(Rupdates < opts.lowRankPer*nDual){
            //use 1 rank updates instead
            LSinfo = LDL_rankn_solve(S, rhs, sol, Rup, Rtmp, nDual, nPrim);
        } else {
            LSinfo = LDLsolve(PARA, rhs, sol, S);
        }
        if(LSinfo < 0){
            if(opts.verbose == 1){
                print("iter: %d, LDL failed with info: %d exiting solver\n", nIter, LSinfo);
                if(opts.sigma == 0){
                    //warn that this can actually cause such a thing
                    print("The sigma option is set to zero, this can cause the solver to fail. Please retry with a positive sigma\n");
                }
            }
            eflag = -1;
            break; //something went wrong
        } else if(LSinfo < nPrim){
            if(opts.verbose == 1){
                print("iter: %d, Problem seems to be non-convex\n", nIter);
            }
            eflag = -4;
            break; //something went wrong
        }

        timeLoop2 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);

        cblas_copy(nPrim, sol, 1, x, 1);
        
        //update z^k+1, y^k+1
        ADMMfloat zpp;
        for(ADMMint i = 0; i < nDual; i++){
            z[i] = z[i] + (sol[nPrim+i] - y[i])/R[i];
            zpp = z[i] + y[i]/R[i];
            if(zpp > u[i]){ zpp = u[i]; }
            else if(zpp < l[i]){ zpp = l[i]; }
            tmp_y[i] = z[i]-zpp; //keep for rPrim check
            y[i] = y[i] + R[i]*tmp_y[i];
            z[i] = zpp;
        }
        timeLoop3 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);

        //do Condition check before convergence -- allows me to re-use temps
        cblas_scal(nPrim+nDual, -1.0, rhs, 1); //rhs = -rhs
        cs_syaxpy(PARA, sol,rhs); //rhs = PARA*sol - rhs
        cond = cblas_amax(nPrim+nDual, rhs, 1);
        

        //check convergence
        rPrim = cblas_amax(nDual, tmp_y, 1); //easy :)
        //dual convergence P*x + q + A'y;
        cblas_copy(nPrim, q, 1, tmp_q, 1); //tmpq = q
        cs_syaxpy(&P, x, tmp_q); //tmpq = q + P*x
        cs_gatxpy(&A, y, tmp_q); //tmpq = q + P*x + AT*y
        rDual = cblas_amax(nPrim, tmp_q, 1);

        if(rPrim < opts.eps_abs && rDual < opts.eps_abs){ //no relative check yet
            //converged
            //Problem -- Since Ax-z is computed differently, the tolerance *could*
            // not be reached yet. Perhaps recompute at this point?
            for(ADMMint i = 0; i < nDual; i++){
                tmp_y[i] = -z[i];
            }
            cs_gaxpy(&A, x, tmp_y);
            rPrim = cblas_amax(nDual, tmp_y, 1);
            if(rPrim < opts.eps_abs){
                break;
            }   
        }
        if(isnan(rPrim) || isnan(rDual)){
            //protection
            eflag = -1;
            break;
        }
         //check if time limit is exceeded -- after the convergence is checked
        if(opts.timeLimit > 0){
            if(dTime(&ttotal) > opts.timeLimit){
                eflag = -3;
                break;
            }
        }
        if((nIter+1) % infeasInterval == 0){ //Look -- if the problem was feasible, I want it to solve fast. if it isnt, I just want to know why
            //primal infeas checks
            //tmp_y = Ax-z, so technically R*tmp_y must be dy.
            for(ADMMint i = 0; i < nDual; i++){
                tmp_y[i] *= R[i];
            }
            cblas_axpy(nPrim, -1.0, x, 1, xp, 1); //xp = -x+xp;
            cblas_scal(nPrim, -1.0, xp, 1);       //xp = x-xp;

            ADMMfloat dy_inf = cblas_amax(nDual, tmp_y, 1);
            vec_dset(nPrim, 0.0, tmp_q);
            cs_gatxpy(&A, tmp_y, tmp_q); //tmp_q = A'*dy
            if(cblas_amax(nPrim, tmp_q, 1) < opts.eps_inf*dy_inf){ 
                ADMMfloat primInfeas = 0;
                for(ADMMint i = 0; i < nDual; i++){
                    if(tmp_y[i] > 0){
                        primInfeas += u[i]*tmp_y[i];
                    }
                    else if(tmp_y[i] < 0){
                        primInfeas += l[i]*tmp_y[i];
                    }
                }
                if(primInfeas < 0){
                    //problem is primal infeasible
                    cblas_copy(nPrim, xp, 1, x, 1);
                    cblas_copy(nDual, tmp_y, 1, y,1);
                    eflag = -2;
                    break;
                }
            }

            //dual infeas
            bool isDualInfeas = true;
            ADMMfloat dx_inf = cblas_amax(nPrim, xp, 1);
            if(cblas_dot(nPrim, q, 1, xp, 1) < 0){
                vec_dset(nPrim, 0.0, tmp_x);
                cs_syaxpy(&P, xp, tmp_x);           //tmpx = P*dx
                if(cblas_amax(nPrim, tmp_x, 1) < opts.eps_inf*dx_inf){
                    //only Adx left
                    vec_dset(nDual, 0.0, tmp_yy);
                    cs_gaxpy(&A, xp, tmp_yy);           //tmp_yy = A*dx
                    for(ADMMint i = 0; i < nDual && isDualInfeas; i++){
                        if(!isinf(l[i]) && !isinf(u[i])){
                            if(fabs(tmp_yy[i]) > opts.eps_inf*dx_inf){
                                isDualInfeas = false;
                                break;
                            }
                        } else if(!isinf(l[i]) && isinf(u[i])) {
                            if(tmp_yy[i] < -opts.eps_inf*dx_inf) {
                                isDualInfeas = false;
                                break;
                            }
                        } else if(isinf(l[i]) && !isinf(u[i])) {
                            if(tmp_yy[i] >  opts.eps_inf*dx_inf) {
                                isDualInfeas = false;
                                break;
                            }
                        }
                    }
                    
                    if(isDualInfeas){
                        //if infeasible, overwrite x and y with the primal and dual certificates
                        cblas_copy(nPrim, xp, 1, x, 1);
                        cblas_copy(nDual, tmp_y, 1, y,1);
                        eflag = -2;
                        break;
                    }
                }
            }
        }//in feas check interval

        //update R
        Rupdates = 0;
        if(cond - rPrim > -1e-2*cond || cond - rDual > -1e-2*cond){ //cond will approach rPrim if rPrim is to machine precision
            bound *= opts.tau;
            Rupdates = nDual;
            if(bound < 1){
                eflag = 2;
                break;
            }
        }
        if(opts.verbose == 1 && nIter % opts.repInterval == 0){
            print("|  %3d | %.3e | %.3e | %.3e | %.3e |\n", nIter, rPrim, rDual, cond, bound);
            repHeadercnt++;
            if(repHeadercnt >= 25){
                print("| Iter | rPrim     | rDual     | rCond     | Bound     | MESSAGE \n");
                repHeadercnt = 0;
            }
        }
        timeLoop4 += dTime(&tLoopTask);
        clock_gettime(CLOCK_MONOTONIC, &tLoopTask);

        ADMMint idx;
        ADMMfloat Rt;
        if(opts.alpha > 1){ //if alpha = 1, dont update, keep Rupdates at zero and benefit from fast LDL
            for(ADMMint i = 0; i < nDual; i++){
                Rup[i] = 0;
                Rt = R[i];
                if(z[i] == u[i] || z[i] == l[i]){
                    if(R[i] < bound){
                        R[i] = min(bound, R[i]*opts.alpha);
                        Rupdates++;
                        Rup[i] = (-1/R[i]) - (-1/Rt);
                    } else {
                        R[i] = bound;
                    }
                } else {
                    if(R[i] > 1/bound){
                        R[i] = max(1/bound, R[i]/opts.alpha); 
                        Rupdates++;
                        Rup[i] = (-1/R[i]) - (-1/Rt);
                    } else {
                        R[i] = 1/bound;
                    }
                }
                idx = PARA->p[nPrim+i+1]-1; //R[i] is always the last column index
                PARA->x[idx] = -1/R[i];
            }
        }
        // print("iter: %d, Rupdates: %d \n", nIter, Rupdates);
        timeLoop5 += dTime(&tLoopTask);

    }
    if(opts.verbose == 2){
        print("Loop---------------------\n");
        print("--Time KKT:      %2.6f s\n", timeLoop1);
        print("--Time Solve:    %2.6f s\n", timeLoop2);
        print("--Time updatezy: %2.6f s\n", timeLoop3);
        print("--Time cond:     %2.6f s\n", timeLoop4);
        print("--Time PARA:     %2.6f s\n", timeLoop5);
        print("Loop Total:      ");
        printTime(&tstart);
    }
    clock_gettime(CLOCK_MONOTONIC, &tLoopTask);
    if(nIter == opts.maxIter){
        eflag = 0;
    }
    ADMMfloat objVal = 0;
    vec_dset(nPrim, 0.0, tmp_x);
    cs_syaxpy(&P, x, tmp_x);
    objVal = 0.5*cblas_dot(nPrim, tmp_x, 1, x, 1); //0.5*x^T*P*x
    objVal += cblas_dot(nPrim, x, 1, q, 1);
    
    switch(eflag){
        case 0:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | MAXIMUM ITERATIONS REACHED\n", nIter, rPrim, rDual, cond, bound);
            info->status = "max iteration reached";
            break;
        case 1:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | COMPLETED\n", nIter, rPrim, rDual, cond, bound);
            info->status = "solved";
            break;
        case 2:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | REQUESTED TOLERANCE NOT ACHIEVABLE \n", nIter, rPrim, rDual, cond, bound);
            info->status = "solved inaccurate";
            break;
        case -1:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | LDL^T FAILED, MATRIX SINGULAR \n", nIter, rPrim, rDual, cond, bound);
            info->status = "failed";
            break;
        case -2:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | PROBLEM INFEASIBLE OR UNBOUNDED \n", nIter, rPrim, rDual, cond, bound);
            info->status = "infeasible or unbounded";
            break;
        case -3:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | TIME LIMIT EXCEEDED \n", nIter, rPrim, rDual, cond, bound);
            info->status = "time limit exceeded";
            break;
        case -4:
            if(opts.verbose == 1) print("|  %3d | %.3e | %.3e | %.3e | %.3e | PROBLEM NON-CONVEX \n", nIter, rPrim, rDual, cond, bound);
            info->status = "problem non-convex";
            break;
    }
    if(opts.verbose == 1){
        print("+------+-----------+-----------+-----------+-----------+---------\n");
        
        //print the final stuff
        print("\n");
        print("Iterations:     %d\n", nIter);
        print("Run Time:       ");
        printTime(&tstart);
        print("Objective:      %.4f\n", objVal);
    } 

    //clean-up
    free(R);
    free(Rup);
    free(Rtmp);
    cs_spfree(PARA);
    cs_ldlfree(S);
    free(z);
    free(tmp_x);
    free(tmp_q);
    free(tmp_y);
    free(tmp_yy);
    free(xp);
    free(rhs);
    free(sol);

    info->runtime = dTime(&ttotal);
    info->nIter = nIter;
    info->rPrim = rPrim;
    info->rDual = rDual;
    info->objVal = objVal;

    if(opts.verbose == 2){
        print("Termination:     ");
        printTime(&tLoopTask);
        print("Total time:      ");
        printTime(&ttotal);
        print("\n");
    }

    return eflag;
}

/*
int LDLtest(ADMMfloat* Pdata, ADMMint *Prowptr, ADMMint *Pcolidx, ADMMint Pnnz, ADMMint nPrim, ADMMfloat* q, ADMMfloat *x){
    ADMMint LSi = 0;
    css* Si;
    const cs Pi = {Pnnz, nPrim, nPrim, Prowptr, Pcolidx, Pdata, -1};
    Si = LDL_symb(&Pi, 0);
    ADMMfloat* tmp_xi; //temporary vector of size x
    tmp_xi = (ADMMfloat*) malloc(nPrim * sizeof(ADMMfloat));
    cblas_copy(nPrim, q, 1, x, 1);
    struct timespec tstar;
    clock_gettime(CLOCK_MONOTONIC, &tstar);
    LSi = LDLsolve(&Pi, x, Si);
    print("native LDL took:  ");
    printTime(&tstar);
    cblas_copy(nPrim, q, 1, tmp_xi, 1);
    cblas_scal(nPrim, -1.0, tmp_xi, 1);
    cs_gaxpy(&Pi, x, tmp_xi); //tmpq = q + P*x
    double rBK = cblas_amax(nPrim, tmp_xi, 1);
    print("LDL err = %.9e \n", rBK);

    Si = LDL_symb(&Pi, 0);
    cblas_copy(nPrim, q, 1, x, 1);
    clock_gettime(CLOCK_MONOTONIC, &tstar);
    LSi = BK_LDLsolve(&Pi, x, Si);
    print("BK LDL took:  ");
    printTime(&tstar);
    cblas_copy(nPrim, q, 1, tmp_xi, 1);
    cblas_scal(nPrim, -1.0, tmp_xi, 1);
    cs_gaxpy(&Pi, x, tmp_xi); //tmpq = q + P*x
    rBK = cblas_amax(nPrim, tmp_xi, 1);
    print("BK err = %.9e\n", rBK);

    cs_sfree(Si);
    free(tmp_xi);
    print("Error = %d \n", LSi);

    return 0;
}
*/

int main() {

    printf("Hello World");
    return 0;
}



