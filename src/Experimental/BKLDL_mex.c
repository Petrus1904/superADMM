/*Matlab interface for the superADMM solver*/

#include "mex.h"
#include "BK_ldl.h"
#include "csparse.h"
#include "ldl.h"
#include "matrix.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define abs(a) (a < 0 ? -a : a)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /*input args are P (sparse matrix), q (vector), A (sparse matrix), l, u
    optional -- x, y, opts*/
    //the function always returns one value: x, even if no return values are expected.
    // MATLAB will handle this by pushing this to "ans" in the workspace. It also allows usage in the command window
    // where it simply prints the solution. 
    const mwSize *dims_A, *dims_b;
    //do error checking
    if(nrhs < 2){
        mexErrMsgIdAndTxt("BKLDL:matchdims", "not enough inputs, please specify A and b");
    }

    dims_A = mxGetDimensions(prhs[0]);
    dims_b = mxGetDimensions(prhs[1]);

    if(!mxIsSparse(prhs[0])){
        mexErrMsgIdAndTxt("BKLDL:invalidtype", "A must be sparse");
    }
    if(dims_A[0] != dims_A[1]){
        mexErrMsgIdAndTxt("BKLDL:sizemismatch", "A has an invalid size, must be (n x n) matrix");
    }
    if(dims_b[1] != 1){ //mxGetDimensions always returns 2 dimensions, since arrays dont exist in MATLAB
        mexErrMsgIdAndTxt("BKLDL:sizemismatch", "b must be a column vector of (n x 1)");
    }
    if(dims_b[0] != dims_A[0]){
        mexErrMsgIdAndTxt("superADMM:matchdims", "b must have the same number of rows as A");
    }
    
    ADMMint n = (ADMMint)dims_A[0];
    const ADMMfloat *b = (ADMMfloat*)mxGetDoubles(prhs[1]);
    ADMMfloat *x;
    ADMMint eflag;

    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    x = (ADMMfloat*)mxGetDoubles(plhs[0]);

    //============================== EXECUTION =============================

    //extract A
    ADMMfloat *Ax = (ADMMfloat*)mxGetDoubles(prhs[0]);
    ADMMint *Ap = (ADMMint*)mxGetJc(prhs[0]);
    ADMMint *Ai = (ADMMint*)mxGetIr(prhs[0]);

    ADMMint *P0 = NULL;
    if(nrhs == 3){
        P0 = (ADMMint*)mxGetInt64s(prhs[2]);
    }
    
    // eflag = ExportBKsolve(Ax, Ap, Ai, b, x, n, P0);
    // if(eflag != 1){
    //     mexErrMsgIdAndTxt("BKLDL:failed", "solver failed, exit code: %d", eflag);
    // }

    ADMMint blkcnt = 0;
    ADMMint *P;
    if(!P0){
        const cs A = {Ap[n],n,n,Ap,Ai,Ax,-1};
        P = cs_amd(&A, 0);
    } else {
        P = cs_malloc(n, sizeof(ADMMint));
        for(ADMMint i = 0; i < n; i++){
            P[i] = P0[i];
        }
    }
    ADMMint *Pinv = cs_pinv(P, n);

    ADMMint *b2     = (ADMMint*)malloc((n+1)*sizeof(ADMMint));
    ADMMint *Lp     = (ADMMint*)malloc((n+1)*sizeof(ADMMint));
    ADMMint *Parent = (ADMMint*)malloc((n)*sizeof(ADMMint));
    ADMMint *Lnz    = (ADMMint*)malloc((n)*sizeof(ADMMint));
    ADMMint *Flag   = (ADMMint*)malloc((n)*sizeof(ADMMint));

    blkcnt = bk_permute(n, Ax, Ai, Ap, P, Pinv, b2);
    bk_symbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag, b2, P, Pinv);

    ADMMint lnz = Lp[n];
    ADMMfloat *Lx   = (ADMMfloat*) malloc(lnz*sizeof(ADMMfloat));
    ADMMfloat *D    = (ADMMfloat*) malloc(b2[n]*sizeof(ADMMfloat));
    ADMMint *Li     = (ADMMint*) malloc(lnz*sizeof(ADMMint));
    ADMMfloat *Y    = (ADMMfloat*) malloc(2*n*sizeof(ADMMfloat));
    ADMMint *Pattern = (ADMMint*) malloc(n*sizeof(ADMMint));
    memset(Li, -1, lnz * sizeof(ADMMint));

    ADMMint res = bk_numeric(n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern, Flag, b2, P, Pinv);

    ADMMint reallnz = 0;
    ADMMint cntzero = 0;
    for(int i =0; i < n; i++){
        reallnz += Lnz[i]-Lp[i];
        for(int j = Lp[i]; j < Lnz[i]; j++){
            if(abs(Lx[j]) > 0){
                cntzero ++;
            }
        }
    }
    mexPrintf("BK status: nblks: %d, proj_nnz: %d, nnz: %d, actual nnz: %d \n", blkcnt, lnz, reallnz, cntzero);
    if(res == n){
        //succes
        ADMMfloat *xx = (ADMMfloat*)malloc(n*sizeof(ADMMfloat));
        ldl_perm(n, xx, b, P);           /* x = P*b */
        bk_lsolve(n, xx, Lp, Li, Lx, Lnz); /* x = L\x */
        bk_dsolve(n, xx, D, b2);            /* x = D\x*/
        bk_ltsolve(n, xx, Lp, Li, Lx, Lnz);/* x = L'\x*/
        ldl_permt(n, x, xx, P);          /* b = P'*x */
        free(xx);
    }

    if(nlhs == 4){
        //return L, D, P too
        plhs[1] = mxCreateSparse(n, n, lnz, mxREAL);
        ADMMfloat *Lex = (ADMMfloat*)mxGetDoubles(plhs[1]);
        ADMMint *Lep = (ADMMint*)mxGetJc(plhs[1]);
        ADMMint *Lei = (ADMMint*)mxGetIr(plhs[1]);

        //copy
        ADMMint jj = 0;
        Lep[0] = 0;
        for(ADMMint i = 0; i < n; i++){
            ADMMint st = Lp[i];
            ADMMint ed = Lnz[i];
            for(ADMMint ii = st; ii < ed; ii++){
                Lex[jj] = Lx[ii];
                Lei[jj] = Li[ii];
                jj++;
            }
            Lep[i+1] = jj;
        }
        //D
        plhs[2] = mxCreateSparse(n, n, n+(b2[n]-n)*2, mxREAL);
        ADMMfloat *Dex = (ADMMfloat*)mxGetDoubles(plhs[2]);
        ADMMint *Dep = (ADMMint*)mxGetJc(plhs[2]);
        ADMMint *Dei = (ADMMint*)mxGetIr(plhs[2]);

        //copy
        jj = 0;
        Lep[0] = 0;
        ADMMint cc = 0;
        ADMMint k = 0;
        for(ADMMint i = 0; i < n; i++){
            jj = b2[i];
            if(cc == 0){ //this is much cheaper to compute
                if(jj == b2[i+1]){
                    Dex[k] = D[jj];
                    Dex[k+1] = D[jj+1];
                    Dex[k+2] = D[jj+1];
                    Dex[k+3] = D[jj+2];

                    Dei[k] = i;
                    Dei[k+1] = i+1;
                    Dei[k+2] = i;
                    Dei[k+3] = i+1;

                    Dep[i+1] = k+2;
                    Dep[i+2] = k+4;
                    k = k + 4;
                    
                    cc = 1; //this skips the next outer loop
                } else {
                    Dex[k] = D[jj];
                    Dei[k] = i;
                    k++;
                    Dep[i+1] = k;
                }
            } else {
                cc = 0;
            }
        }

        plhs[3] = mxCreateNumericMatrix(n, 1, mxINT64_CLASS, mxREAL);
        ADMMint* P1 = (ADMMint*)mxGetInt64s(plhs[3]);
        for(ADMMint i = 0; i < n; i++){
            P1[i] = P[i];
        }
    }

    cs_free(P);
    cs_free(Pinv);
    free(b2);
    free(Lp);
    free(Parent);
    free(Lnz);
    free(Flag);
    free(Lx);
    free(D);
    free(Li);
    free(Y);
    free(Pattern);
}


