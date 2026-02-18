/*Matlab interface for the superADMM solver*/

#include "mex.h"
#include "superADMM.h" //for ADMMfloat and ADMMint

ADMMint* mwIndextoInt(mwIndex* in, ADMMint len){
    ADMMint* out = malloc(len*sizeof(ADMMint));
    for(int i = 0; i < len; i++){
        out[i] = (ADMMint)in[i];
    }
    return out;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /*input args are P (sparse matrix), q (vector), A (sparse matrix), l, u
    optional -- x, y, opts*/
    //the function always returns one value: x, even if no return values are expected.
    // MATLAB will handle this by pushing this to "ans" in the workspace. It also allows usage in the command window
    // where it simply prints the solution. 
    const mwSize *dims_P, *dims_A, *dims_q, *dims_l, *dims_u;
    //do error checking
    if(nrhs < 5){
        mexErrMsgIdAndTxt("superADMM:matchdims", "not enough inputs, please specify P, q, A, l, u");
    }

    dims_P = mxGetDimensions(prhs[0]);
    dims_q = mxGetDimensions(prhs[1]);
    dims_A = mxGetDimensions(prhs[2]);
    dims_l = mxGetDimensions(prhs[3]);
    dims_u = mxGetDimensions(prhs[4]);
    
    //=====================================START OF ERROR CHECKING==============================
    for(int i = 0; i < 5; i++){
        if(mxGetNumberOfDimensions(prhs[i]) > 2){
            mexErrMsgIdAndTxt("superADMM:sizemismatch", "tensor or higher dimensional matrices cannot be inputs");
        }
        if(!mxIsDouble(prhs[i])){
            mexErrMsgIdAndTxt("superADMM:invalidtype", "P, q, A, l and u must be of doubles (can be sparse)");
        }
        if(mxIsComplex(prhs[i])){
            mexErrMsgIdAndTxt("superADMM:invalidtype", "All inputs must be real");
        }
    }
    if(dims_P[0] != dims_P[1]){
        mexErrMsgIdAndTxt("superADMM:sizemismatch", "P has an invalid size, must be (n x n) matrix");
    }
    if(dims_q[1] != 1){ //mxGetDimensions always returns 2 dimensions, since arrays dont exist in MATLAB
        mexErrMsgIdAndTxt("superADMM:sizemismatch", "q must be a column vector of (n x 1)");
    }
    if(dims_q[0] != dims_P[0]){
        mexErrMsgIdAndTxt("superADMM:matchdims", "P and q mismatch in size");
    }
    if(mxIsSparse(prhs[1])){
        mexErrMsgIdAndTxt("superADMM:typemismatch", "q must be a dense vector");
    }
    if(dims_l[0] != dims_u[0]){
        mexErrMsgIdAndTxt("superADMM:matchdims", "l and u do not have the same length");
    }
    if(dims_A[0] != dims_l[0]){
        mexErrMsgIdAndTxt("superADMM:matchdims", "A must have the same number of rows as in l and u");
    }
    if(dims_A[1] != dims_P[0]){
        mexErrMsgIdAndTxt("superADMM:matchdims", "A must have the same number of columns as P");
    }
    if(dims_l[1] != 1){
        mexErrMsgIdAndTxt("superADMM:sizemismatch", "l must be a column vector of (m x 1)");
    }
    if(dims_u[1] != 1){
        mexErrMsgIdAndTxt("superADMM:sizemismatch", "u must be a column vector of (m x 1)");
    }
    if(mxIsSparse(prhs[3])){
        mexErrMsgIdAndTxt("superADMM:typemismatch", "l must be a dense vector");
    }
    if(mxIsSparse(prhs[4])){
        mexErrMsgIdAndTxt("superADMM:typemismatch", "u must be a dense vector");
    }
    if(nrhs == 8 && !mxIsStruct(prhs[7])){
        mexErrMsgIdAndTxt("superADMM:typemismatch", "options argument must be a struct");
    }
    //verify that P and A are both sparse or dense
    if(mxIsSparse(prhs[0]) != mxIsSparse(prhs[2])){ //XOR condition
        mexErrMsgIdAndTxt("superADMM:typemismatch", "P and A must be both sparse or dense");
    }
    
    ADMMint nPrim = (ADMMint)dims_P[0];
    ADMMint nDual = (ADMMint)dims_l[0];

    const ADMMfloat *q = (ADMMfloat*)mxGetDoubles(prhs[1]);
    const ADMMfloat *l = (ADMMfloat*)mxGetDoubles(prhs[3]);
    const ADMMfloat *u = (ADMMfloat*)mxGetDoubles(prhs[4]);

    ADMMfloat *x;
    ADMMfloat *y;
    mxArray *mxPamd;

    ADMMint eflag;

    plhs[0] = mxCreateDoubleMatrix(nPrim, 1, mxREAL);
    x = (ADMMfloat*)mxGetDoubles(plhs[0]);

    //define options
    ADMMopts opts = {DEFVERBOSE, DEFMAXITER, DEFSIGMA,
                     DEFRHO0, DEFTAU, DEFALPHA,
                     DEFRBOUND, DEFEPSABS, DEFEPSREL, DEFEPSINF,
                     DEFREPIVAL, DEFTIMELIM, DEFLRPER, NULL};

    ADMMint copyPamd = 0;

    if(nrhs == 8){
        //overwrite options
        mxArray* opt;
        opt = mxGetField(prhs[7], 0, "verbose");
        if(opt && mxIsScalar(opt)) opts.verbose     = (ADMMint)mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "maxIter");
        if(opt && mxIsScalar(opt)) opts.maxIter     = (ADMMint)mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "sigma");
        if(opt && mxIsScalar(opt)) opts.sigma       = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "rho_0");
        if(opt && mxIsScalar(opt)) opts.rho_0       = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "tau");
        if(opt && mxIsScalar(opt)) opts.tau         = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "alpha");
        if(opt && mxIsScalar(opt)) opts.alpha       = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "RBound");
        if(opt && mxIsScalar(opt)) opts.RBound      = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "eps_abs");
        if(opt && mxIsScalar(opt)) opts.eps_abs     = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "eps_rel");
        if(opt && mxIsScalar(opt)) opts.eps_rel     = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "eps_inf");
        if(opt && mxIsScalar(opt)) opts.eps_inf     = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "repInterval");
        if(opt && mxIsScalar(opt)) opts.repInterval = (ADMMint)mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "timeLimit");
        if(opt && mxIsScalar(opt)) opts.timeLimit   = mxGetScalar(opt);
        opt = mxGetField(prhs[7], 0, "lowRankPer");
        if(opt && mxIsScalar(opt)){
            opts.lowRankPer  = mxGetScalar(opt);
        } else {
            if(!mxIsSparse(prhs[0])){
                //increase this value for dense matrices
                //appears to significantly improve performance
                opts.lowRankPer = 0.5;
            }
        }
        opt = mxGetField(prhs[7], 0, "Pamd");
        if(opt && !mxIsEmpty(opt)){
            //verify length
            if(!mxIsInt64(opt)){
                mexErrMsgIdAndTxt("superADMM:typemismatch", "Permutation vector Pamd must be of type Int64. Please use 'info.Pamd' for the correct permutation. Type 'help superADMM_preCompute' for more information.");
            } else {
                const mwSize *dims_Pamd = mxGetDimensions(opt);
                if(dims_Pamd[0] != nPrim+nDual+1) mexErrMsgIdAndTxt("superADMM:sizemismatch", "Permutation vector Pamd is of incorrect length, please keep the problem size unchanged when re-using Pamd. Type 'help superADMM_preCompute' for more information.");
                copyPamd = 1; //flag to copy the Pamd after all the error checks
            }
        }
    }
    //check options
    if(opts.sigma < 0) mexErrMsgIdAndTxt("superADMM:valueError", "sigma must be nonnegative");
    if(opts.alpha < 1) mexErrMsgIdAndTxt("superADMM:valueError", "alpha must be greater or equal to 1");
    if(opts.alpha == 1) mexWarnMsgIdAndTxt("superADMM:valueWarning", "alpha = 1 disables the superlinear convergence rate, set alpha > 1 for fast convergence");
    if(opts.rho_0 <= 0) mexErrMsgIdAndTxt("superADMM:valueError", "rho_0 must be greater than zero");
    if(opts.tau > 1 || opts.tau <= 0) mexErrMsgIdAndTxt("superADMM:valueError", "tau must be between 0 and 1");
    if(opts.tau == 1) mexWarnMsgIdAndTxt("superADMM:valueWarning", "tau = 1 disables the numerical bounding method, this can lead to nonconverging effects");
    if(opts.RBound < 1) mexErrMsgIdAndTxt("superADMM:valueError", "RBound must be greater than 1");
    if(opts.eps_abs < 0) mexErrMsgIdAndTxt("superADMM:valueError", "absolute convergence parameter eps_abs must be greater or equal to 0");
    if(opts.eps_rel < 0) mexErrMsgIdAndTxt("superADMM:valueError", "relative convergence parameter eps_rel must be greater or equal to 0");
    if(opts.eps_abs == 0 && opts.eps_rel == 0) mexWarnMsgIdAndTxt("superADMM:valueWarning", "setting both convergence criteria to 0 prevents exiting with exitflag 1 (Solved correctly), please check the solution tolerances to ensure a satisfactory solution");
    if(opts.eps_inf <= 0) mexErrMsgIdAndTxt("superADMM:valueError", "infeasibility tolerance eps_inf must be greater than 0");
    if(opts.repInterval < 1) mexErrMsgIdAndTxt("superADMM:valueError", "repInterVal must be greater or equal to 1");
    if(opts.lowRankPer > 1 || opts.lowRankPer < 0) mexErrMsgIdAndTxt("superADMM:valueError", "Low rank update percentage (lowRankPer) must be between 0 and 1 (inclusive)");
    
    //basic verification of u and l --> early detect infeasibilities
    for(int i = 0; i < nDual; i++){
        if(u[i] < l[i]){
            mexErrMsgIdAndTxt("superADMM:infeasible", "constraint at index %d is infeasible (%.3f > %.3f)", i+1, l[i], u[i]);
        }
    }

    //check if x0, y0 exists, and initialize x and y.
    if(nrhs >= 6){
        if(!mxIsEmpty(prhs[5])){
            const mwSize *dims_x = mxGetDimensions(prhs[5]);
            if(dims_x[0] != nPrim){
                mexErrMsgIdAndTxt("superADMM:matchdims", "x must be vector of size (n x 1)");
            }
            if(dims_x[1] != 1){
                mexErrMsgIdAndTxt("superADMM:matchdims", "x must be vector of size (n x 1)");
            }
            if(!mxIsDouble(prhs[5]) || mxIsSparse(prhs[5])){
                mexErrMsgIdAndTxt("superADMM:typemismatch", "x must be a dense double vector");
            }
            if(mxIsComplex(prhs[5])){
                mexErrMsgIdAndTxt("superADMM:invalidtype", "All inputs must be real");
            }
            const ADMMfloat *x0 = (ADMMfloat*)mxGetDoubles(prhs[5]);
            for(int i = 0; i < nPrim; i++){
                //copy
                x[i] = x0[i];
            }
        }
    }
    if(nrhs >= 7){
        if(!mxIsEmpty(prhs[6])){
            const mwSize *dims_y = mxGetDimensions(prhs[6]);
            if(dims_y[0] != nDual){
                mexErrMsgIdAndTxt("superADMM:matchdims", "y must have the same number of elements as l or u");
            }
            if(dims_y[1] != 1){
                mexErrMsgIdAndTxt("superADMM:matchdims", "y must be vector of size (m x 1)");
            }
            if(!mxIsDouble(prhs[6]) || mxIsSparse(prhs[6])){
                mexErrMsgIdAndTxt("superADMM:typemismatch", "y must be a dense double vector");
            }
            if(mxIsComplex(prhs[6])){
                mexErrMsgIdAndTxt("superADMM:invalidtype", "All inputs must be real");
            }
        }
    }

    //===============================END OF ERROR CHECKING=================================

    if(nlhs > 1){
        plhs[1] = mxCreateDoubleMatrix(nDual, 1, mxREAL);
        y = (ADMMfloat*)mxGetDoubles(plhs[1]);
    } else {
        y = (ADMMfloat*)calloc(nDual, sizeof(ADMMfloat)); //calloc here, free later
    }
    if(nrhs >= 7){
        if(!mxIsEmpty(prhs[6])){
            const ADMMfloat *y0 = (ADMMfloat*)mxGetDoubles(prhs[6]);
            for(int i = 0; i < nDual; i++){
                //copy
                y[i] = y0[i];
            }
        }
    }

    ADMMint *Pamd = NULL;
    if(!mxIsSparse(prhs[0])){
        //replace Pamd with an empty matrix -- prevent errors on re-runs.
        if(nlhs == 4){
            mxPamd = mxCreateNumericMatrix(0, 0, mxINT64_CLASS, mxREAL);
        }
    } else {
        if(nlhs == 4){
            //only do this if the output is requested!
            mxPamd = mxCreateNumericMatrix(nDual+nPrim+1, 1, mxINT64_CLASS, mxREAL);
            Pamd = (ADMMint*)mxGetInt64s(mxPamd);
        } else {
            Pamd = (ADMMint*)malloc((nDual+nPrim+1)*sizeof(ADMMint)); //free later
        }
        Pamd[0] = -1; //flag to indicate data is unusable
        if(copyPamd == 1){
            mxArray* opt = mxGetField(prhs[7], 0, "Pamd");
            const ADMMint* Pamd0 = (ADMMint*)mxGetInt64s(opt);
            for(int i = 0; i < nDual + nPrim; i++){
                //numerical checks
                if(Pamd0[i] >= nDual + nPrim || Pamd0[i] < 0){
                    mexErrMsgIdAndTxt("superADMM:valueError", "Permutation vector Pamd incorrectly defined.");
                }
                Pamd[i] = Pamd0[i]; //copy
            }
            Pamd[nPrim+nDual] = Pamd0[nPrim+nDual];
        }
    }
    opts.Pamd = Pamd; //link
    
    
    ADMMinfo res = {-1, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, NULL, NULL};
    //============================== EXECUTION =============================
    if(!mxIsSparse(prhs[0])){
        //Dense routine
        ADMMfloat *Px = (ADMMfloat*)mxGetDoubles(prhs[0]);
        ADMMfloat *Ax = (ADMMfloat*)mxGetDoubles(prhs[2]);
        eflag = superADMMsolverDense(Px, q, Ax, l, u, x, y, nPrim, nDual, opts, &res);

        
    } else {
        //Sparse routine

        //extract P
        ADMMfloat *Px = (ADMMfloat*)mxGetDoubles(prhs[0]);
        //We always cast -- even if the types are equal
        ADMMint *Pp = mwIndextoInt(mxGetJc(prhs[0]), (nPrim+1));
        ADMMint *Pi = mwIndextoInt(mxGetIr(prhs[0]), Pp[nPrim]);
        
        //and A
        ADMMfloat *Ax = (ADMMfloat*)mxGetDoubles(prhs[2]);
        ADMMint *Ap = mwIndextoInt(mxGetJc(prhs[2]), (nPrim+1));
        ADMMint *Ai = mwIndextoInt(mxGetIr(prhs[2]), Ap[nPrim]);
        
        eflag = superADMMsolverSparse(Px, Pp, Pi, q, Ax, Ap, Ai, l, u, x, y, nPrim, nDual, opts, &res);

        //Free because these are copies of the actual ones
        free(Pi);
        free(Pp);
        free(Ai);
        free(Ap);
    }
    
    if(nlhs == 1){
        free(y); //only free if it is not a return value
    }
    if(nlhs >= 3){
        plhs[2] = mxCreateDoubleScalar(eflag);
    }
    if(nlhs == 4){
        const char* resfn[] = {"nIter", "prim_res", "dual_res", "prim_tol", "dual_tol", "runtime", "objVal", "status", "Pamd"};
        mxArray *mxresptr = mxCreateStructMatrix(1,1,9, resfn);
        mxSetField(mxresptr, 0, "nIter", mxCreateDoubleScalar(res.nIter));
        mxSetField(mxresptr, 0, "prim_res", mxCreateDoubleScalar(res.prim_res));
        mxSetField(mxresptr, 0, "dual_res", mxCreateDoubleScalar(res.dual_res));
        mxSetField(mxresptr, 0, "prim_tol", mxCreateDoubleScalar(res.prim_tol));
        mxSetField(mxresptr, 0, "dual_tol", mxCreateDoubleScalar(res.dual_tol));
        mxSetField(mxresptr, 0, "runtime", mxCreateDoubleScalar(res.runtime));
        mxSetField(mxresptr, 0, "objVal", mxCreateDoubleScalar(res.objVal));
        mxSetField(mxresptr, 0, "status", mxCreateString(res.status));
        mxSetField(mxresptr, 0, "Pamd", mxPamd);
        plhs[3] = mxresptr;
    } else {
        free(Pamd);
    }
}


