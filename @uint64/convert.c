/*
 * =============================================================
 Convert usigned 64 bit representation of GF(q^2) element to signed representation
 (double)

 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

double gf2double(mod_t x) {

	if( (x & SIGN_BIT) !=0 )
		return( (double) (((int64) x) - M61 ));
	else
		return( (double) (int64) x);
}

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, n, ndim, *dims;
    mod_t *pr, *pi;
	double *zr, *zi;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 1) {
	mexErrMsgTxt("One input argument required.");
    } 
    if(nlhs > 1){
	mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsUint64(prhs[0]))){
	mexErrMsgTxt("Input argument must be of type uint64.");
    }	

    pr = (mod_t *) mxGetData(prhs[0]);
    pi = (mod_t *) mxGetImagData(prhs[0]);
    n = mxGetNumberOfElements(prhs[0]);
    ndim=mxGetNumberOfDimensions(prhs[0]);
    dims=mxGetDimensions(prhs[0]);
	
	if(pi==NULL) {
	    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	}
	else {
		plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
	}
	zr=mxGetPr(plhs[0]);
	zi=mxGetPi(plhs[0]);

    /* Perform task... */
    for(i=0; i < n; i++) {
		zr[i]=gf2double(pr[i]);
    }

	if(pi!=NULL) {
    for(i=0; i < n; i++) {
		zi[i]=gf2double(pi[i]);
    }
	}
}
