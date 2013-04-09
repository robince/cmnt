/*
 * =============================================================
 Convert usigned 64 bit representation of GF(q^2) element to signed representation
 (double)

 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

mod_t double2gf(double x) {

	if( x < 0 )
		return( (mod_t) (((int64) x) + M61 ));
	else
		return( (mod_t) x);
}

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, n, ndim, *dims;
    double *pr, *pi;
    mod_t *zr, *zi;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 1) {
	mexErrMsgTxt("One input argument required.");
    } 
    if(nlhs > 1){
	mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input argument  */
    if (!(mxIsDouble(prhs[0]))){
	mexErrMsgTxt("Input argument must be of type double.");
    }	

    pr = mxGetData(prhs[0]);
    pi = mxGetImagData(prhs[0]);
    n = mxGetNumberOfElements(prhs[0]);
    ndim=mxGetNumberOfDimensions(prhs[0]);
    dims=mxGetDimensions(prhs[0]);
	
	if(pi==NULL) {
	    plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);
	}
	else {
		plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
	}
	zr=(mod_t *) mxGetPr(plhs[0]);
	zi=(mod_t *) mxGetPi(plhs[0]);

    /* Perform task... */
    for(i=0; i < n; i++) {
		zr[i]=double2gf(pr[i]);
    }

	if(pi!=NULL) {
    for(i=0; i < n; i++) {
		zi[i]=double2gf(pi[i]);
    }
	}
}
