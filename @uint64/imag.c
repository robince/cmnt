/*
 * =============================================================
 Overloading the imag() function for uint64 elements of GF(q^2)


 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, n, ndim;
	int *dims;
    mod_t *pi, *z;
    
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

	ndim=mxGetNumberOfDimensions(prhs[0]);
	dims=mxGetDimensions(prhs[0]);
	pi=(mod_t*) mxGetImagData(prhs[0]);

    plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);

    z=(mod_t*) mxGetData(plhs[0]);
    n = mxGetNumberOfElements(plhs[0]);

    /* If there is an imaginary part of the input */
    if(pi != NULL) { 
	for(i=0; i < n; i++) {
	    z[i] = pi[i];
	    }
	}
	else {
	for(i=0; i < n; i++) {
	    z[i] = 0;
	    }
	}
	
}
