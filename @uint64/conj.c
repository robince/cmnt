/*
 * =============================================================
 Overloading the conj() function for uint64 elements of GF(q^2)


 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, n;
    mod_t *pr, *pi;
    
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

	/* Create output array */
    plhs[0] = mxDuplicateArray(prhs[0]);

    pr = (mod_t *) mxGetData(plhs[0]);
    pi = (mod_t *) mxGetImagData(plhs[0]);
    n = mxGetNumberOfElements(plhs[0]);

    /* If there is an imaginary part of the input, negate it (bitwise NOT) */
    if(pi != NULL) { 
	for(i=0; i < n; i++) {
	    pi[i] = ((~ pi[i]) & M61);
	    }
	}
	
}
