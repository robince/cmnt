/*
 * =============================================================
 Sum - just sum's all the elements in the argument (should be a vector)
  
 (61 bit ones' compliment)

 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, n, ndimx, dims[2], k;
    int *dimsx;
    mod_t *xr, *xi, *zr, *zi;
    int ndim=2;
    dims[0]=1;
    dims[1]=1;
    
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

    xr = (mod_t *) mxGetData(prhs[0]);
    xi = (mod_t *) mxGetImagData(prhs[0]);
    n = mxGetNumberOfElements(prhs[0]);

        /* Real Input... */
    if(xi == NULL) {
	plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);
	zr=(mod_t*) mxGetData(plhs[0]);	
	zr[0] = 0;
	for(i=0; i < n; i++) {
	    zr[0] += xr[i];
	    k++;
	    /* only need to reduce every 3 additions */
	    if(k > 2) {
		    zr[0] = (zr[0] & M61) + (zr[0] >> 61);
		    k=0;
	    }
	    }
	if(k != 0)	/* reduction outstanding */
	    zr[0] = (zr[0] & M61) + (zr[0] >> 61);
	}

    else {	/* Complex Input... */
	plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
	zr=(mod_t*) mxGetData(plhs[0]);
	zi=(mod_t*) mxGetImagData(plhs[0]);
	zr[0] = 0;
	zi[0] = 0;
	for(i=0; i < n; i++) {
	    zr[0] += xr[i];
	    zi[0] += xi[i];
	    k++;
	    /* only need to reduce every 3 additions */
	    if(k > 2) {
		    zr[0] = (zr[0] & M61) + (zr[0] >> 61);
		    zi[0] = (zi[0] & M61) + (zi[0] >> 61);
		    k=0;
	    }
	    }
	if(k != 0)	/* reduction outstanding */
	    zr[0] = (zr[0] & M61) + (zr[0] >> 61);
	    zi[0] = (zi[0] & M61) + (zi[0] >> 61);
	}	
}
