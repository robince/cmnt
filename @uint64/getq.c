/*
 * =============================================================
 Return modulus - Q = 0x1FFFFFFFFFFFFFFF


 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mod_t *z;
  int ndim=2;
  int dims[2];
  dims[0]=1;
  dims[1]=1;

  if (nlhs > 1) 
	mexErrMsgTxt("Too many outputs");
 
  
  /* Set the output pointer to the output matrix. */
  plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);
  
  /* Create a C pointer to a copy of the output matrix. */
  z = (mod_t *) (mxGetData(plhs[0]));

  // output q
  *z = M61;
}
