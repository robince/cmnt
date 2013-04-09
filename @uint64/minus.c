/*
 * =============================================================
 Overloading the - operator for uint64 elements of GF(q^2)

 Robin Ince
 * =============================================================
 */



#include "i64.h"
#include "mex.h"

mod_t submod(mod_t x, mod_t y)
{
  // add x + -y
  mod_t temp = x + ((~ y) & M61);
  // modulo reduction
  return( (temp & M61) + (temp >> 61));
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mod_t *xr, *xi, *yr, *yi, *zr, *zi;
  int ndimx, ndimy, i, n;
  int *dimsx, *dimsy;

  /*  Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using
     mexErrMsgTxt within an if statement. It will never
     get to the else statement if mexErrMsgTxt is executed.
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */
  if (nrhs != 2) 
    mexErrMsgTxt("Two inputs required.");
  if (nlhs > 1) 
    mexErrMsgTxt("Too many output arguments");

  if (!mxIsUint64(prhs[0]) || !mxIsUint64(prhs[1])) {
    mexErrMsgTxt("Inputs must be of type uint64");
  }
  
  /* Check to make sure the input arguments have same dimensions */
  ndimx=mxGetNumberOfDimensions(prhs[0]);
  ndimy=mxGetNumberOfDimensions(prhs[1]);
  dimsx=mxGetDimensions(prhs[0]);
  dimsy=mxGetDimensions(prhs[1]);

  if ( ndimx != ndimy )
	 mexErrMsgTxt("Number of array dimensions must match");

  /* */
  for (i=0; i<ndimx; i++) {
	  if (dimsx[i] != dimsy[i])
		  mexErrMsgTxt("Matrix dimensions must agree.");
  }
	  
  xr = (mod_t*) mxGetData(prhs[0]);
  xi = (mod_t*) mxGetImagData(prhs[0]);
  yr = (mod_t*) mxGetData(prhs[1]);
  yi = (mod_t*) mxGetImagData(prhs[1]);
	
  if( (xi==NULL) && (yi==NULL) ) {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxREAL);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++)
		  zr[i]=submod(xr[i],yr[i]);
  }
  else if( (xi==NULL) ) {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++) {
		  zr[i]=submod(xr[i],yr[i]);
		  zi[i]=((~ yi[i]) & M61);
	  }
  }
  else if( (yi==NULL) ) {
	  plhs[0] = mxDuplicateArray(prhs[0]);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++)
		  zr[i]=submod(xr[i],yr[i]);
  }
  else {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  for(i=0;i<n;i++) {
		  zr[i]=submod(xr[i],yr[i]);
		  zi[i]=submod(xi[i],yi[i]);
	  }
  }
  
}
    
