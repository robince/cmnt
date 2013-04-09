/*
 * =============================================================
 Special function for multiplication by r_8 (2^30)
 
 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"



mod_t addmod(mod_t x, mod_t y)
{
  // addition
  mod_t temp = x + y;
  // modulo reduction
  return( (temp & M61) + (temp >> 61));
}

mod_t submod(mod_t x, mod_t y)
{
  // add x + -y
  mod_t temp = x + ((~ y) & M61);
  // modulo reduction
  return( (temp & M61) + (temp >> 61));
}

mod_t timesr(mod_t x)
{
  /* multiply by 2^30 and reduce */
  mod_t temp = ((x << 30) & M61);
  temp += (x >> 31);
  /*return( (temp & M61) + (temp >> 61));*/
  return(temp);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mod_t *xr, *xi, *yr, *yi, *zr, *zi;
  mod_t temp1, temp2;
  int ndimx, ndimy, i, n, ndim;
  int *dimsx, *dimsy, *dims;
  mxArray *temp;

  /*  Check for proper number of arguments. */
  /* NOTE: You do not need an else statement when using
     mexErrMsgTxt within an if statement. It will never
     get to the else statement if mexErrMsgTxt is executed.
     (mexErrMsgTxt breaks you out of the MEX-file.) 
  */
  if (nrhs != 1) 
    mexErrMsgTxt("One input required.");
  if (nlhs > 1) 
    mexErrMsgTxt("Too many output arguments");

  if (!mxIsUint64(prhs[0])) {
    mexErrMsgTxt("Inputs must be of type uint64");
  }
  	  
  xr = (mod_t*) mxGetData(prhs[0]);
  xi = (mod_t*) mxGetImagData(prhs[0]);

  ndimx=mxGetNumberOfDimensions(prhs[0]);
  dimsx=mxGetDimensions(prhs[0]);

  /* create output array */
  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
  
  if( (xi==NULL) ) {	/* real input */
	  zr= (mod_t*) mxGetData(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++) {
		  zr[i] = timesr(xr[i]);
		  zi[i] = ((~zr[i]) & M61);
	  }
  }
  else { 		/* complex input */
	  zr= (mod_t*) mxGetData(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++) {
		  temp1 = timesr(xr[i]);
		  temp2 = timesr(xi[i]);
		  zr[i] = addmod( temp1,temp2 );
		  zi[i] = submod( temp2,temp1 );
	  }
  }
  
}
    
