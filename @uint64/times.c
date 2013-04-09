/*
 * =============================================================
 Overloading the .* operator for uint64 elements of GF(q^2)

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


mod_t rmultmod(mod_t x, mod_t y)
{
	// uint68 * uint64 -> uint128 multiplication
	uint32 a[2],b[2];
	uint64 temp1, temp2, temp3, temp4, sum1, sum2;
	mod_t out2[2];
  
	a[0] = x;
	a[1] = (x >> 32);
	b[0] = y;
	b[1] = (y >> 32);

	// schoolboy, no knuth trick
	temp1 = ((uint64) a[0]) * b[0];
	
	temp2 = ((uint64) a[1]) * b[0];
	
	temp3 = ((uint64) a[0]) * b[1];
	
	temp4 = ((uint64) a[1]) * b[1];

	sum1 = (temp1 >> 32) + (temp2 & M32) + (temp3 & M32);
	sum2 = (sum1 >> 32) + (temp2 >> 32) + (temp3 >> 32);

	out2[0] = (temp1 & M32) + (sum1 << 32);
	out2[1] = sum2 + temp4;

	// reduce 128 bit result modulo M61
	temp1= (out2[0] & M61) + (out2[0] >> 61) + ((out2[1] << 3) & M61) + (out2[1] >> 58);
	temp1= (temp1 & M61) + (temp1 >> 61);
	temp1= (temp1 & M61) + (temp1 >> 61);

	return(temp1);
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  mod_t *xr, *xi, *yr, *yi, *zr, *zi;
  int ndimx, ndimy, i, n, ndim;
  int *dimsx, *dimsy, *dims;
  mxArray *temp;

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

  /* Special case if either of the inputs is a scalar */
  if( ((ndimy==2) && (dimsy[0]*dimsy[1]==1)) || ((ndimx==2) && (dimsx[0]*dimsx[1]==1)) ) /* either a scalar */
  {
	  if( (ndimy==2) && (dimsy[0]*dimsy[1]==1) )	/* y a scalar */
	  {
		  xr = (mod_t*) mxGetData(prhs[1]);	
		  xi = (mod_t*) mxGetImagData(prhs[1]);
		  yr = (mod_t*) mxGetData(prhs[0]);
		  yi = (mod_t*) mxGetImagData(prhs[0]);
		  ndim=ndimx;
		  dims=dimsx;
	  }
	  else											/* x a scalar */
	  {
		  xr = (mod_t*) mxGetData(prhs[0]);
		  xi = (mod_t*) mxGetImagData(prhs[0]);
		  yr = (mod_t*) mxGetData(prhs[1]);
		  yi = (mod_t*) mxGetImagData(prhs[1]);
		  ndim=ndimy;
		  dims=dimsy;
	  }
	  
	  if( (xi==NULL) && (yi==NULL) ) {
		  plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);
		  zr= (mod_t*) mxGetData(plhs[0]);
		  n = mxGetNumberOfElements(plhs[0]);
		  for(i=0;i<n;i++)
			  zr[i]=rmultmod(xr[0],yr[i]);
	  }
	  else if( (xi==NULL) ) {
		  plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
		  zr= (mod_t*) mxGetData(plhs[0]);
		  zi= (mod_t*) mxGetImagData(plhs[0]);
		  n = mxGetNumberOfElements(plhs[0]);
		  for(i=0;i<n;i++) {
			  zr[i]=rmultmod(xr[0],yr[i]);
			  zi[i]=rmultmod(xr[0],yi[i]);
		  }
	  }
	  else if( (yi==NULL) ) {
		  plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
		  zr= (mod_t*) mxGetData(plhs[0]);
		  zi= (mod_t*) mxGetImagData(plhs[0]);
		  n = mxGetNumberOfElements(plhs[0]);
		  for(i=0;i<n;i++)	{
			  zr[i]=rmultmod(xr[0],yr[i]);
			  zi[i]=rmultmod(xi[0],yr[i]);
		  }
	  }
	  else {
		  plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
		  zr= (mod_t*) mxGetData(plhs[0]);
		  n = mxGetNumberOfElements(plhs[0]);
		  zi= (mod_t*) mxGetImagData(plhs[0]);
		  for(i=0;i<n;i++) {
			  zr[i] = submod( rmultmod(xr[0],yr[i]),rmultmod(xi[0],yi[i]) );
			  zi[i] = addmod( rmultmod(xr[0],yi[i]),rmultmod(xi[0],yr[i]) );
		  }
	  }
  }

  else {			/* pointwise multiplication... vectors must be the same size */

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
		  zr[i]=rmultmod(xr[i],yr[i]);
  }
  else if( (xi==NULL) ) {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++) {
		  zr[i]=rmultmod(xr[i],yr[i]);
		  zi[i]=rmultmod(xr[i],yi[i]);
	  }
  }
  else if( (yi==NULL) ) {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  for(i=0;i<n;i++)	{
		  zr[i]=rmultmod(xr[i],yr[i]);
		  zi[i]=rmultmod(xi[i],yr[i]);
	  }
  }
  else {
	  plhs[0] = mxCreateNumericArray(ndimx, dimsx, mxUINT64_CLASS, mxCOMPLEX);
	  zr= (mod_t*) mxGetData(plhs[0]);
	  n = mxGetNumberOfElements(plhs[0]);
	  zi= (mod_t*) mxGetImagData(plhs[0]);
	  for(i=0;i<n;i++) {
		  zr[i] = submod( rmultmod(xr[i],yr[i]),rmultmod(xi[i],yi[i]) );
		  zi[i] = addmod( rmultmod(xr[i],yi[i]),rmultmod(xi[i],yr[i]) );
	  }
  }
  
  }
}
    
