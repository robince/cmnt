/*
 * =============================================================
 Overloading the .^ operator for uint64 elements of GF(q^2)
 (binary exponentiation)

 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

bool isOdd(mod_t x)
{
  return (x % 2) == 1;
}


mod_t addmod(mod_t x, mod_t y)
{
  // addition
  mod_t temp = x + y;
  // modulo reduction
  return( (temp & M61) + (temp >> 61) );
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
  mod_t *xr, *xi, *yr, *yi, *r;
  mod_t k, br, bi, tr, ti, cr, ci;
  int i, n, ndim;
  int *dims;
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
  if ( mxIsComplex(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1]) != 1) {
    mexErrMsgTxt("Input y must be a real scalar.");
  }

  ndim= mxGetNumberOfDimensions(prhs[0]);
  dims= mxGetDimensions(prhs[0]);
  xr= (mod_t*) mxGetData(prhs[0]);
  xi= (mod_t*) mxGetImagData(prhs[0]);
  n=    mxGetNumberOfElements(prhs[0]);
  r=    mxGetData(prhs[1]);

  if( xi==NULL) {				/* input is real */
	  plhs[0]=mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxREAL);
	  yr= (mod_t*) mxGetData(plhs[0]);

	  for(i=0;i<n;i++) {
		  k=r[0];
		  cr=xr[i];
		  br=1;
		  while(k>=1)
		  {
			  if( isOdd(k) )
			  {
				  br = rmultmod(cr,br);
				  k = k - 1;
			  }
			  cr = rmultmod(cr,cr);
			  k = k / 2;
		  }
		  yr[i] = br;
	  }
  }
  else {						/* input is complex */
	  plhs[0]=mxCreateNumericArray(ndim, dims, mxUINT64_CLASS, mxCOMPLEX);
	  yr= (mod_t*) mxGetData(plhs[0]);
	  yi= (mod_t*) mxGetImagData(plhs[0]);

	  for(i=0;i<n;i++) {
		  k=r[0];
		  cr=xr[i];
		  ci=xi[i];
		  br=1;
		  bi=0;

		  while(k>=1)
		  {
			  if( isOdd(k) )
			  {		/* b = c * b */
				  tr = submod( rmultmod(cr,br),rmultmod(ci,bi) );
				  ti = addmod( rmultmod(cr,bi),rmultmod(ci,br) );
				  br=tr;
				  bi=ti;
				  k = k - 1;
			  }
					/* c = c^2 */
			  tr = submod( rmultmod(cr,cr),rmultmod(ci,ci) );
			  ti = addmod( rmultmod(cr,ci),rmultmod(ci,cr) );
			  cr=tr;
			  ci=ti;
			  k = k / 2;
		  }
		  yr[i] = br;
		  yi[i] = bi;
	  }
  }
}




