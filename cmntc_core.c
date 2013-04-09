/*
 * =============================================================
 Fast-CMNT following Numerical Recipes Algorithm

 Robin Ince
 * =============================================================
 */


#include "i64.h"
#include "mex.h"

#define		SWAP(a,b)			tempr=(a);(a)=(b);(b)=tempr

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
	uint64 temp1, temp2, temp3, sum1;
	mod_t out2[2];
  
	a[0] = x;
	a[1] = (x >> 32);
	b[0] = y;
	b[1] = (y >> 32);

	// schoolboy, no knuth trick
	temp1 = ((uint64) a[0]) * b[0];
	
	temp2 = ((uint64) a[1]) * b[0];
	
	temp3 = ((uint64) a[0]) * b[1];
	
	sum1 = (temp1 >> 32) + (temp2 & M32) + (temp3 & M32);

	out2[0] = (temp1 & M32) + (sum1 << 32);
	out2[1] = ((sum1 >> 32) + (temp2 >> 32) + (temp3 >> 32)) + (((uint64) a[1]) * b[1]);

	/* reduce 128 bit result modulo M61*/
	temp1= (out2[0] & M61) + (out2[0] >> 61) + ((out2[1] << 3) & M61) + (out2[1] >> 58);
	temp1= (temp1 & M61) + (temp1 >> 61);
	temp1= (temp1 & M61) + (temp1 >> 61);

	return(temp1);
}

/* The gateway routine */
void 
mexFunction( int nlhs, mxArray *plhs[],  int nrhs, const mxArray *prhs[])
{
    int i, j, k, n, m, istep, tind;
    mod_t *xr, *xi, *rr, *ri, *tr, *ti, tempr, tempi, A, B, C;
	mxArray *twiddle;
    
    /* Check for proper number of input and output arguments */    
    if (nrhs != 2) {
		mexErrMsgTxt("Two input arguments required.");
    } 
    if(nlhs > 1){
		mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Check data type of input arguments  */
    if (!(mxIsUint64(prhs[0])) || !(mxIsUint64(prhs[1]))){
		mexErrMsgTxt("Input arguments must be of type uint64.");
    }	
	if (!(mxIsComplex(prhs[0])) || !(mxIsComplex(prhs[1]))){
		mexErrMsgTxt("Input arguments must be complex");
    }

	/* Create output arrays */
	plhs[0] = mxDuplicateArray(prhs[0]);
	twiddle = mxDuplicateArray(prhs[0]);

	xr = (mod_t *) mxGetData(plhs[0]);
	xi = (mod_t *) mxGetImagData(plhs[0]);

	rr = (mod_t *) mxGetData(prhs[1]);
	ri = (mod_t *) mxGetImagData(prhs[1]);

	tr = (mod_t *) mxGetData(twiddle);
	ti = (mod_t *) mxGetImagData(twiddle);

	n = mxGetNumberOfElements(plhs[0]);
	k=1;

	/* Table of Twiddle Factors */
	tr[0]=1;
	ti[0]=0;
	for(i=1; i<(n/2); i++)
	{
		j=i-1;
		/* complex multiplication as three real additions */	
		A = rmultmod(rr[0],tr[j]);
		B = rmultmod(ri[0],ti[j]);
		C = rmultmod(rr[0]+ri[0],tr[j]+ti[j]);
		
		tr[i] = submod(A,B);
		ti[i] = submod(C,addmod(A,B));
	}
	
	j=0;
	/* Bit Reversal Data Reordering */
	for(i=0; i<n; i++)
	{
		if (i <= j)
		{
			SWAP(xr[j],xr[i]);
			SWAP(xi[j],xi[i]);
		}
		m = n/2;
		while (j > m-1)
		{
			j = j - m;
			m = m/2;
			if (m < 1)
				break;
		}
		j = j + m;
	}

	/* Transform */	
	do
	{
		istep = 2*k;
		for (m=0; m<k; m++)
		{
			/* calculates twiddle(nm.2k) */
			tind = m*(n/istep);
			/* need r^mn/istep */
			for (i=m; i<n; i+=istep)
			{
				j=i+k;
				/* complex multiplication as three real multiplications */	
				A = rmultmod(tr[tind],xr[j]);
				B = rmultmod(ti[tind],xi[j]);
				C = rmultmod(tr[tind]+ti[tind],xr[j]+xi[j]);
				
				tempr = submod(A,B);
				tempi = submod(C,addmod(A,B));
				
				xr[j] = submod(xr[i], tempr);
				xi[j] = submod(xi[i], tempi);
				xr[i] = addmod(xr[i], tempr);
				xi[i] = addmod(xi[i], tempi);
			}
		}
		k = istep;
	}
	while (k < n);

	/* clean up */
	mxDestroyArray(twiddle);
}

