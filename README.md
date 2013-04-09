cmnt
====

Complex Mersenne Number Transform - filtering in Finite Fields

This is a Matlab package from my 2005 MSc thesis on error free signal processing. 

This package overloads the Matlab `uint64` type with operations module M61 (Mersenne prime `2^61-1`).

The `uint64` type and overloaded operations allow implementation of number theoretic transforms (analogous to the FFT) but in the finite field defined modulo the Mersenne prime. Being in the finite field all operations are exact (no round off or truncation error). See the PDF of the thesis for more details. 

* `cmnt.m` : radix-2 Cooley-Tukey algoirthm.
* `cmnt2.m` : radix-2 Cooley-Tukey algoirthm - special case implementation for length 4 and 8.
* `cmnt3.m` : length N real transform as N/2 complex transform.
* `cmnt4.m` : two real signals of length N as a complex transform.
* `cmntc.m` : faster mex implementation 
* `icmnt.m` : inverse transform

* `cmntconv.m` : lossless convolution with the transform.
* `cmntfilt.m` : lossless filtering with the transform.

* `cmntrb.m` : Rader-Brenner algorithm (I think there wasn't any performance improvement)

* `fig_*.m` : some scripts to plot some demonstration figures
