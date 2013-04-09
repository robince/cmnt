
mex getq.c

cd @double
mex -I../ ui64.c
cd ..

%%
cd @uint64

mex -I../ conj.c
mex -I../ convert.c
mex -I../ imag.c
mex -I../ minus.c
mex -I../ plus.c
mex -I../ power.c
mex -I../ r8.c
mex -I../ r8conj.c
mex -I../ sum.c
mex -I../ times.c
mex -I../ uminus.c