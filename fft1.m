function y = fft1(x)
% Fast Complex Mersenne Number Transform
%
% Simple Radix-2 Transform

% Power of two length sequences only
d=length(x);
if( d ~= 2.^nextpow2(d) )
    error('Length of input vector must be a power of 2');
end

% Generate r - primtive root of order d
r = exp(-2*pi*sqrt(-1)/d);

% Table of twiddle factors
twiddle=zeros(1,d);
twiddle(1) = 1;
for i=2:d
    twiddle(i) = r .* twiddle(i-1);
end

% Output
y = cmntrecur(x, twiddle, d);

% -----------------------------------------------------
% cmntrecur - Recursive Subroutine to Implement radix-2
% -----------------------------------------------------

function y = cmntrecur(x, table, d)

n=length(x);
if (n==1)
    y=x;
else
    m = n/2;
    yT=cmntrecur(x(1:2:n),table, d);
    yB=cmntrecur(x(2:2:n),table, d);
    fiddle = d/n;
    z = table(1:fiddle:fiddle*m).*yB;
    y = [(yT+z) (yT-z)];
end
