function y = cmnt(x,y)
% Fast Complex Mersenne Number Transform
%
% Simple Radix-2 Transform
% including n=4 optimization
% computes two real transforms as complex function

% Power of two length sequences only
nx=length(x);
ny=length(y);

if( nx ~= ny)
    error('Vectors must be same length');
end
n=nx;
    
if( n ~= 2.^nextpow2(n) )
    error('Length of input vectors must be a power of 2');
end

% Generate alpha - primitive root of order 2^(p+1), (p=61)
alpha = complex( (uint64(2).^(uint64(2).^uint64(59))),ui64(-3).^(uint64(2).^uint64(59)) );

% Generate r - primtive root of order n
r = alpha .^ ((uint64(2).^uint64(62))./uint64(n));

% Split first and second half into complex and real part
x_comp = complex( x , y );

% Table of twiddle factors
twiddle=uint64(zeros(1,n));
twiddle(1) = 1;
for i=2:n
    twiddle(i) = r .* twiddle(i-1);
end

% Transform
y_comp = cmntrecur(x_comp, twiddle, n);

% Reconstruct Real Output
% Declare arrays (most important; dynamic allocation kills performance)
X = uint64(zeros(1,n));
Y = X;

i2 = uint64(1) ./ uint64(2); % Precompute 1/2 (multiplication cheaper than division)
% Want real and imaginary part of y_comp(1) here
% done this way to get multiplication by i2 out of the loop
% Quicker to do vector .* scalar than to have scalar mult in the loop.
X(1) = y_comp(1) + conj(y_comp(1));
Y(1) = conj(y_comp(1)) - y_comp(1); 

for m = 2:n
	a=y_comp(m);
	b=conj(y_comp((n+2)-m));
	X(m) = (a + b);
	Y(m) = (b - a);
end

y = [(X.*i2) (Y.*complex(uint64(0),i2))];

% -----------------------------------------------------
% cmntrecur - Recursive Subroutine to Implement radix-2
% -----------------------------------------------------

function y = cmntrecur(x, table, d)

n=length(x);
if (n==1)
    y=x;
elseif (n==8)
    C = x(1) - x(5);
    D = x(2) - x(6);
    Er = x(3) - x(7);
    E = complex(-imag(Er),real(Er));
    F = x(8) - x(4);
    G = x(1) + x(5);
    H = x(3) + x(7);
    I = x(4) + x(8);
    J = x(2) + x(6);
    A = G + H;
    B = I + J;
    K = G - H;
    Lr = I - J;
    L = complex(-imag(Lr),real(Lr));
    M = r8conj(F);
    N = r8(D);
    O = r8conj(D);
    P = r8(F);
    
    y = [ (A+B) (C + N - E + M) (K + L) (C - P + E - O) ...
          (A-B) (C - N - E - M) (K - L) (C + P + E + O)];
else
    m = n/2;
    yT=cmntrecur(x(1:2:n),table, d);
    yB=cmntrecur(x(2:2:n),table, d);
    fiddle = d/n;
    z = table(1:fiddle:fiddle*m).*yB;
    y = [(yT+z) (yT-z)];
end
