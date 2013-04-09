% Fast Complex Mersenne Number Transform
%
% Rader-Brenner Algorithm
% always returns 1 x n, regardless of input

function y = cmntrb(x)

% Power of two length sequences only
d=length(x);
if( d ~= 2.^nextpow2(d) )
    error('Length of input vector must be a power of 2');
end

% Generate alpha - primitive root of order 2^(p+1), (p=61)
alpha = complex( (uint64(2).^(uint64(2).^uint64(59))),ui64(-3).^(uint64(2).^uint64(59)) );

% Generate r - primtive root of order d
r = alpha .^ ((uint64(2).^uint64(62))./uint64(d));

% Table of twiddle factors
twiddle=uint64(zeros(1,d));
twiddle(1) = uint64(1);
for k=2:d
    twiddle(k) = r .* twiddle(k-1);
end

twiddle = imag( twiddle ./ ((twiddle.*twiddle) - uint64(ones(1,d)) ));

% Output
y = cmntrecur(x, twiddle, d);

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
	fiddle=d/n;
	m = n/2;		% length of next level transform
	t = uint64(zeros(1,m));	% setup arrays (dynamic allocation slow)

	Q = sum(x(2:2:n));	% constant term
	
	% Setup first subsequence	
	t = [ x(1) x(3:2:n) ]; 
	
	% Transform first sub-sequences
	B = cmntrecur(t, table, d);
	
	% Setup + transform second subsequence
	t = [ (x(2) - x(n)) (x(4:2:n) - x(2:2:n-2)) ]; %c
	C = cmntrecur(t, table, d);

	t = table(1:fiddle:d/2);
	
	% Recreate Output
	D = complex(uminus(t .* imag(C)), t.*real(C));
	y = [(B-D) (B+D)]; 
	y(1) = B(1) + Q;
	y(m+1) = B(1) - Q;
end
