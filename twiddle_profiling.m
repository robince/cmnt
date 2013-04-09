% Fast Complex Mersenne Number Transform
%
% Calculation of twiddle factors
% Profiling Test

function y = cmnt(x)
% Power of two length sequences only
d=length(x);
if( d ~= 2.^nextpow2(d) )
    error('Length of input vector must be a power of 2');
end

% Generate alpha's
% Primitive root of order 2^(p+1), p=61
alpha = complex( (uint64(2).^(uint64(2).^uint64(59))) , ui64(-3).^(uint64(2).^uint64(59)) );

%generate r - primtive root of order d
r = alpha .^ ((uint64(2).^uint64(62))./uint64(d));

% Original Method
for i=1:d
	twiddle1(i)=r.^uint64(i);
end

% Recursive Method
twiddle2(i)=r;
for i=2:d
	twiddle2(i)=r .* twiddle2(i-1);
end

% Original Method with Explicit Declaration
twiddle3=uint64(zeros(1,d));
for i=1:d
    twiddle3(i) = r.^uint64(i);
end

% Recursive Method with Explicit Declaration
twiddle4=uint64(zeros(1,d));
twiddle4(1)=r;
for i=2:d
	twiddle4(i) = r .* twiddle4(i-1);
end


