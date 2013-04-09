function z=rdivide(x,y)
% Division in GF(M61^2) (multiplication by inverse)
% Real only
% Only needed to renormalise inverse transform & calculate primitve roots
q=getq;

if	( isreal(x) & isreal(y) )	% Both Real
	z = x .* (y.^(q-uint64(2)));
elseif	( isreal(y) )			% X Complex, Y Real
	z = x .* (y.^(q-uint64(2)));
elseif	( isreal(x) ) 			% X Real, Y Complex
	z = complex( x.*real(y),uminus(x.*imag(y)) );
	z = z .* ((real(y).*real(y) + imag(y).*imag(y)) ).^(q-uint64(2));
else					% Both Complex
	z = complex( (real(x).*real(y) + imag(x).*imag(y)), (imag(x).*real(y) - real(x).*imag(y)));
	z = times(z, (real(y).*real(y) + imag(y).*imag(y)).^(q-uint64(2)));
end

