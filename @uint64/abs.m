function y=abs(x)

y = ( (real(x).^uint64(2)) + (imag(x).^uint64(2)) ).^(uint64(1)./uint64(2))
