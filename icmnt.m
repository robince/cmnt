% Inverse Complex Mersenne Number Transform
function y=icmnt(x)

d=length(x);
y=(uint64(1)./uint64(d)).*conj(cmnt(conj(x)));
