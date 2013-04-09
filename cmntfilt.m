function y = cmntfilt(b_in,x_in,ncmnt)

error(nargchk(2,3,nargin));

% Check inputs are vectors, and ensure column
[mx,n] = size(x_in);
if mx == 1
    x_in = x_in(:);    % turn row into a column
elseif n ~=1
    error('Function only accepts vector inputs')
end
nx = size(x_in,1);

[mb,n] = size(b_in);
if mb == 1
    b_in = b_in(:);    % turn row into a column
elseif n ~=1
    error('Function only accepts vector inputs')
end
nb = size(b_in,1);

if nargin < 3
% figure out which ncmnt and L to use - same as fftfilt
    if nb >= nx     % take a single FFT in this case
        ncmnt = 2^nextpow2(nb+nx-1);
        L = nx;
    else
        fftflops= [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
       149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];
        n = 2.^(1:20);
        validset = find(n>(nb-1));   % must have ncmnt > (nb-1)
        n = n(validset); 
        fftflops = fftflops(validset);
        % minimize (number of blocks) * (number of flops per cmnt)
        L = n - (nb - 1);
        [dum,ind] = min( ceil(nx./L) .* fftflops );
        ncmnt = n(ind);
        L = L(ind);
    end

else  % ncmnt is given
    if ncmnt < nb
        ncmnt = nb;
    end
    ncmnt = 2.^(ceil(log(ncmnt)/log(2))); % force this to a power of 2 for speed
    L = ncmnt - nb + 1;
end

% Convert to Finite Field
x = ui64(x_in);
b = ui64([b_in; zeros(ncmnt-nb,1)]); %convert and pad

% Precompute twiddle factors for transform
% Generate alpha - primitive root of order 2^(p+1), (p=61)
alpha = complex( (uint64(2).^(uint64(2).^uint64(59))),ui64(-3).^(uint64(2).^uint64(59)) );
% Generate r - primtive root of order d
r = alpha .^ ((uint64(2).^uint64(62))./uint64(ncmnt));
% Table of twiddle factors
twiddle=uint64(zeros(ncmnt,1));
twiddle(1) = 1;
for i=2:ncmnt
    twiddle(i) = r .* twiddle(i-1);
end

B = cmntrecur(b, twiddle, ncmnt);

y = uint64(zeros(size(x)));

istart = 1;
while istart <= nx
    iend = min(istart+L-1,nx);
    if (iend - istart) == 0
        X = x(istart(ones(ncmnt,1)),:);  % need to fft a scalar
    else
%	block = [x(istart:iend,:); zeros(ncmnt-(iend-istart),1)];
        X = cmntrecur([x(istart:iend,:); zeros(ncmnt-(iend-istart+1),1)], twiddle, ncmnt);
    end
    Y=(uint64(1)./uint64(ncmnt)).*conj(cmntrecur(conj(X.*B), twiddle, ncmnt));
    yend = min(nx,istart+ncmnt-1);
    y(istart:yend,:) = y(istart:yend,:) + Y(1:(yend-istart+1),:);
    istart = istart + L;
end

if ~any(imag(b_in)) & ~any(imag(x_in))
	y = real(y);
end

y=convert(y);

if (mx == 1)&(size(y,2) == 1)
    y = y(:).';    % turn column back into a row
end


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
    
    y = [ (A+B); (C + N - E + M); (K + L); (C - P + E - O); ...
          (A-B); (C - N - E - M); (K - L); (C + P + E + O)];
else
    m = n/2;
    yT=cmntrecur(x(1:2:n),table, d);
    yB=cmntrecur(x(2:2:n),table, d);
    fiddle = d/n;
    z = table(1:fiddle:fiddle*m).*yB;
    y = [(yT+z); (yT-z)];
end
