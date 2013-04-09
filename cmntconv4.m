function y = conv4(signal,filter)
%
% CONV4(x,h)
%	Simple Overlap Add Block Convolution Algorithm Using CMNT.
% 	Filter Length must be a factor of two, due to requirements of CMNT.
%	Version 3: 	implements transform internally to reuse table of twiddle factors
%			includes n=4 optimisation
%			computes two real blocks as one complex transform


Lx = length(signal); %signal length
Lh = length(filter); %filter length
Lt = Lx+Lh-1;        %output length

if( Lh ~= 2.^nextpow2(Lh) )
    error('Length of impulse response must be a power of 2');
end

padding = zeros(1,Lh); %zero padding equal to filter length
h=ui64([filter padding]); %pad filter with zeros, same length as filter
d=Lh*2; %length of transform

% Precompute Twiddle Factors
% Generate alpha - primitive root of order 2^(p+1), (p=61)
alpha = complex( (uint64(2).^(uint64(2).^uint64(59))),ui64(-3).^(uint64(2).^uint64(59)) );

% Generate r - primtive root of order d
r = alpha .^ ((uint64(2).^uint64(62))./uint64(d));

% Table of twiddle factors
twiddle=uint64(zeros(1,d));
twiddle(1) = 1;
for i=2:d
    twiddle(i) = r .* twiddle(i-1);
end

H = cmntrecur(h,twiddle,d); %Transform filter - only need do this once!!

num = floor(Lx/Lh); % number of complete times to filter
num2 = floor(num/2); % number of pairs of blocks
leftover = mod(Lx,Lh); % partial remainder left to filter

y = uint64(zeros(1,Lt)); %create empty output length of signal + length of filter

for n = 1:2:num2*2 	% calculate for pairs of blocks
    x1 = ui64([signal(1+(n-1)*Lh:Lh*n) padding]); % extract 2 blocks
    x2 = ui64([signal(1+(n)*Lh:Lh*(n+1)) padding]);
    X = cmntrecur(complex(x1,x2),twiddle,d); % transform as complex signal
    XH = X.*H; % convolution
    xh = (uint64(1)./uint64(d)).*conj(cmntrecur(conj(XH),twiddle,d));
    %xh1 = real((uint64(1)./uint64(d)).*conj(cmntrecur(conj(XH),twiddle,d))); % result block 1
    %xh2 = imag((uint64(1)./uint64(d)).*conj(cmntrecur(conj(XH),twiddle,d))); % result block 2
    y(1+(n-1)*Lh:((n-1)*Lh)+(d-1)) = y(1+(n-1)*Lh:((n-1)*Lh)+(d-1))+real(xh(1:d-1)); % add block 1
    y(1+n*Lh:(n*Lh)+(d-1)) = y(1+n*Lh:(n*Lh)+(d-1))+imag(xh(1:d-1)); % add block 2
end

if (num/2) ~= num2 	% one more complete block
    n=num;
    x = ui64([signal(1+(n-1)*Lh:Lh*n) padding]); %extract section of signal and pad with zeros
    X = cmntrecur(x,twiddle,d); %transform section of signal
    XH = X.*H; %multiply i.e. convolution in time domain = mult in freq domain
    xh = (uint64(1)./uint64(d)).*conj(cmntrecur(conj(XH),twiddle,d)); % inverse transform back to time domain
    y(1+(n-1)*Lh:((n-1)*Lh)+(d-1)) = y(1+(n-1)*Lh:((n-1)*Lh)+(d-1))+xh(1:d-1); % add new section into output, adding at the overlaps    
end

if leftover ~= 0
    n=num;
    start=1+(n*Lh);
    x = ui64([signal(start:Lx) zeros(1,(Lh-leftover)) padding]); %extract section of signal and pad with zeros
    X = cmntrecur(x,twiddle,d); %transform section of signal
    XH = X.*H; %multiply i.e. convolution in time domain = mult in freq domain
    xh = (uint64(1)./uint64(d)).*conj(cmntrecur(conj(XH),twiddle,d)); % inverse transform back to time domain
    y(start:Lt) = y(start:Lt)+xh(1:(Lt-start)+1); % add new section into output, adding at the overlaps    
end

y = convert(y);
 
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


