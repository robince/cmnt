
function y = conv1(signal,filter)
%
% CONV1(x,h)
%	Simple Overlap Add Block Convolution Algorithm Using CMNT.
% 	Filter Length must be a factor of two, due to requirements of CMNT.
% 	Version 1: uses cmnt function

%rows/columns
Lx = length(signal); %signal length
Lh = length(filter); %filter length
Lt = Lx+Lh-1;        %output length

if( Lh~= 2.^nextpow2(Lh) )
    error('Length of impulse response must be a power of 2');
end
padding = zeros(1,Lh); %zero padding equal to filter length
h=ui64([filter padding]); %pad filter with zeros, same length as filter
H = cmnt(h); %Transform filter - only need do this once!!

num = floor(Lx/Lh); %number of complete times to filter
leftover = mod(Lx,Lh); % partial remainder left to filter

y = uint64(zeros(1,Lt)); %create empty output length of signal + length of filter-1

for n = 1:num
    x = ui64([signal(1+(n-1)*Lh:Lh*n) padding]); %extract section of signal and pad with zeros
    X = cmnt(x); %transform section of signal
    XH = X.*H; %multiply i.e. convolution in time domain = mult in freq domain
    xh = icmnt(XH); % inverse transform back to time domain
    y(1+(n-1)*Lh:((n-1)*Lh)+((2*Lh)-1)) = y(1+(n-1)*Lh:((n-1)*Lh)+((2*Lh)-1))+xh(1:Lh+Lh-1); % add new section into output, adding at the overlaps    
end

if leftover ~= 0
    start=1+(n*Lh);
    x = ui64([signal(start:Lx) zeros(1,(Lh-leftover)) padding]); %extract section of signal and pad with zeros
    X = cmnt(x); %transform section of signal
    XH = X.*H; %multiply i.e. convolution in time domain = mult in freq domain
    xh = icmnt(XH); % inverse transform back to time domain
    y(start:Lt) = y(start:Lt)+xh(1:(Lt-start)+1); % add new section into output, adding at the overlaps    
end

y = convert(y);
    
