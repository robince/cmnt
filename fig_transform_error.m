% fig_transform_error.m
% Generation of Plot
% Comparision of CMNT and FFT Error

% Load Data
[in,fs,bits] = wavread('audio/test1_hi_norm_24.wav',[110001 634288]);
% Scale
x=2^(bits-1) .* in';
x64=ui64(x);
d=length(x64);
% FFT
t = cputime; 	x_fft = fft(x);			fft_time   = cputime - t
t = cputime; 	x_fft_out = ifft(x_fft);	ifft_time  = cputime - t
% CMNT
t = cputime;	x_cmnt = cmnt(x64);		cmnt_time  = cputime - t
t = cputime;	x_cmnt_out = icmnt(x_cmnt);	icmnt_time = cputime - t
% CMNT2
t = cputime;	temp1 = cmnt2(x64);		cmnt2_time  = cputime - t
t = cputime;	temp2 = (uint64(1)./uint64(d)).*conj(cmnt2(conj(temp1)));	icmnt2_time = cputime - t
% CMNT3
t = cputime;	temp1 = cmnt3(x64);		cmnt3_time  = cputime - t
% CMNT4
t = cputime;	temp1 = cmnt4(x64,x64);		cmnt4_time  = cputime - t
% CMNTRB
t = cputime;	temp1 = cmntrb(x64);		cmntrb_time  = cputime - t
t = cputime;	temp2 = (uint64(1)./uint64(d)).*conj(cmntrb(conj(temp1)));	icmntrb_time = cputime - t
% FFT1
t = cputime;	temp1 = fft1(x);		fft1_time  = cputime - t
t = cputime;	temp2 = (1./d).*conj(fft1(conj(temp1)));	ifft1_time = cputime - t

subplot(3,1,1);
plot( convert(x_cmnt_out) - x );
title('Error in reconstructed signal after CMNT');
xlabel('Samples');
ylabel('Error');
xlim([0 524288]);

subplot(3,1,2)
plot( real(x_fft_out) - x );
title('Error in reconstructed signal after FFT')
xlabel('Samples');
ylabel('Error');
xlim([0 524288]);

subplot(3,1,3)
plot( in );
title('Original Signal');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 524288]);
