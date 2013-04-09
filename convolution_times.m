% fig_convolution1.m
% Generation of Plot
% Comparision of CMNT and FFT Error

% Load Signal
[in,fs,bits] = wavread('audio/test1_hi_norm_24.wav',[120001 644288]);
% Scale
x=2^(bits-1) .* in';
x=x(1:32768);

% Load Impulse Response
[in,fs,bits] = wavread('audio/flutter_echo_ir_stereo_24.wav');
h=2^(bits-1) .* in(1:4096,1)';

% fft convolution
t = cputime; 	temp = conv(x,h);			conv_time = cputime - t
% conv1
t = cputime; 	temp = cmntconv1(x,h);			conv1_time = cputime - t
% conv2
t = cputime; 	temp = cmntconv2(x,h);			conv2_time = cputime - t
% conv3
t = cputime; 	temp = cmntconv3(x,h);			conv3_time = cputime - t
% conv4
t = cputime; 	temp = cmntconv4(x,h);			conv4_time = cputime - t
% fftfilt
t = cputime; 	temp = fftfilt(h,x);			fftfilt_time = cputime - t
% cmntfilt
t = cputime; 	temp = cmntfilt(h,x);			cmntfilt_time = cputime - t
