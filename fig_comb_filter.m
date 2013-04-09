% fig_convolution1.m
% Generation of Plot
% Comparision of CMNT and FFT Error

% Load Signal
[in,fs,bits] = wavread('audio/test1_hi_norm_24.wav',[120000 644288]);
% Scale
x=2^(bits-1) .* in';
x=x(1:32768);
% Create Comb Filter
h=zeros(1,65536);
h(1)=128;
h(8000)=64;
h(16000)=32;
h(24000)=16;
% Zero Pad Signal
x=[x zeros(1,32768)];

% CMNT Calculations
X_64 = cmnt(ui64(x));
H_64 = cmnt(ui64(h));
Y_64 = X_64 .* H_64;
y_cmnt = convert(real(icmnt(Y_64)));

% FFT Calculations
X_fft = fft(x);
H_fft = fft(h);
Y_fft = X_fft .* H_fft;
y_fft = real(ifft(Y_fft));

figure(1)
subplot(3,1,1)
plot(h)
title('Impulse Response of Comb Filter');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 65536]);

subplot(3,1,2)
plot(x)
title('Audio Signal');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 65536]);

subplot(3,1,3)
plot(y_cmnt)
title('Filtered Signal');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 65536]);

% CMNT
Y_64 = cmnt(ui64(y_cmnt));
X_64_2 = Y_64 ./ H_64;
x_cmnt = convert(real(icmnt(X_64_2)));
% FFT
Y_fft = fft(y_fft);
X_fft = Y_fft ./ H_fft;
x_fft = real(ifft(X_fft));

% Plot Errors
figure(2)
subplot(2,1,1)
plot(abs(x_fft - x))
title('Absolute difference from FFT filtering');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 65536]);

subplot(2,1,2)
plot(abs(x_cmnt - x))
title('Absolute difference from CMNT filtering');
xlabel('Samples');
ylabel('Amplitude');
xlim([0 65536]);

