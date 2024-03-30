% Clear command window and close all figures
clc
close all

% Define parameters
tfinal = 0.01; % Final time for the analog signal
t = 0:0.00001:tfinal; % Time vector with sampling interval

% Create a new figure for plotting
figure;

% Generate analog signal composed of two cosine waves
xanalog = cos(2*pi*400*t) + cos(2*pi*700*t);

% Plot analog signal in time domain
subplot(4,2,2);
plot(t, xanalog , 'r-');
xlabel("Time");
ylabel("Amplitude");
title("Analog signal");

% Compute Fourier transform of the analog signal
x_DFT = abs(fftshift(fft(xanalog)));
l = length(x_DFT);
x_length = -l/2:l/2-1;

% Plot analog signal in frequency domain
subplot(4,2,1);
stem(100*x_length+50, x_DFT);
axis([-800 800 0 600]);
xlabel("Frequency");
ylabel("Magnitude");
title("Analog signal in Frequency Domain");

% Critical Sampling (fs=2fm)
fs = 1400;
tsamp = 0:1/fs:14/fs; % Time vector for critical sampling
xsampled = cos(2*pi*400*tsamp) + cos(2*pi*700*tsamp); % Sampled signal
xsampled_DFT = abs(fftshift(fft(xsampled))); % Compute DFT of the sampled signal
xsampled_length = -length(xsampled_DFT)/2:length(xsampled_DFT)/2-1;

% Plot sampled signal in frequency domain
subplot(4,2,3);
stem(100*xsampled_length+50, xsampled_DFT);
axis([-800 800 0 15])
xlabel("Frequency");
ylabel("Magnitude");
title("Critical Sampling(fs=2fm)");

% Reconstruct the sampled signal
xrecon = ifft(abs(xsampled_DFT));
subplot(4,2,4);
plot(tsamp, xrecon);
xlabel("Time");
ylabel("Amplitude");
title("Critical Sampling(fs=2fm)");

% Under Sampling (fs<2fm)
fs = 700;
tsamp = 0:1/fs:6/fs; % Time vector for under sampling
xsampled = cos(2*pi*400*tsamp) + cos(2*pi*700*tsamp); % Sampled signal
xsampled_DFT = abs(fftshift(fft(xsampled))); % Compute DFT of the sampled signal
xsampled_length = -length(xsampled_DFT)/2:length(xsampled_DFT)/2-1;

% Plot sampled signal in frequency domain
subplot(4,2,5);
stem(100*xsampled_length+50, xsampled_DFT);
axis([-800 800 0 15])
xlabel("Frequency");
ylabel("Magnitude");
title("Under Sampling(fs<2fm)");

% Reconstruct the sampled signal
xrecon = ifft(abs(xsampled_DFT));
subplot(4,2,6);
plot(tsamp, xrecon);
xlabel("Time");
ylabel("Amplitude");
title("Under Sampling(fs<2fm)");

% Over Sampling (fs>2fm)
fs = 5000;
tsamp = 0:1/fs:0.01; % Time vector for over sampling
xsampled = cos(2*pi*400*tsamp) + cos(2*pi*700*tsamp); % Sampled signal
xsampled_DFT = abs(fft(xsampled)); % Compute DFT of the sampled signal
xsampled_length = -length(xsampled_DFT)/2:length(xsampled_DFT)/2-1;

% Plot sampled signal in frequency domain
subplot(4,2,7);
stem(100*xsampled_length, xsampled_DFT);
axis([-800 800 0 15])
xlabel("Frequency");
ylabel("Magnitude");
title("Over Sampling(fs>2fm)");

% Reconstruct the sampled signal
xrecon = ifft(abs(xsampled_DFT));
subplot(4,2,8);
plot(tsamp, xrecon);
axis([0 0.01 -2 2])
xlabel("Time");
ylabel("Amplitude");
title("Over Sampling(fs>2fm)");
