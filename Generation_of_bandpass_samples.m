N = 1000;  % number of samples

m = 0;
sgma = 1;

% Initialize arrays to store Gaussian random variables
Xl = zeros(1, N);
X2 = zeros(1, N);

% Generate Gaussian random variables
for i = 1:2:N
    [Xl(i), Xl(i + 1)] = gngauss(m, sgma);
    [X2(i), X2(i + 1)] = gngauss(m, sgma);
end

% Initialize filter coefficients
A = [1, -0.9];
B = 1;

% Filter the Gaussian random variables
Xc = filter(B, A, Xl);
Xs = filter(B, A, X2);

% Bandpass modulation
fc = 1000 / pi;
band_pass_process = zeros(1, N);

for i = 1:N
    band_pass_process(i) = Xc(i) * cos(2 * pi * fc * i) - Xs(i) * sin(2 * pi * fc * i);
end

% Determine the autocorrelation and the spectrum of the bandpass process
M = 50;
bpp_autoco = Rx_est(band_pass_process, M);
bpp_spectrum = fftshift(abs(fft(bpp_autoco)));

% Plotting commands
figure;

subplot(2, 1, 1);
plot(0:M, bpp_autoco);
title('Autocorrelation of Bandpass Process');
xlabel('Lag');
ylabel('Autocorrelation');

subplot(2, 1, 2);
frequencies = linspace(-0.5, 0.5, M + 1);
plot(frequencies, bpp_spectrum);
title('Power Spectral Density of Bandpass Process');
xlabel('Normalized Frequency');
ylabel('Power Spectral Density');


function [gsrvl, gsrv2] = gngauss(m, sgma)
    % Generate two Gaussian random variables with mean m and standard deviation sgma
    gsrvl = m + sgma * randn;
    gsrv2 = m + sgma * randn;
end


% Autocorrelation estimation function
function [Rx] = Rx_est(X, M)
    N = length(X);
    Rx = zeros(1, M + 1);
    
    for m = 1:M + 1
        for n = 1:N - m + 1
            Rx(m) = Rx(m) + X(n) * X(n + m - 1);
        end
        Rx(m) = Rx(m) / (N - m + 1);
    end

end
