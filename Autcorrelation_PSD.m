N = 1000; % Number of samples
M = 75; % Maximum lag for autocorrelation

% Initialize arrays for storing ensemble averages
Rx_av = zeros(1, M + 1);
Sx_av = zeros(1, M + 1);

for j = 1:10 % Take the ensemble average over ten realizations
    % Generate N i.i.d. uniformly distributed random variables between -0.5 and 0.5
    X = rand(1, N) - 1/2;
    
    % Autocorrelation of the realization
    Rx = Rx_est(X, M);
    
    % Power spectrum of the realization
    Sx = fftshift(abs(fft(Rx)));
    
    % Sum of the autocorrelations
    Rx_av = Rx_av + Rx;
    
    % Sum of the spectrums
    Sx_av = Sx_av + Sx;
end

% Average the autocorrelations and spectrums
Rx_av = Rx_av / 10;
Sx_av = Sx_av / 10;

% Plot Ensemble Average Spectrum
figure;
subplot(2, 1, 1);
plot(Sx_av);
title('Ensemble Average Spectrum');
xlabel('Frequency');
ylabel('Magnitude');

% Plot Ensemble Average Autocorrelation
subplot(2, 1, 2);
plot(Rx_av);
title('Ensemble Average Autocorrelation');
xlabel('Lag');
ylabel('Autocorrelation');

% Plot FFT of Ensemble Average Spectrum
figure;
subplot(2, 1, 1);
Sxy = fftshift(abs(fft(Sx_av)));
plot(Sxy);
title('FFT of Ensemble Average Spectrum');
xlabel('Frequency');
ylabel('Magnitude');

% Plot FFT of Autocorrelation
subplot(2, 1, 2);
xy = fftshift(abs(fft(Rx_av)));
plot(xy);
title('FFT of Ensemble Average Autocorrelation');
xlabel('Frequency');
ylabel('Magnitude');

% Ensemble average autocorrelation function estimation
% Function to estimate autocorrelation
function [Rx] = Rx_est(X, M)
    N = length(X);
    Rx = zeros(1, M + 1);
    for m = 1:M+1
        for n = 1:N-m+1
            Rx(m) = Rx(m) + X(n) * X(n+m-1);
        end
        Rx(m) = Rx(m) / (N - m + 1);
    end
end
