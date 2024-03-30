N = 1000;   % The maximum value of n
M = 50;

Rxav = zeros(1, M + 1);
Ryav = zeros(1, M + 1);
Sxav = zeros(1, M + 1);
Syav = zeros(1, M + 1);

% Perform ensemble averaging over ten realizations
for i = 1:10
    X = rand(1, N) - (1/2);
    Y(1) = 0;
    
    % Generate Y(n) using a first-order autoregressive process
    for n = 2:N
        Y(n) = 0.9 * Y(n - 1) + X(n);
    end
    
    % Compute autocorrelations of X(n) and Y(n)
    Rx = Rx_est(X, M);
    Ry = Rx_est(Y, M);
    
    % Compute power spectra of X(n) and Y(n)
    Sx = fftshift(abs(fft(Rx)));
    Sy = fftshift(abs(fft(Ry)));
    
    % Accumulate results for ensemble averaging
    Rxav = Rxav + Rx;
    Ryav = Ryav + Ry;
    Sxav = Sxav + Sx;
    Syav = Syav + Sy;
end

% Average over the ten realizations
Rxav = Rxav / 10;
Ryav = Ryav / 10;
Sxav = Sxav / 10;
Syav = Syav / 10;

% Display results
disp('Ensemble Average Autocorrelation of X:');
disp(Rxav);
disp('Ensemble Average Autocorrelation of Y:');
disp(Ryav);
disp('Ensemble Average Power Spectrum of X:');
disp(Sxav);
disp('Ensemble Average Power Spectrum of Y:');
disp(Syav);

% Plot Ensemble Average Autocorrelations and Power Spectra
figure;

% Plot Ensemble Average Autocorrelations
subplot(2, 1, 1);
plot(0:M, Rxav, 'DisplayName', 'X(n)');
title('Ensemble Average Autocorrelations');
xlabel('Lag');
ylabel('Autocorrelation');
legend;
subplot(2, 1, 2);
plot(0:M, Ryav, 'DisplayName', 'Y(n)');
title('Ensemble Average Autocorrelations');
xlabel('Lag');
ylabel('Autocorrelation');
legend;

figure;

% Plot Ensemble Average Power Spectra
subplot(2, 1, 1);
frequencies = linspace(-0.5, 0.5, M + 1);
plot(frequencies, Sxav, 'DisplayName', 'X(n)');
title('Ensemble Average Power Spectra');
xlabel('Normalized Frequency');
ylabel('Power Spectral Density');
legend;
subplot(2, 1, 2);
plot(frequencies, Syav, 'DisplayName', 'Y(n)');
title('Ensemble Average Power Spectra');
xlabel('Normalized Frequency');
ylabel('Power Spectral Density');
legend;

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
