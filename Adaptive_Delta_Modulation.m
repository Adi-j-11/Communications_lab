close all
clear all
clc

% Define sampling period of input signal and desired output sampling period
td = 0.01;  % Original sampling period
ts = 0.02;  % Desired output sampling period

% Generate input signal (sine wave)
t = 0:td:5;         % Time vector
x = 8*sin(2*pi*t);  % Sine wave signal with frequency of 0.5 Hz and amplitude of 8

% Define minimum step size for adaptive delta modulation
delta = 0.1;

% Plot the original input signal
figure(1)
plot(t, x);
xlabel('Time')
ylabel('Amplitude')
title('Original Input Signal')

% Apply adaptive delta modulation
ADMout = adeltamod(x, delta, td, ts);

% Generate time vector for the adaptive delta modulation output
t = linspace(0, length(x)*td, length(ADMout));

% Plot the adaptive delta modulation output
figure;
plot(t, ADMout);
xlabel('Time');
ylabel('Amplitude');
title('Adaptive Delta Modulation Output');

% Explanation of adaptive delta modulation function
% The working of the Advanced Delta Modulator is similar to the regular
% Delta Modulator. The only difference is that the amplitude step
% size is variable and it keeps getting doubled if the previous output/s
% don't seem to 'catch up' with the input signal. This problem is
% referred to as 'Slope overload' in textbooks.
% Usage
% ADMout = adeltamod(sig_in, Delta, td, ts);
% Delta -- min. step size. This will be multiplied 2nX if required
% sig_in -- the signal input, should be a vector
% td -- the original sampling period of the input signal, sig_in
% ts -- the required sampling period for ADM output. Note that it
% should be an integral multiple of the input signal's period.
% If not, it will be rounded up to the nearest integer.
% Function output: ADMout
function ADMout = adeltamod(sig_in, Delta, td, ts)
    % Check if the sampling period ratio is less than 2
    if round(ts/td) < 2
        error('Sampling period ratio is less than 2, cannot perform adaptive delta modulation.');
    end
    
    % Calculate the downsampling factor
    Nfac = round(ts/td);
    
    % Downsample the input signal
    xsig = downsample(sig_in, Nfac);
    
    % Get lengths of downsampled input and original input
    Lxsig = length(xsig);
    Lsig_in = length(sig_in);

    % Initialize the output of adaptive delta modulation
    ADMout = zeros(1, Lsig_in); % Initialize output

    % Initialize variables for adaptive step size
    sum = 0;
    cnt1 = 0;
    cnt2 = 0;

    % Iterate over the downsampled input signal
    for i = 1:Lxsig
        % Check if the downsampled input is greater or less than the current sum
        if xsig(i) > sum
            cnt1 = cnt1 + 1;
            cnt2 = 0;
        elseif xsig(i) < sum
            cnt2 = cnt2 + 1;
            cnt1 = 0;
        end

        % Update the sum based on the current and previous states
        if cnt1 == 0
            sum = sum - 8 * Delta * (cnt2 >= 4);
            sum = sum - 4 * Delta * (cnt2 == 3);
            sum = sum - 2 * Delta * (cnt2 == 2);
            sum = sum - Delta * (cnt2 == 1);
        else
            sum = sum + 8 * Delta * (cnt1 >= 4);
            sum = sum + 4 * Delta * (cnt1 == 3);
            sum = sum + 2 * Delta * (cnt1 == 2);
            sum = sum + Delta * (cnt1 == 1);
        end

        % Assign the sum to the corresponding indices in the output
        ADMout(((i-1)*Nfac + 1):(i*Nfac)) = sum;
    end
end
