% Define time vector
t = 0:0.01:10;

% Generate original signal (sine wave)
a = sin(t);

% Perform uniform PCM for different number of levels
[sqnr8, aquan8, code8] = uniform_pcm(a, 8);
[sqnr16, aquan16, code16] = uniform_pcm(a, 16);

% Display SQNR for N = 8 and N = 16
disp('SQNR for N = 8:');
disp(sqnr8);

disp('SQNR for N = 16:');
disp(sqnr16);

% Plot original signal and quantized versions
plot(t, a, '-', t, aquan8, '-.', t, aquan16, '-', t, zeros(1, length(t)));
legend('Original Signal', '8 Level Quantized Signal', '16 Level Quantized Signal');
title('Signal and Quantized Versions');
xlabel('Time');
ylabel('Amplitude');

% Function for uniform PCM
function [sqnr, a_quan, code] = uniform_pcm(a, n)
    % Normalize the signal
    amax = max(abs(a));
    a_quan = a / amax;
    b_quan = a_quan;
    
    % Define quantization levels
    d = 2 / n;
    q = d * (0:n-1) - (n-1)/2*d;

    % Quantization process
    for i = 1:n
        indices = (q(i) - d/2 <= a_quan) & (a_quan <= q(i) + d/2);
        a_quan(indices) = q(i) * ones(1, sum(indices));
        b_quan(a_quan == q(i)) = (i-1) * ones(1, sum(a_quan == q(i)));
    end

    % Denormalize the quantized signal
    a_quan = a_quan * amax;
    
    % Convert quantized signal to binary code
    num_bits = ceil(log2(n));
    code = zeros(length(a), num_bits);
    for i = 1:length(a)
        for j = num_bits:-1:0
            if fix(b_quan(i) / 2^j) == 1
                code(i, (num_bits - j) + 1) = 1;
                b_quan(i) = b_quan(i) - 2^j;
            end
        end
    end

    % Calculate SQNR
    sqnr = 20 * log10(norm(a) / norm(a - a_quan));
end
