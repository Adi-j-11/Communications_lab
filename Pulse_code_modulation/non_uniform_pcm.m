% Generate time vector from 0 to 10 with step 0.01
t = 0:0.01:10;

% Generate a sinusoidal signal
a = sin(t);

% Apply mu-law PCM and plot original and quantized signals
[sqnr, aquan] = mula_pcm(a, 16, 255);
display(sqnr);
plot(a);
hold on;
plot(aquan);
legend('Original Signal', 'Quantized Signal');

% Function to perform mu-law PCM compression followed by uniform PCM quantization
function [sqnr, a_quan, code] = mula_pcm(a, n, mu)
    % Apply mu-law compression to the input signal
    [y, maximum] = mulaw(a, mu);
    
    % Perform uniform PCM quantization on the compressed signal
    [sqnr, y_q, code] = u_pcm(y, n);
    
    % Apply inverse mu-law compression to obtain the quantized signal
    a_quan = invmulaw(y_q, mu);
    
    % Rescale the quantized signal
    a_quan = maximum * a_quan;
    
    % Calculate the SQNR of the quantized signal
    sqnr = 20 * log10(norm(a) / norm(a - a_quan));
end

% Function to perform mu-law compression
function [y, a] = mulaw(x, mu)
    % Normalize input signal
    a = max(abs(x));
    % Apply mu-law compression
    y = (log(1 + mu * abs(x/a)) / log(1 + mu)) .* sign(x);
end

% Function to perform inverse mu-law compression
function x = invmulaw(y, mu)
    % Apply inverse mu-law compression
    x = (((1 + mu).^(abs(y)) - 1) ./ mu) .* sign(y);
end

% Function to perform uniform PCM quantization
function [sqnr, a_quan] = u_pcm(a, n)
    % Normalize input signal
    amax = max(abs(a));
    a_quan = a / amax;
  
    % Calculate quantization step size
    d = 2 / n;
    % Calculate quantization levels
    q = d * (0:n-1) - (n-1) / 2 * d;

    % Quantize the input signal
    for i = 1:n
        % Find indices of elements within the quantization range
        indices = (q(i) - d/2 <= a_quan) & (a_quan <= q(i) + d/2);
        % Replace elements within the quantization range with the quantization level
        a_quan(indices) = q(i) * ones(1, sum(indices));
    end
    
    % Rescale the quantized signal
    a_quan = a_quan * amax;
    % Calculate the SQNR of the quantized signal
    sqnr = 20 * log10(norm(a) / norm(a - a_quan));
end
