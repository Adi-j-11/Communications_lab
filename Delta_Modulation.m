clc;
clear all;
close all;

a = 2;
t = 0:2*pi/50:2*pi; % Signal Generation
x = a * sin(t);
l = length(x);

% Plot the original analog signal
subplot(4,1,1);
plot(x, 'r');
title('Analog Signal');
xlabel('Sample');
ylabel('Amplitude');

% Delta modulation with step size = delta
subplot(4,1,2);
delta = 0.2;
xn1 = delta_modulation(x, delta);
stairs(xn1);
title(['Delta Modulation - Step Size: ', num2str(delta)]);
xlabel('Sample');
ylabel('Quantized Level');

% Delta modulation with step size = 2 * delta
subplot(4,1,3);
delta2 = 2 * delta;
xn2 = delta_modulation(x, delta2);
stairs(xn2);
title(['Delta Modulation - Step Size: ', num2str(delta2)]);
xlabel('Sample');
ylabel('Quantized Level');

subplot(4,1,4);
delta3 = 3 * delta;
xn2 = delta_modulation(x, delta3);
stairs(xn2);
title(['Delta Modulation - Step Size: ', num2str(delta3)]);
xlabel('Sample');
ylabel('Quantized Level');
sgtitle('Effect of Different Step Sizes in Delta Modulation');
function xn = delta_modulation(x, delta)
    l = length(x);
    xn = zeros(1, l);

    for i = 1:l
        if x(i) > xn(i)
            xn(i+1) = xn(i) + delta;
        else
            xn(i+1) = xn(i) - delta;
        end
    end
end
