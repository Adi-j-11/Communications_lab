v = 1; % voltage level of a bit
R = 1; % Bitrate
T = 1/R; % Bit period
f = 0:0.001*R:2*R; % frequency vector in terms of bit rate
f = f + 1e-10; % Otherwise, sin(0)/0 is undefined

% PSD curves are plotted for Bitrate=1bps and Pulse amplitude=1V

% Unipolar NRZ
s_unipolar_nrz = ((v^2*T/4).*(sin(pi.*f*T)./(pi.*f*T)).^2);
s_unipolar_nrz(1) = s_unipolar_nrz(1) + (v^2/4); % corresponds to an impulse function of weight v^2/4 at f=0 added to s(f) at f=0

% Manchester code
s_manchester = (v^2*T).*((sin(pi.*f*T/2)./(pi.*f*T/2)).^2).*(sin(pi.*f*T/2).^2);

% Polar NRZ
s_polar_nrz = ((v^2*T).*(sin(pi.*f*T)./(pi.*f*T)).^2);

% Bipolar RZ
s_bipolar_rz = (v^2*T/4).*((sin(pi.*f*T/2)./(pi.*f*T/2)).^2).*(sin(pi.*f*T).^2);

% Plotting
figure;
plot(f, (s_unipolar_nrz), '-r', 'LineWidth', 2);
hold on;
plot(f, (s_manchester), '--g', 'LineWidth', 2);
plot(f, (s_polar_nrz), '--b', 'LineWidth', 2);
plot(f, (s_bipolar_rz), '--k', 'LineWidth', 2);

legend('Unipolar NRZ', 'Manchester code', 'Polar NRZ', 'Bipolar RZ/ RZ-AMI');
xlabel('Normalized frequency');
ylabel('Power spectral density (dB)');
title('Power Spectral Densities for Different Modulation Schemes');
grid on;
