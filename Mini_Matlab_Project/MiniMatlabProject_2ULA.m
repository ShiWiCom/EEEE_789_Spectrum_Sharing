% ######################### Feb 23, 2025
% ######################### Mini MATLAB Project: Beampattern and Beamwidth Analysis of Two Uniform Linear Arrays
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Beam Focusing of Two ULAs

clear; clc;

% Constants
N = 64;
F = 20; % Focal Point
c = 3 * 10^8; % Speed of light
f_c = 15 * 10^9; % Carrier frequency
wave_length = c / f_c;
space_Ant_elements = wave_length * 0.5; % Spacing between two adjacent elements
B = [128, 72, 32]; % Distance between two ULAs in multiples of space_Ant_elements
B_bar = ((B + (N - 1)) / 2) * space_Ant_elements; % The center shift of the ULAs from 0 to -B_bar and +B_bar

%% Fraunhofer distance of each ULA
dF = 2 * (space_Ant_elements * N).^2 / wave_length;

%% Antenna axis
xn = ((1:N) - ((N + 1) / 2)) * space_Ant_elements;

%% Define spatial coordinates
x = linspace(-5, 5, 1000);
z = linspace(0, 50, 1000);  

%% Initialize Normalized Array Gain (NAG) matrix
NAG = zeros(length(z), length(x));

% Compute the beampattern 
for k = 1:length(B_bar)
    for i = 1:length(x)
        for j = 1:length(z)
            ULA1 = exp(-1i * (2 * pi / wave_length) * ((xn - B_bar(k) - x(i)).^2 / (2 * z(j)))) .* ...
                   exp(1i * (2 * pi / wave_length) * ((xn - B_bar(k)).^2 / (2 * F)));
            ULA2 = exp(-1i * (2 * pi / wave_length) * ((xn + B_bar(k) - x(i)).^2 / (2 * z(j)))) .* ...
                   exp(1i * (2 * pi / wave_length) * ((xn + B_bar(k)).^2 / (2 * F)));
            NAG(j, i, k) = (1 / ((2 * N)^2)) * (abs(sum(ULA1) + sum(ULA2)).^2);
        end
    end
end

%% Compute Beamwidth with calculated equation
BW_NF = zeros(length(x), length(B_bar));
for i = 1:length(B_bar)
    BW_NF(:, i) = abs(sinc(N * x / (2 * F)) .* (cos(2 * pi * B_bar(i) * x / (wave_length * F)))).^2;
end
BW_UB = (sinc(N * x / (2 * F))).^2;
BW_3dB = 1.77 * F / N;

%% Plot the results
figure;
tiledlayout(3, 2)

% Plot for B = 128
nexttile
imagesc(x, z, NAG(:, :, 1)); 
set(gca, 'YDir', 'normal'); 
colorbar;                 
xlabel('x-axis (m)');
ylabel('z-axis (m)');
title(['Two ULAs separated by B = ', num2str(B(1))]);
hold on;
scatter(0, F, 50, 'r', 'filled'); 
plot([xn(1) - B_bar(1), xn(end) + B_bar(1)], [0, 0], '-m', 'LineWidth', 5); 

nexttile
hold on
plot(x, BW_UB, '-k', 'LineWidth', 2.5, 'DisplayName', 'Upper Bound');
plot(x, BW_NF(:, 1), '-r', 'LineWidth', 2, 'DisplayName', 'BW_2ULA');
plot([-BW_3dB / 2, BW_3dB / 2], [0.5, 0.5], '-.r', 'LineWidth', 3, 'DisplayName', 'Beamwidth');
title('Beampattern along x-axis with N = 64 and B = 128');
xlim([-2, 2]);
xlabel('x_t');
ylabel('Normalized Array Gain');
box on;
legend('show', 'FontSize', 10, 'Location', 'northeast');

% Plot for B = 72
nexttile
imagesc(x, z, NAG(:, :, 2)); 
set(gca, 'YDir', 'normal'); 
colorbar;                 
xlabel('x-axis (m)');
ylabel('z-axis (m)');
title(['Two ULAs separated by B = ', num2str(B(2))]);
hold on;
scatter(0, F, 50, 'r', 'filled'); 
plot([xn(1) - B_bar(2), xn(end) + B_bar(2)], [0, 0], '-m', 'LineWidth', 5); 

nexttile
hold on
plot(x, BW_UB, '-k', 'LineWidth', 2.5, 'DisplayName', 'Upper Bound');
plot(x, BW_NF(:, 2), '-r', 'LineWidth', 2, 'DisplayName', 'BW_2ULA');
plot([-BW_3dB / 2, BW_3dB / 2], [0.5, 0.5], '-.r', 'LineWidth', 3, 'DisplayName', 'Beamwidth');
title('Beampattern along x-axis with N = 64 and B = 72');
xlim([-2, 2]);
xlabel('x_t');
ylabel('Normalized Array Gain');
box on;
legend('show', 'FontSize', 10, 'Location', 'northeast');

% Plot for B = 32
nexttile
imagesc(x, z, NAG(:, :, 3)); 
set(gca, 'YDir', 'normal'); 
colorbar;                 
xlabel('x-axis (m)');
ylabel('z-axis (m)');
title(['Two ULAs separated by B = ', num2str(B(3))]);
hold on;
scatter(0, F, 50, 'r', 'filled'); 
plot([xn(1) - B_bar(3), xn(end) + B_bar(3)], [0, 0], '-m', 'LineWidth', 5); 

nexttile
hold on
plot(x, BW_UB, '-k', 'LineWidth', 2.5, 'DisplayName', 'Upper Bound');
plot(x, BW_NF(:, 3), '-r', 'LineWidth', 2, 'DisplayName', 'BW_2ULA');
plot([-BW_3dB / 2, BW_3dB / 2], [0.5, 0.5], '-.r', 'LineWidth', 3, 'DisplayName', 'Beamwidth');
title('Beampattern along x-axis with N = 64 and B = 32');
xlim([-2, 2]);
xlabel('x_t');
ylabel('Normalized Array Gain');
box on;
legend('show', 'FontSize', 10, 'Location', 'northeast');
