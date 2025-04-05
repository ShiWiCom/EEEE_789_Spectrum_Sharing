% ######################### Dec 11, 2024
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Implementation of curved and planar wavefronts
clear; clc;
% Constants
antenna_Number =2;
Focal_point = 12;
c = 3 * 10^8;
f_c =30 * 10^9;
wave_length = c / f_c
space_Ant_elemets = wave_length/2;
% Antenna axis
antenna_axise = ((1:antenna_Number) - ((antenna_Number + 1) / 2))* space_Ant_elemets;
% antenna_axise = ((1:antenna_Number) )* space_Ant_elemets;

A = space_Ant_elemets^2;%wave_length^2 / (4);
x = linspace(-5,5,150);
z_Planar = linspace(1010,1020,150);
z_Curved = linspace(10,15,150);
WaveFront_Planar = zeros(length(x), length(z_Planar));
WaveFront_Curved = zeros(length(x), length(z_Curved));
for i = 1:length(x)
    for j = 1:length(z_Planar)
        WaveFront_Planar(i,j) = sqrt(1 / ((antenna_Number)^2))*sum(exp(-1i * (2 * pi / wave_length) * ((antenna_axise - x(i)*wave_length).^2/(2*z_Planar(j)*wave_length) + z_Planar(j)*wave_length)));
    end
    for j =1:length(z_Curved)
        WaveFront_Curved(i,j) = sqrt(1 / ((antenna_Number)^2))*sum(exp(-1i * (2 * pi / wave_length) * ((antenna_axise - x(i)*wave_length).^2/(2*z_Curved(j)*wave_length) + z_Curved(j)*wave_length)));
    end
end
 

%% Plot the results
%% Planar WaveFront
figure(1);
set(gcf, 'Color', 'w')

% imagesc(x, (z), real(G_results')); % Heatmap
h = surf(x, z_Planar, real(WaveFront_Planar'), ...
    'EdgeColor', [0.3 0.3 0.3], ...       % Darker edge for better contrast
    'FaceAlpha', 0.9);                    % Slight transparency for depth

shading interp                          % Smooth color transitions
% colormap parula                         % Better color contrast
lighting phong                          % Add smooth lighting
camlight headlight                      % Position the light source
material dull                           % Less shiny surface

xlim([-5,5])
ylim([1010,1015])
set(gca, 'YDir', 'normal'); % Ensure correct orientation
% 
% colorbar;                  % Add colorbar
xlabel('x','FontSize', 14);
ylabel('z','FontSize', 14);
box on;
%% Curved WaveFront
figure(2);
set(gcf, 'Color', 'w')

% imagesc(x, (z), real(G_results')); % Heatmap
h = surf(x, z_Curved, real(WaveFront_Curved'), ...
    'EdgeColor', [0.3 0.3 0.3], ...       % Darker edge for better contrast
    'FaceAlpha', 0.9);                    % Slight transparency for depth

shading interp                          % Smooth color transitions
% colormap parula                         % Better color contrast
lighting phong                          % Add smooth lighting
camlight headlight                      % Position the light source
material dull                           % Less shiny surface
set(h, 'LineWidth', 1.5);

xlim([-5,5])
ylim([10,15])
set(gca, 'YDir', 'normal'); % Ensure correct orientation
xlabel('x','FontSize', 14);
ylabel('z','FontSize', 14);
box on;
