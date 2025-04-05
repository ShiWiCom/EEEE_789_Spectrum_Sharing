% ######################### Apr, 5, 2025
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Farfield beamforming
clear; clc;
% Constants
antenna_Number = 50;
Focal_point = 20;
c = 3 * 10^8;
f_c = 15 * 10^9;
lambda = c / f_c;
space_Ant_elemets =0.5*lambda;
% Antenna axis
antenna_axise = ((1:antenna_Number) - ((antenna_Number + 1) / 2))* space_Ant_elemets;

dF = 2*(space_Ant_elemets*(antenna_Number)).^2/lambda;

A = space_Ant_elemets^2;%lambda^2 / (4);

x = linspace(-5,5,2000);
z = linspace(0,100,1000); 

G_results_0 = zeros(length(x),length(z));
G_results_45 = zeros(length(x),length(z));
for i = 1:length(x)
    for j = 1:length(z)
         f1 = exp(1i * ((2 * pi / lambda) * (((antenna_axise)*sin(pi*0))))).* ...
                exp(-1i * (2 * pi / lambda) * ((antenna_axise - x(i)).^2 / (2 * z(j))));

         G_results_0(i,j) = (1 / ((antenna_Number)^2)) * (abs(sum(f1)).^2);
         
         
         f2 = exp(1i * ((2 * pi / lambda) * (((antenna_axise)*sin(pi/50))))).* ...
                exp(-1i * (2 * pi / lambda) * ((antenna_axise - x(i)).^2 / (2 * z(j))));

         G_results_45(i,j) = (1 / ((antenna_Number)^2)) * (abs(sum(f2)).^2);
    end
end
 

%% Plot the results
figure(1);
hold on;
set(gcf, 'Color', 'w');
imagesc(x,z,G_results_0');
colorbar; 
xlabel('x-axis (m)','FontSize', 14);
ylabel('z-axis (m)','FontSize', 14);

xlim([-5,5])
ylim([0,100])
plot([-antenna_Number*space_Ant_elemets/2,antenna_Number*space_Ant_elemets/2] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)

%%
figure(2);
hold on;
set(gcf, 'Color', 'w');
imagesc(x,z,G_results_45');
colorbar; 
xlabel('x-axis (m)','FontSize', 14);
ylabel('z-axis (m)','FontSize', 14);

xlim([-5,5])
ylim([0,100])
plot([-antenna_Number*space_Ant_elemets/2,antenna_Number*space_Ant_elemets/2] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)




















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