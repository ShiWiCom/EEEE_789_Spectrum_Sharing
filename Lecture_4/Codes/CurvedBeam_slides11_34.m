% ######################### Apr, 5, 2025
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Curved beam - slides 11 - 
clear; clc;
% Constants
antenna_Number = 250;
Focal_point = 12;
c = 3 * 10^8;
f_c = 15 * 10^9;
wave_length = c / f_c;
space_Ant_elemets = wave_length/2;
% Antenna axis
antenna_axise = ((1:antenna_Number) - ((antenna_Number + 1) / 2))* space_Ant_elemets;

A = space_Ant_elemets^2;%wave_length^2 / (4);
x = linspace(-5,5,1000);
z = linspace(1,20,1000);  
a = -0.0125;
b = 0.0025;
c = 1;
d=antenna_axise(1);
G_results = zeros(length(x), length(z));
% gz =-0.0125*(z-1).^2+0.0025;
% dgz = -0.025*(z-1);
Zg =sqrt(abs((antenna_axise -(a*c^2+b+d))/(-a)));
Xg = a*(Zg-c).^2+b+d;

for i = 1:length(x)
    for j = 1:length(z)
        G_results(i, j) =(1 / ((antenna_Number)^2)) * (abs(sum(exp(1i * ((2 * pi / wave_length) * ((antenna_axise-Xg).^2 ./ (2 * Zg)))).* ...
        exp(-1i * (2 * pi / wave_length) * ((antenna_axise - x(i)).^2 / (2 * z(j))))))^2);
  
    end
end 

%% Plot the results
figure;
set(gcf, 'Color', 'w');
imagesc(x, (z), G_results'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation

colorbar;                  % Add colorbar
xlabel('x-axis (m)','FontSize', 24);
ylabel('z-axis (m)','FontSize', 24);
%yscale log;
title('Curved Beam','FontSize', 20);
% Add markers (example: red dot and magenta square)
hold on;
% scatter(0, 30, 50, 'r', 'filled'); % Red dot at (0, 50)
plot([-antenna_Number*wave_length/4,antenna_Number*wave_length/4] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)
