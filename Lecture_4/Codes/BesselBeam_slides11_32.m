% ######################### Apr, 5, 2025
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Bessel beam and Bessel beam with obstruction- slides 11 30 31 32
% Constants
antenna_Number = 250;
Focal_point = 20;
c = 3 * 10^8;
f_c = 15 * 10^9;
lambda = c / f_c;
space_Ant_elemets = lambda/2;

% Antenna axis
antenna_axise = ((1:antenna_Number) - ((antenna_Number + 1) / 2)) * space_Ant_elemets;

A = space_Ant_elemets^2;
x = linspace(-5,5,1000); % X-axis range
z = linspace(1,80,1000); % Z-axis range             
Bessel_obstruction = zeros(length(x), length(z));
Bessel = zeros(length(x), length(z));
% Obstruction parameters
block_x_min = -0.15;  % -10 cm
block_x_max = 0.15;   % +10 cm
block_Loc_z = 25;    % Obstruction location at z = 25

for i = 1:length(x)
    for j = 1:length(z)
        % Create a line for the obstruction to block the signal
        line = block_Loc_z*((x(i)-antenna_axise)/z(j))+antenna_axise;
        % Check if the point is in the blocked region
        f = exp(1i * ((2 * pi / lambda) * ((abs(antenna_axise) * sin(pi * 0.007))))) .* ...
            exp(-1i * (2 * pi / lambda) * ((antenna_axise - x(i)).^2 / (2 * z(j)) )); 
        for n=1:length(f)
            if (line(n) >= block_x_min && line(n) <= block_x_max && z(j) >= block_Loc_z)
                f_b(n) = 0; % Set field to zero in blocked region
            else
                f_b(n) = f(n);
            end
        end
        % Nearfield-Bessel Beam 
        Bessel(i, j) = (1 / ((antenna_Number)^2)) * (abs(sum(f))^2);
        % Nearfield-Bessel Beam with Obstruction
        Bessel_obstruction(i, j) = (1 / ((antenna_Number)^2)) * (abs(sum(f_b))^2);
    end
end

%% Plot results
figure(1);
hold on;
set(gcf, 'Color', 'w');
imagesc(x, z, Bessel'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation
colorbar;                  
xlabel('x-axis (m)','FontSize', 14);
ylabel('z-axis (m)','FontSize', 14);
xlim([-5,5])
ylim([1,80])
title('Bessel Beam','FontSize', 14);
% plot the antenna aperture
plot([-antenna_Number*wave_length/4,antenna_Number*wave_length/4] , [1,1],'-m','LineWidth',5); 

%%
figure(2);
hold on;
set(gcf, 'Color', 'w');
imagesc(x, z, Bessel_obstruction'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation
colorbar;                  
xlabel('x-axis (m)','FontSize', 14);
ylabel('z-axis (m)','FontSize', 14);
xlim([-5,5])
ylim([1 ,80])
title('Bessel Beam with Obstruction','FontSize', 14);
% Add obstruction marker
rectangle('Position', [block_x_min, block_z, (block_x_max - block_x_min), 1], ...
          'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', [1, 0, 0, 0.3]); % Red transparent block
% plot the antenna aperture
plot([-antenna_Number*wave_length/4,antenna_Number*wave_length/4] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)

