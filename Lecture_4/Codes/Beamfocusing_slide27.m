% ######################### Apr, 5, 2025
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Nearfield beamfocusing-slide 27
clear; clc;
% Constants
antenna_Number = [150, 250, 500];
F = 20;% Focal Point
c = 3 * 10^8;
f_c = 15 * 10^9;
lambda = c / f_c;
space_Ant_elemets =0.5*lambda;
% antenna_axise = ((1:antenna_Number) )* space_Ant_elemets;
dF = 2*(space_Ant_elemets*(antenna_Number)).^2/lambda;

x = linspace(-3,3,2000);
z = linspace(0,50,1000);             
G_results = zeros(length(antenna_Number),length(x), length(z));
for m = 1: length(antenna_Number)
    % Antenna axis
    antenna_axise = ((1:antenna_Number(m)) - ((antenna_Number(m) + 1) / 2))* space_Ant_elemets;
    for i = 1:length(x)
        for j = 1:length(z)
             f =exp(1i * (2 * pi / lambda) * (antenna_axise).^2 / (2 * F)).* ...
                exp(-1i * (2 * pi / lambda) * (antenna_axise - x(i)).^2 / (2 * z(j)));
 
             G_results(m, i, j) = (1 / ((antenna_Number(m))^2)) * (abs(sum(f)).^2);
        end
    end
end
 

%% Plot the results

figure(1);
set(gcf, 'Color', 'w');
hold on;
imagesc(x,z, reshape(G_results(1,:,:),[length(x), length(z)])'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation

colorbar;                  % Add colorbar
xlabel('x-axis (m)','FontSize', 20);
ylabel('z-axis (m)','FontSize', 20);
xlim([-3,3])
ylim([1,50])
title('Beam focoucing M=150','FontSize', 14);
plot([-antenna_Number(1)*space_Ant_elemets/2,antenna_Number(1)*space_Ant_elemets/2] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)
%%
figure(2);
set(gcf, 'Color', 'w');
hold on;
imagesc(x,z, reshape(G_results(2,:,:),[length(x), length(z)])'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation

colorbar;                  % Add colorbar
xlabel('x-axis (m)','FontSize', 20);
ylabel('z-axis (m)','FontSize', 20);
xlim([-3,3])
ylim([1,50])
title('Beam focoucing M=250','FontSize', 14);
plot([-antenna_Number(2)*space_Ant_elemets/2,antenna_Number(2)*space_Ant_elemets/2] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)
%%
figure(3);
set(gcf, 'Color', 'w');
hold on;
imagesc(x,z, reshape(G_results(3,:,:),[length(x), length(z)])'); % Heatmap
set(gca, 'YDir', 'normal'); % Ensure correct orientation

colorbar;                  % Add colorbar
xlabel('x-axis (m)','FontSize', 20);
ylabel('z-axis (m)','FontSize', 20);
xlim([-3,3])
ylim([1,50])
title('Beam focoucing M=500','FontSize', 14);
plot([-antenna_Number(3)*space_Ant_elemets/2,antenna_Number(3)*space_Ant_elemets/2] , [1,1],'-m','LineWidth',5); % Magenta square at (0, 0)
