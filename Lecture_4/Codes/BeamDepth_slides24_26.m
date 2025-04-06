% ######################### Apr, 5, 2025
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### EEEE-789 Spectrum Sharing
% ######################### Beamdepth in nearfield- slide 24 and 26
clc; clear; close all;

% Define parameter ranges
x_values = 0:0.001:10; % x = d_F / (8*z_eff)

% Define fixed constants
lambda = 1;  % Normalization (not affecting the core result)

results = zeros(1, length(x_values));
% Define function handles for fresnel integrands
cos_func = @(t) cos((pi/2) * t.^2);
sin_func = @(t) sin((pi/2) * t.^2);
%% square array
 for j = 1:length(x_values)
       
        x = x_values(j);
        C = integral(cos_func, 0, sqrt(x));
        S = integral(sin_func, 0, sqrt(x));
        results(j) = (1 / x)^2 * (C^2 + S^2)^2;
 end

 %% linear array
 % for j = 1:length(F_values)
 %        x = x_values(j);
 % 
 %        x = x_values(j);
 %        C = integral(cos_func, 0, sqrt(x));
 %        S = integral(sin_func, 0, sqrt(x));
 %        % Compute the final equation
 %        results(j) = (1 / x) * (C^2 + S^2);
 %        
 % end
%% Plot the results
figure(1);
hold on;
plot(x_values, results,'-k','LineWidth',2.5);
plot(1.25, 0.5, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% surf(F_values, results, 'EdgeColor', 'none');
text(1.25, 0.5, '(1.25, 0.5)', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);

xlabel('x');
ylabel('A(x)');
box on;
%% slide 26
% Constants
antenna_Number = [150,250,500];
Focal_Point = [20,30,60]; % Focal point
c = 3 * 10^8;
f_c = 15 * 10^9;
lambda = c / f_c;
space_Ant_elemets = lambda/2;
z = 0:0.1:100;%linspace(1,50,500);             
G_results = zeros(length(Focal_Point),length(antenna_Number), length(z));

for n = 1:length(Focal_Point)
    F = Focal_Point(n);
    for i = 1:length(antenna_Number)
        antenna_axise = ((1:antenna_Number(i)) - ((antenna_Number(i) + 1) / 2))* space_Ant_elemets;
        for j = 1:length(z)
    
            f =exp(1i * (2 * pi / lambda) * ((antenna_axise).^2 / (2 * F))).* ...
              exp(-1i * (2 * pi / lambda) * ((antenna_axise).^2 / (2 * z(j))));
        
            G_results(n,i, j) = (1 / ((antenna_Number(i))^2)) * (abs(sum(f)).^2);
          
        end
    end
end
%%
% dF = 2*(space_Ant_elemets*(antenna_Number)).^2/lambda;
% BD = (dF.*Focal_point./(dF-6.952*Focal_point))-(dF.*Focal_point./(dF+6.952*Focal_point));

%% Plot the results
figure(2);
hold on;
plot(z, reshape(G_results(1,1,:),[1, length(z)]),'-k','LineWidth',2.5,'DisplayName','N=150'); % Heatmap
plot(z,  reshape(G_results(1,2,:),[1, length(z)]),'-b','LineWidth',2.5,'DisplayName','N=250'); 
plot(z, reshape(G_results(1,3,:),[1, length(z)]),'-r','LineWidth',2.5,'DisplayName','N=500'); 
xlabel('z-axis (m)','FontSize', 20);
ylabel('Normalized Array Gain','FontSize', 20);
%yscale log;
title('Beam Depth, Focal Point = 20m','FontSize', 14);
% Add markers (example: red dot and magenta square)
legend;
box on;
%%
figure(3);
hold on;
plot(z, reshape(G_results(1,1,:),[1, length(z)]),'-k','LineWidth',2.5,'DisplayName','N=150, F=20'); % Heatmap
plot(z, reshape(G_results(2,1,:),[1, length(z)]),'-b','LineWidth',2.5,'DisplayName','N=150, F=30'); 
plot(z, reshape(G_results(3,1,:),[1, length(z)]),'-r','LineWidth',2.5,'DisplayName','N=150, F=60'); 
xlabel('z-axis (m)','FontSize', 20);
ylabel('Normalized Array Gain','FontSize', 20);
%yscale log;
title('Beam Depth','FontSize', 14);
% Add markers (example: red dot and magenta square)
box on;
legend;