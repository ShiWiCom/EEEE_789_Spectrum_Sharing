% ######################### Dec 13, 2024
% ######################### Written by Shima Mashhadi, PhD student at Rochester Institute of Technology, supervised by Professor Alireza Vahid.
% ######################### Near-field beamfocucing design 
% ######################### EEEE-789 Spectrum Sharing
% ######################### Beamwidth in nearfield- slide 19
clear; clc;
% Constants
antenna_Number =[150; 250];
Focal_point = 20;
l = [1,2];
c = 3* 10^8;
f_c = 15* 10^9;
lambda = c / f_c;
space_Ant_elemets = lambda/2;
x = linspace(-0.6,0.6,16*1000);

z = Focal_point;  
%% Beam Depth
dF = 2*(space_Ant_elemets.*antenna_Number).^2/lambda;
BD = dF*Focal_point/(dF-10*Focal_point)- dF*Focal_point/(dF+10*Focal_point);
%% Fresnel approximation 
BW_FA = abs(sinc(antenna_Number.*x/(2*z))).^2;

%% Upper Bound
BW_UB = (sinc(antenna_Number.*x/(2*z))).^2;

%%
BW_3dB = 1.77*Focal_point./antenna_Number;
%% Plot results

figure(1);
hold on;

plot(x,BW_UB(1,:),'-k','LineWidth',2.5);%,'DisplayName');%'Upper bound');
plot([-BW_3dB(1)/2 , BW_3dB(1)/2], [0.5 , 0.5],'-.r','LineWidth',3, 'DisplayName', 'BeamWidth');
title('The Beamwidth of a ULA along the x-axis with N = 150');
%xlabel('x-axis(m)', 'FontSize', 18);
xlabel('x_t')
ylabel('Normalized Antenna gain')
box on;
%%
figure(2);
hold on;
plot(x,BW_UB(2,:),'-k','LineWidth',2.5);%,'DisplayName');%'Upper bound');
plot([-BW_3dB(2)/2 , BW_3dB(2)/2], [0.5 , 0.5],'-.r','LineWidth',3, 'DisplayName', 'BeamWidth');
title('The Beamwidth of a ULA along the x-axis with N = 250');
xlabel('x_t')
ylabel('Normalized Antenna gain')
box on;
% legend('show', 'FontSize', 10, 'Location', 'northeast');
%ylabel('Normalized gain' ,'FontSize', 18);
