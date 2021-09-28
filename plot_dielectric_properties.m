%%plot dielectric properties
clc;
close all;
clearvars; 
%% initialize
spinal=load('spinalcord.txt');
thyroid=load('thyroid.txt');
fat=load('fat.txt');
muscle=load('muscle.txt');
skin=load('skin.txt');
freq=skin(:,1);
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize', 14)

figure;
plot(freq.*1e-9,skin(:,2));

hold on;
plot(freq.*1e-9,fat(:,2));plot(freq.*1e-9,muscle(:,2));plot(freq.*1e-9,thyroid(:,2));plot(freq.*1e-9,spinal(:,2))
legend('Skin','Fat','Muscle','Thyroid','Spinal Cord')
xlabel('Frequency [GHz]')
ylabel('Dielectric constant')
title('Dielectric constant values for different benign tissues')
grid on;

figure;


plot(freq.*1e-9,skin(:,3));
hold on;
plot(freq.*1e-9,fat(:,3));plot(freq.*1e-9,muscle(:,3));plot(freq.*1e-9,thyroid(:,3));plot(freq.*1e-9,spinal(:,3))
legend('Skin','Fat','Muscle','Thyroid','Spinal Cord')
xlabel('Frequency [GHz]')
ylabel('Electrical conductivity (S/m)')
title('Electrical conductivity values for different benign tissues')
grid on;