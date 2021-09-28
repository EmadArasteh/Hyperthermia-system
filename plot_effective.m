%%plot effective properties
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

f=[0.024 0.047 0.15 0.47 0.02 0.02 0.19 0.04 0.039]; % air bone fat muscle osephageous spinal_cord skin Thyroid Trachea
eps_spinal=spinal(:,2);
eps_thyroid=thyroid(:,2);
eps_fat=fat(:,2);
eps_muscle=muscle(:,2);
eps_skin=skin(:,2);
eps_air=ones(numel(eps_spinal),1);
eps_bone=eps_air.*15.28;
eps_ose=eps_air.*77.9;
eps_trachea=eps_air.*52.96;
eps=[eps_air eps_bone eps_fat eps_muscle eps_ose eps_spinal eps_skin eps_thyroid eps_trachea];


sigma_spinal=spinal(:,3);
sigma_thyroid=thyroid(:,3);
sigma_fat=fat(:,3);
sigma_muscle=muscle(:,3);
sigma_skin=skin(:,3);
sigma_air=ones(size(sigma_spinal));
sigma_bone=sigma_air.*0.064;
sigma_ose=sigma_air.*0.9;
sigma_trachea=sigma_air.*0.548;
sigma_air=sigma_air.*0;
sigma=[sigma_air sigma_bone sigma_fat sigma_muscle sigma_ose sigma_spinal sigma_skin sigma_thyroid sigma_trachea];


for i=1:numel(freq)
eps_effective(i)=sum(eps(i,:).*f);
sigma_effective(i)=sum(sigma(i,:).*f);
end

set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize', 14)

figure;
plot(freq.*1e-9,eps_effective)
xlabel('Frequency [GHz]')
ylabel('Dielectric constant')
title('Dielectric constant values for effective tissue')
grid on;
axis([1 2 40 42])


figure;
plot(freq.*1e-9,sigma_effective)
axis([1 2 0.7 1.1])
xlabel('Frequency [GHz]')
ylabel('Electrical conductivity (S/m)')
title('Electrical conductivity values for effective tissue')
grid on;

effective=zeros(3,numel(freq));
effective(1,:)=freq;
effective(2,:)=eps_effective;
effective(3,:)=sigma_effective;
save('effective.mat')