%%%undefined gelatine for Thyroid produceing
clc;
close all;
clearvars;
clc;
clearvars;
close all;
%% General initialization

eps_0=8.854e-12;
freq=linspace(.5e9,5e9,47*91*51);
w=freq.*2*pi;
w=w';
max_error_eps=.05;
 max_error_sigma=.06;
max_error_Cp=.07;
max_error_K=.07;

 max_error_tot=(max_error_eps+max_error_sigma+max_error_Cp+max_error_K)*.5;
%%%%%%%%%%%%%%%%Density
Rho_DW=1; %g/ml
Rho_Oil=4.58/5; %g/ml
Rho_salt=0.52;
Rho_Gel=3.56/5; %g/ml 
Rho_Surfactant=1.1; %gram/mole
%%%%%%%%%%%%%%%%%%%%%%Molar_mass  gram/m
n_DW=18.015; %g/mole H20
n_cannola=18.015;% g/mole
n_Gelatine=288.38; %g/mole NaC12H25SO4

%%%%%%%%%%%%%%%Cp J/g^oc
Cp_Oil=1.67;
Cp_DW=4.184; %J/g^oc
Cp_salt=0;
Cp_Gelat=4.5;
Cp_Thyroid=3.6;
Cp_Surfactant=0;
%%%%%%%%%%%k
K_DW=.6;
K_spinal_cord=.51;
K_tumor=.5;
K_Oil=.17;
K_Gel=.6;
K_Thyroid=.52;
K_Surfactant=0;




%% start----gelatine as unknown
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55Gelatine soluted
% syms x sigma_Gelatine%ratio of gel at 1.5 GHz;
% ratio_of_gelat=solve('(x*10.842568513702741)+((1-x)*(77.986522345783911))-42.9042=0;',x);
% eps_imag_gel_sol=11.9069;
% f=1.5e9;
% sigma_gel_sol=2*pi*f*eps_0*eps_imag_gel_sol;
% sigma_Gelat=solve('(0.522*sigma_Gelatine)+((1-0.522)*(0.494007621824038))-0.9935;',sigma_Gelatine);
sigma_Gelat=1.4508895723527008352490421455939;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
sigma_Surfactant=25;
eps=[77.986522345783911   2.754934882341327  10.842568513702741 0];
sigm=[0.494007621824038   0.012806404700502   sigma_Gelat sigma_Surfactant];
epsilon_prime_Thyroid=58.59683581386818;
sigma_Thyroid=1.322141371156771;

Cp=[Cp_DW Cp_Oil Cp_Gelat Cp_Surfactant];
K=[K_DW K_Oil K_Gel K_Surfactant];

M_DW=23;
change=.2;
M_Gel=4.3-change:.1:4.3+change;
M_Oil=3.1-change:.1:3.1+change;

for M_Surfactant=1:.1:1;

M_Surfactant

for k=1:numel(M_Oil);
    
  for j=1:numel(M_Gel)
 
V_Gel=M_Gel(j)/Rho_Gel;
V_DW=M_DW/Rho_DW;
V_Oil=M_Oil(k)/Rho_Oil;       
V_Surfactant=M_Surfactant/Rho_Surfactant;
V=[V_DW V_Oil V_Gel];
M=[M_DW M_Oil(k) M_Gel(j) M_Surfactant];
f=[V_DW/sum(V) V_Oil/sum(V) V_Gel/sum(V) V_Surfactant/sum(V)];




error_eps(j,k)=sum(f.*eps)-epsilon_prime_Thyroid(i);
% sum(f.*(eps-epsilon_prime_Thyroid(i))./(eps+(3.*epsilon_prime_Thyroid(i))));
error_sigma(j,k)=sum(f.*sigm)-sigma_Thyroid(i);
% (sigm-sigma_Thyroid(i))./(sigm+(3.*sigma_Thyroid(i))));
error_Cp(j,k)=(Cp_Thyroid)-sum((M.*Cp)./sum(M));
error_K(j,k)=K_Thyroid-sum((M.*K)./sum(M));



  end
  
   
   
end

[idx1,idx2]=find(abs(error_eps)<=max_error_eps.*epsilon_prime_Thyroid(i) & abs(error_sigma)<=max_error_sigma.*sigma_Thyroid(i) & abs(error_Cp)<=max_error_Cp.*Cp_Thyroid & abs(error_K)<=max_error_K.*K_Thyroid & (abs(error_eps./epsilon_prime_Thyroid(i))+abs(error_sigma./sigma_Thyroid(i))+abs(error_Cp./Cp_Thyroid)+abs(error_K./K_Thyroid))<=max_error_tot);
(M_Gel(idx1));
 M_Oil(idx2);
 eerror_of_cp=error_Cp(idx1,idx2)
 eeror_of_K=error_K(idx1,idx2)
end
