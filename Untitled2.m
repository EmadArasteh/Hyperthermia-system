%%3d fdtd with truncation for fat
clc;
close all;
clearvars;
%% initialization
c_0=3e8;

f=1e9;
T0=1/f;
T_simulation=6.25*T0;
lambda=c_0/f;
E_0=500;
%%%%%%%%%%%%%%%%%%%%%%MESH
len_body_x=300e-3;
len_body_y=len_body_x;
len_body_z=len_body_x;
factor_air=10;
len_air_x=lambda*factor_air/4;
len_air_y=len_air_x;
len_air_z=len_air_x;

len_tot_x=2*len_air_x+len_body_x;
len_tot_y=2*len_air_z+len_body_y;
len_tot_z=2*len_air_z+len_body_z;

delta_x=lambda/20;
delta_y=delta_x;
delta_z=delta_x;
delta=delta_x;

Nx=round(len_tot_x/delta_x);
Ny=round(len_tot_y/delta_y);
Nz=round(len_tot_z/delta_z);


num_node_air_x=round(len_air_x/delta_x);
num_node_air_y=round(len_air_y/delta_y);
num_node_air_z=round(len_air_z/delta_z);




mu_r=zeros(Nx,Ny,Nz);
sigma=zeros(Nx,Ny,Nz);
Rho=zeros(Nx,Ny,Nz);
eps_r=zeros(Nx,Ny,Nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%Dielectric



Rho_fat=911;
Rho_min=Rho_fat;
eps_r_fat=10;
eps_r_max=eps_r_fat;
sigma_fat=0.17;
mu_fat=1;

mu_air=1;
sigma_air=0;
Rho_air=1;
eps_r_air=1;


mu_r(:,:,:)=mu_fat;
mu_r(1:num_node_air_x/2,:,:)=mu_air;
mu_r(end:end-num_node_air_x/2+1,:,:)=mu_air;
mu_r(:,1:num_node_air_y/2,:)=mu_air;
mu_r(:,end:end-num_node_air_y/2+1,:)=mu_air;
mu_r(:,:,1:num_node_air_z/2)=mu_air;
mu_r(:,:,end:end-num_node_air_z/2+1)=mu_air;



eps_r(:,:,:)=eps_r_fat;
eps_r(1:num_node_air_x/2,:,:)=eps_r_air;
eps_r(end:end-num_node_air_x/2+1,:,:)=eps_r_air;
eps_r(:,1:num_node_air_y/2,:)=eps_r_air;
eps_r(:,end:end-num_node_air_y/2+1,:)=eps_r_air;
eps_r(:,:,1:num_node_air_z/2)=eps_r_air;
eps_r(:,:,end:end-num_node_air_z/2+1)=eps_r_air;



sigma(:,:,:)=sigma_fat;
sigma(1:num_node_air_x/2,:,:)=sigma_air;
sigma(end:end-num_node_air_x/2+1,:,:)=sigma_air;
sigma(:,1:num_node_air_y/2,:)=sigma_air;
sigma(:,end:end-num_node_air_y/2+1,:)=sigma_air;
sigma(:,:,1:num_node_air_z/2)=sigma_air;
sigma(:,:,end:end-num_node_air_z/2+1)=sigma_air;



Rho(:,:,:)=Rho_fat;
Rho(1:num_node_air_x/2,:,:)=Rho_air;
Rho(end:end-num_node_air_x/2+1,:,:)=Rho_air;
Rho(:,1:num_node_air_y/2,:)=Rho_air;
Rho(:,end:end-num_node_air_y/2+1,:)=Rho_air;
Rho(:,:,1:num_node_air_z/2)=Rho_air;
Rho(:,:,end:end-num_node_air_z/2+1)=Rho_air;



%%%%%%%%%%%%Time
u_max=c_0/sqrt(eps_r_max);
m=ceil(sqrt(eps_r_max));
delta_t1=(sqrt(eps_r_max)*delta_x)/(m*c_0);
delta_t2=delta/(sqrt(3)*u_max);
delta_t=min(delta_t1,delta_t2);
Time_vec=0:delta_t:T_simulation;


%%%%%%%%ending initialization
epsilon_0=(10^(-9)/(36*pi));
mu_0=4.*pi*(10^(-7));
epsilon=epsilon_0.*eps_r;
mu=mu_0.*mu_r;
alfa=delta_t./(delta_x.*mu);
beta=1-((sigma.*delta_t)./epsilon);
gamma=delta_t./(delta_x.*epsilon);
Hx=zeros(Nx,Ny,Nz);
Hy=zeros(Nx,Ny,Nz);
Hz=zeros(Nx,Ny,Nz);
Ex=zeros(Nx,Ny,Nz);
Ey=zeros(Nx,Ny,Nz);
Ez=zeros(Nx,Ny,Nz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FDTD process


i=2:Nx-1;
j=2:Ny-1;
k=2:Nz-1;
       ie=1:Nx-1;
       ih=2:Nx;
            je=1:Ny-1;
            jh=2:Ny;
                ke=1:Nz-1;
                kh=2:Nz;
for n=1:numel(Time_vec)

   n;
    Hx(ih,jh,kh)=Hx(ih,jh,kh)-alfa(ih,jh,kh).*(Ey(ie,je,ke+1)-Ey(ie,je,ke)-Ez(ie,je+1,ke)+Ez(ie,je,ke));
    Hy(ih,jh,kh)=Hy(ih,jh,kh)-alfa(ih,jh,kh).*(Ez(ie+1,je,ke)-Ez(ie,je,ke)-Ex(ie,je,ke+1)+Ex(ie,je,ke));
    Hz(ih,jh,kh)=Hz(ih,jh,kh)-alfa(ih,jh,kh).*(Ex(ie,je+1,ke)-Ex(ie,je,ke)-Ey(ie+1,je,ke)+Ey(ie,je,ke));
    
    
    Ex(ie,je,ke)=(Ex(ie,je,ke).*beta(ie,je,ke))-(gamma(ie,je,ke).*(Hz(ih,jh,kh)-Hz(ih,jh-1,kh)-Hy(ih,jh,kh)+Hy(ih,jh,kh-1)));
    Ey(ie,je,ke)=(Ey(ie,je,ke).*beta(ie,je,ke))-(gamma(ie,je,ke).*(Hz(ih,jh,kh)-Hz(ih-1,jh,kh)-Hx(ih,jh,kh)+Hx(ih,jh,kh-1)));
    Ez(ih,je,ke)=(Ez(ih,je,ke).*beta(ih,je,ke))-(gamma(ih,je,ke).*(Hy(ih,jh,kh)-Hy(ih-1,jh,kh)-Hx(ih,jh,kh)+Hx(ih,jh-1,kh)));
%%%%source

Ex(ie,je,num_node_air_z/2)=Ex(ie,je,num_node_air_z/2)+E_0*cos(2*pi*f*(n)*delta_t);
%% truncation


%%%%%%%%%%%at z=0 and end
Ex_save1{n}=Ex((i-1),j,[2 end-1]);
Ex_save2{n}=Ex(i,j,[2 end-1]);
Ex_save3{n}=Ex(i+1,j,[2 end-1]);
Ey_save1{n}=Ey(i,j-1,[2 end-1]);
Ey_save2{n}=Ey(i,j,[2 end-1]);
Ey_save3{n}=Ey(i,j+1,[2 end-1]);
Ex1{n}=Ex(end-1,j,[2 end-1]);
Ex11{n}=Ex(1,j,[2 end-1]);
Ex2{n}=Ex(end,j,[2 end-1]);
Ex22{n}=Ex(2,j,[2 end-1]);
Ey1{n}=Ey(i,end,[2 end-1]);
Ey11{n}=Ey(i,1,[2 end-1]);
Ey2{n}=Ey(i,end-1,[2 end-1]);
Ey22{n}=Ey(i,2,[2 end-1]);




%%%%%%%%%%%at y=0 and y=end
Ex_save4{n}=Ez(i-1,[2 end-1],k);
Ex_save5{n}=Ez(i,[2 end-1],k);
Ex_save6{n}=Ez(i+1,[2 end-1],k);
Ez_save1{n}=Ez(i,[2 end-1],k-1);
Ez_save2{n}=Ez(i,[2 end-1],k);
Ez_save3{n}=Ez(i,[2 end-1],k+1);
Ex3{n}=Ex(end-1,[2 end-1],k);
Ex33{n}=Ex(1,[2 end-1],k);
Ex4{n}=Ex(end-1,[2 end-1],k);
Ex44{n}=Ex(2,[2 end-1],k);
Ez1{n}=Ez(i,[2 end-1],end-1);
Ez11{n}=Ez(i,[2 end-1],1);
Ez2{n}=Ez(i,[2 end-1],end);
Ez22{n}=Ez(i,[2 end-1],2);


%%%%%%%%%%%%%%at x=0 and end
Hy_save1{n}=Hy([1 end],j-1,k);
Hy_save2{n}=Hy([1 end],j,k);
Hy_save3{n}=Hy([1 end],j+1,k);
Hz_save1{n}=Hz([1 end],j,k-1);
Hz_save2{n}=Hz([1 end],j,k);
Hz_save3{n}=Hz([1 end],j,k+1);
Hy1{n}=Hy([2 end-1],1,k);
Hy11{n}=Hy([2 end-1],end,k);
Hy2{n}=Hy([2 end-1],1+1,k);
Hy22{n}=Hy([2 end-1],end-1,k);
Hz2{n}=Hz([2 end-1],j,1);
Hz22{n}=Hz([2 end-1],j,end);
Hz1{n}=Hz([2 end-1],j,1+1);
Hz11{n}=Hz([2 end-1],j,end-1);








%%%%%%%%TRUNCATION SOFT_LATTICE
if n>m
   
    
    %at z=0 and z=end
     Ex(i,j,[1 end])=(1/3)*(Ex_save1{n-m}+ Ex_save2{n-m}+ Ex_save2{n-m});
     Ey(i,j,[1 end])=(1/3)*(Ey_save3{n-m}+ Ey_save2{n-m}+ Ey_save3{n-m});
     Ex(end,j,[1 end])=.5*(Ex1{n-m}+Ex2{n-m});
     Ey(i,end,[1 end])=.5*(Ey1{n-m}+Ey2{n-m});
     Ex(1,j,[1 end])=.5*(Ex11{n-m}+Ex22{n-m});
     Ey(i,1,[1 end])=.5*(Ey11{n-m}+Ey22{n-m});
     %at y=0 and y=end
     Ex(i,[1 end],k)=(1/3)*(Ex_save4{n-m}+Ex_save5{n-m}+Ex_save6{n-m});
     Ez(i,[1 end],k)=(1/3)*(Ez_save1{n-m}+Ez_save2{n-m}+Ez_save3{n-m});
     Ex(end,[1 end],k)=.5*(Ex3{n-m}+Ex4{n-m});
     Ez(i,[1 end],end)=.5*(Ez1{n-m}+Ez2{n-m});
     Ex(1,[1 end],k)=.5*(Ex33{n-m}+Ex44{n-m});
     Ez(i,[1 end],1)=.5*(Ez11{n-m}+Ez22{n-m});     
     %at x=0 and x=end
     Hy([1 end],j,k)=(1/3)*(Hy_save1{n-m}+Hy_save2{n-m}+Hy_save3{n-m});
     Hz([1 end],j,k)=(1/3)*(Hz_save1{n-m}+Hz_save2{n-m}+Hz_save3{n-m});
     Hy([2 end-1],1,k)=.5*(Hy1{n-m}+Hy2{n-m});
     Hz([2 end-1],j,1)=.5*(Hz1{n-m}+Hz2{n-m});
     Hy([2 end-1],end,k)=.5*(Hy11{n-m}+Hy22{n-m});
     Hz([2 end-1],j,end)=.5*(Hz11{n-m}+Hz22{n-m});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of Truncation


%%
clc;
selection=abs(Ex(:,Ny/2,num_node_air_z/2));
n

%%Plotting


plot(selection)
ylabel(' Ez(0,0,0)(V/m)')
grid on;

title([' Ez at n=' num2str(n) ])

pause(.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%