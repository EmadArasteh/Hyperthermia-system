clc;
close all;
clearvars;
%%initialization
c=3e8;
T_simulation=2400;
%%%%%%%%%%%%%%%%%%%%%%MESH
len_body_x=300e-3;
len_body_y=len_body_x;
len_body_z=len_body_x;
len_air_x=30e-3;
len_air_y=len_air_x;
len_air_z=len_air_x;

len_tot_x=len_air_x+len_body_x;
len_tot_y=len_air_z+len_body_y;
len_tot_z=len_air_z+len_body_z;

Nx=40;
Ny=Nx;
Nz=Ny;


delta_x=(len_tot_x/Nx);
delta_y=(len_tot_y/Ny);
delta_z=(len_tot_z/Nz);
delta=delta_x;

num_node_air_x=round(len_air_x/delta_x);
num_node_air_y=round(len_air_y/delta_y);
num_node_air_z=round(len_air_z/delta_z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%import

f=1e9;
Cp_fat=2300;
Cp_min=Cp_fat;
K_fat=.22;
K_max=K_fat;
Rho_fat=920;
Rho_min=Rho_fat;
eps_r_fat=10;
eps_max=eps_r_fat;
sigma_fat=0.17;
b_fat=816;
b_max=b_fat;
Tb=37;
Ta=20;
ha=10.5;
hb=50;
K_air=0;
Cp_air=1;
Rho_air=1;
b_air=0;


%%%%%%%%%%%%MATRIX of Import
SAR=zeros(Nx,Ny,Nz);

T=ones(Nx,Ny,Nz)*Tb;
T(1:num_node_air_x/2,:,:)=Ta;
T(end:end-num_node_air_x/2+1,:,:)=Ta;




K=zeros(Nx,Ny,Nz);
K(:,:,:)=K_fat;
K(1:num_node_air_x/2,:,:)=K_air;
K(end:end-num_node_air_x/2+1,:,:)=K_air;
K(:,1:num_node_air_y/2,:)=K_air;
K(:,end:end-num_node_air_y/2+1,:)=K_air;
K(:,:,1:num_node_air_z/2)=K_air;
K(:,:,end:end-num_node_air_z/2+1)=K_air;




Cp=zeros(Nx,Ny,Nz);
Cp(:,:,:)=Cp_fat;
Cp(1:num_node_air_x/2,:,:)=Cp_air;
Cp(end:end-num_node_air_x/2+1,:,:)=Cp_air;
Cp(:,1:num_node_air_y/2,:)=Cp_air;
Cp(:,end:end-num_node_air_y/2+1,:)=Cp_air;
Cp(:,:,1:num_node_air_z/2)=Cp_air;
Cp(:,:,end:end-num_node_air_z/2+1)=Cp_air;



Rho=zeros(Nx,Ny,Nz);
Rho(:,:,:)=Rho_fat;
Rho(1:num_node_air_x/2,:,:)=Rho_air;
Rho(end:end-num_node_air_x/2+1,:,:)=Rho_air;
Rho(:,1:num_node_air_y/2,:)=Rho_air;
Rho(:,end:end-num_node_air_y/2+1,:)=Rho_air;
Rho(:,:,1:num_node_air_z/2)=Rho_air;
Rho(:,:,end:end-num_node_air_z/2+1)=Rho_air;


b=zeros(Nx,Ny,Nz);
b(:,:,:)=b_fat;
b(1:num_node_air_x/2,:,:)=b_air;
b(end:end-num_node_air_x/2+1,:,:)=b_air;
b(:,1:num_node_air_y/2,:)=b_air;
b(:,end:end-num_node_air_y/2+1,:)=b_air;
b(:,:,1:num_node_air_z/2)=b_air;
b(:,:,end:end-num_node_air_z/2+1)=b_air;


h=zeros(Nx,Ny,Nz);
h(:,:,:)=hb;
h(1:num_node_air_x/2,:,:)=ha;
h(end:end-num_node_air_x/2+1,:,:)=ha;
h(:,1:num_node_air_y/2,:)=ha;
h(:,end:end-num_node_air_y/2+1,:)=ha;
h(:,:,1:num_node_air_z/2)=ha;
h(:,:,end:end-num_node_air_z/2+1)=ha;

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
SAR(i,j,k)=20*exp(-2*(((i)-Nx/2).^2+((j)-Ny/2).^2+((k)-Nz/2).^2));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%time
u_max=c/sqrt(eps_max);
von_neuman=2*Rho_min*Cp_min*(delta_x.^2)/(12*K_max+b_max*(delta_x.^2));
taflove=(((1/delta_x.^2)+(1/delta_y.^2)+(1/delta_z.^2)).^(-.5))/u_max;
delta_t=min(von_neuman,1);

%%%%%%%%%%%%%%

%%%%%%Timing vector


Time_vec=1:delta_t:T_simulation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% FDTD process

T(1,:,:)=Ta;
T(end,:,:)=Ta;
T(:,1,:)=Ta;
T(:,end,:)=Ta;
T(:,:,1)=Ta;
T(:,:,end)=Ta;

num_min_bound_x=num_node_air_x/2+1;
num_max_bound_x=Nx-num_node_air_x/2;
num_min_bound_y=num_node_air_y/2+1;
num_max_bound_y=Ny-num_node_air_y/2;
num_min_bound_z=num_node_air_z/2+1;
num_max_bound_z=Nz-num_node_air_z/2;

selection_max=zeros(numel(Time_vec,1));
%selection_x={};
       % for
       i=2:Nx-1;
       %for 
            j=2:Ny-1;
            %for
                k=2:Nz-1;
               
for n=1:numel(Time_vec)


      
    selection_max(n)=max(max(max(T)));
    
 %%%%%%%%%%%plot in time
 if mod(n,T_simulation/20)==0
    selection_x=T(:,Ny/2,Nz/2);
   
 
   %n
x=linspace(0,len_tot_x,Nx).*1e3;
plot(x,selection_x,x,SAR(:,Ny/2,Nz/2))
xlabel('x(mm)')

title(['X-oriented Temeprature({^o}c) and SAR(W/Kg) distribution at the center of the body at t=' num2str(Time_vec(n)) 'Second'])
pause(.5);   
    grid on;
    
    
 end
        
%%%%%%%%%%%%%%%%%%%%%%
   T(i,j,k)=T(i,j,k)+(SAR(i,j,k).*delta_t./Cp(i,j,k))-(delta_t.*(T(i,j,k)-Tb)./(Rho(i,j,k).*Cp(i,j,k)))+...
       (delta_t.*K(i,j,k)./(Rho(i,j,k).*Cp(i,j,k).*delta_x.^2)).*(T(i+1,j,k)+T(i-1,j,k)+T(i,j+1,k)+T(i,j-1,k)+T(i,j,k+1)+T(i,j,k-1)-6.*T(i,j,k));
    
    
            %end
        %end
    %end
    
T(num_min_bound_x,:,:)=(K(num_min_bound_x-1,:,:).*T(num_min_bound_x-1,:,:)./(K(num_min_bound_x-1,:,:)+h(num_min_bound_x-1,:,:)*delta_x))+(Ta.*h(num_min_bound_x-1,:,:)*delta_x./(K(num_min_bound_x-1,:,:)+h(num_min_bound_x-1,:,:).*delta_x));
T(:,num_min_bound_y,:)=(K(:,num_min_bound_y-1,:).*T(:,num_min_bound_y-1,:)./(K(:,num_min_bound_y-1,:)+h(:,num_min_bound_y-1,:).*delta_x))+(Ta.*h(:,num_min_bound_y-1,:).*delta_x./(K(:,num_min_bound_y-1,:)+h(:,num_min_bound_y-1,:).*delta));
T(:,:,num_min_bound_z)=(K(:,:,num_min_bound_z-1).*T(:,:,num_min_bound_z-1)./(K(:,:,num_min_bound_z-1)+h(:,:,num_min_bound_z-1).*delta_x))+(Ta.*h(:,:,num_min_bound_z-1).*delta_x./(K(:,:,num_min_bound_z-1)+h(:,:,num_min_bound_z-1).*delta));
T(num_max_bound_x,:,:)=(K(num_max_bound_x+1,:,:).*T(num_max_bound_x+1,:,:)./(K(num_max_bound_x+1,:,:)+h(num_max_bound_x+1,:,:).*delta_x))+(Ta.*h(num_max_bound_x+1,:,:).*delta_x./(K(num_max_bound_x+1,:,:)+h(num_max_bound_x+1,:,:).*delta));
T(:,num_max_bound_y,:)=(K(:,num_max_bound_y,:).*T(:,num_max_bound_y,:)./(K(:,num_max_bound_y,:)+h(:,num_max_bound_y,:).*delta_x))+(Ta.*h(:,num_max_bound_y,:).*delta_x./(K(:,num_max_bound_y,:)+h(:,num_max_bound_y,:).*delta));
T(:,:,num_max_bound_z)=(K(:,:,num_max_bound_z).*T(:,:,num_max_bound_z)./(K(:,:,num_max_bound_z)+h(:,:,num_max_bound_z).*delta_x))+(Ta.*h(:,:,num_max_bound_z).*delta_x./(K(:,:,num_max_bound_z)+h(:,:,num_max_bound_z).*delta));
    





end


%%Plotting

figure;
plot(Time_vec,selection_max)
grid on;
xlabel('t(Seconds)')
ylabel(' Temperature(^oc)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%