%===========   Reynolds number Known,Real physical Size 3D======================================= 
% dx=1;                                                           %**********************
% dt=1;                                                           %**********************
% 
% %--------------------------------------------------
% % L_r=1;          %Real phyical length (m)                        **********************
% % u_max_r=10;             %Real physical maximum speed (m/s)      **********************
% % nx=100;                         %NX in LBM        
% % delta_x=L_r/(nx*dx);                    %length  m/lu    resolution
% %----------------------------------------------------
% 
% delta_x=1e-6/dx;                   %Resolution in (m)
% nx=100;                     
% %----------------------------------------------------
% 
% 
% 
% rho_r=1000;             %Real density (kg/m^3)                  **********************
% Re=1000;            %Reynolds number                             **********************  
% G_r=9.8;                      %Real Gravity (m/s^2)             ************************
% sigma_r=0.0238;                 %Real surface tension (N/m)         ********************        
% 
% 
% 
% rho_s=1;                %lattice density
% 
% delta_rho=rho_r/rho_s;
% delta_t=0.1*delta_x/u_max_r;              %time    s/lt
% delta_kg=rho_r*delta_x^3/rho_s;
% 
% 
% l_dx_p=delta_x*dx;                      %length   every dx in real physical unit
% l_dt_p=delta_t*dt;                      %time       every dt in real physical unit
% niu_s=0.1*nx/Re;                        %lattice viscosity
% G_s=G_r/(delta_x/(delta_t^2));          %Gravity in lattice 
% %G_s=G_r/delta_rho;
% delta_pressure=101325/(1/3);
% sigma_s=sigma_r/(delta_pressure*delta_x);


%========================================================================================





%=============    Real viscosity known  ,Real physical Size, Re number====================  

dx=1;                                                           %**********************
dt=1;                                                           %**********************

%--------------------------------------
% L_r=1;          %Real phyical length (m)                        **********************
% nx=100;                         %NX in LBM        
% delta_x=L_r/(nx*dx);                    %length  m/lu    resolution
%----------------------------------------------
delta_x=4.9e-6/dx;                   %Resolution in (m)
%-------------------------------------------

rho_r=1000;             %Real density (kg/m^3)                  ********************** 
G_r=9.8;                      %Real Gravity (m/s^2)             ************************
sigma_r=0.0238;                 %Real surface tension (N/m)         ********************        
niu_r=1e-6;                     %Real viscosity(m/s^2)          *********************
niu_s=0.016;                     %lattice viscosity       ***********************


rho_s=1;                %lattice density

delta_t=niu_s*delta_x*delta_x/niu_r;              %time    s/lt
delta_kg=rho_r*delta_x^3/rho_s;


l_dx_p=delta_x*dx;                      %length   every dx in real physical unit
l_dt_p=delta_t*dt;                      %time       every dt in real physical unit
G_s=G_r/(delta_x/(delta_t^2));          %Gravity in lattice 
delta_pressure=101325/(1/3);

sigma_s=sigma_r/(delta_pressure*delta_x);







