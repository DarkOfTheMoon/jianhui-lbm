c_s_r=1497.532;                %REAL SOUND SPEED  (m/s)
niu_r=1.0e-6;               %REAL VISCOSITY (m^2/s)  KINETIC VISCOSITY
rho_r=1000;                    %REAL DENSITY (kg/m^3)


rho_s=1000.0;                      %Lattice Density
dx=3.0;                         %Lattice size
dt=9.0;                         %Time Size
c_s_s=dx/dt/sqrt(3);                          %lattice sound speed
sigma_rel=0.0238e-3;               %Real surface tension
F_rel=9.8;                      %Real Force (N)

%=================================================================
% niu_s=0.1;                      %Lattice Viscosity
% rho_ref=rho_r/rho_s;
% u_ref=c_s_r/c_s_s;
% L_ref=niu_r/(u_ref*niu_s);       
% 
% t_r=dt*L_ref/u_ref;              %dt in real physics units
%  L_r=dx*L_ref;                     %resolution in physical units
%=================================================================

L_rel=4.9e-6;                    %RESOLUTION IN REAL UNIT (m)
rho_ref=rho_r/rho_s;
u_ref=c_s_r/c_s_s;
L_ref=L_rel/dx;
niu_s=niu_r/(u_ref*L_ref);
t_r=dt*L_ref/u_ref;              %dt in real physics units
%=================================================================

t_ref=t_r/dt;

F_ref=rho_ref*(L_ref^4)/(t_ref*t_ref);
F_s=F_rel*F_ref;

sigma_s=sigma_rel*(F_ref/L_ref);