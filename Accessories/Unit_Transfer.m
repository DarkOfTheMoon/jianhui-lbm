c_s_r=332.532;                %REAL SOUND SPEED  (m/s)
niu_r=9.2e-4;               %REAL VISCOSITY (m^2/s)
rho_r=730;                    %REAL DENSITY (kg/m^3)


rho_s=1.0;                      %Lattice Density
dx=1.0;                         %Lattice size
dt=1.0;                         %Time Size
c_s_s=dx/dt/sqrt(3);                          %lattice sound speed
sigma_rel=0.0238e-3;               %Real surface tension
F_rel=9.8e-3;                      %Real Force (N)

%=================================================================
% niu_s=0.1;                      %Lattice Viscosity
% rho_ref=rho_r/rho_s;
% u_ref=c_s_r/c_s_s;
% L_ref=niu_r/(u_ref*niu_s);       
% 
% t_r=dt*L_ref/u_ref;              %dt in real physics units

%=================================================================

L_rel=24e-6;                    %RESOLUTION IN REAL UNIT (m)
rho_ref=rho_r/rho_s;
u_ref=c_s_r/c_s_s;
L_ref=L_rel;
niu_s=niu_r/(u_ref*L_ref);
t_r=dt*L_ref/u_ref;              %dt in real physics units
%=================================================================

t_ref=t_r/dt;

F_ref=rho_ref*(L_ref^4)/(t_ref*t_ref);
F_s=F_rel*F_ref;

sigma_s=sigma_rel*(F_ref/L_ref);