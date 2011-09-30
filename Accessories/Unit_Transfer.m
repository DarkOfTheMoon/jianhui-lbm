c_s_r=332.532;                %REAL SOUND SPEED  (m/s)
niu_r=1.5e-5;               %REAL VISCOSITY (m^2/s)
rho_r=1000;                    %REAL DENSITY (kg/m^3)


rho_s=1.0;                      %Lattice Density
niu_s=0.1;                      %Lattice Viscosity
dx=1.0;                         %Lattice size
dt=1.0;                         %Time Size
c_s_s=1/sqrt(3);                          %lattice sound speed



%=================================================================
rho_ref=rho_r/rho_s;
u_ref=c_s_r/c_s_s;
L_ref=niu_r/(u_ref*niu_s);       

t_r=dt*L_ref/u_ref;              %dt in real physics units

