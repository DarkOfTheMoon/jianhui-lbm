Ma=0.1;                        %Mach number
Pr=0.71;                        %Prandtl number
H=64;                          %Scale of lattice
Ra=100000;                        %Rayleigh number
dt=1;


tau_f=0.5+(Ma*H*sqrt(3*Pr))/(dt*sqrt(Ra));
niu_f=1/3*(tau_f-0.5)*dt;
yix=niu_f/Pr;



niu_f
yix
