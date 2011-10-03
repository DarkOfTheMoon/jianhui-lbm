%Lattice Boltzmann Code:
%ostringstream name2;
%	name2<<"LBM_velocity_"<<m<<".out";
%	ofstream out2(name2.str().c_str());
%	for (int j=0;j<NY0;j++)
%		out2<<u[2][j][1][0]<<endl;

Pic_Num=4200;
dat=load(strcat('LBM_velocity_',int2str(Pic_Num),'.out'));
dx=0.5;


[nx,ny]=size(dat);
nx_dat=0:dx:nx*dx-dx;

nx=(nx-1)*dx;

s_v=0.04;
nu=1.0/(3*s_v)-1.0/6.0;
nu=0.02;
a=(nx-dx)/2;

g=1.0e-6;
x=-a:0.01:a;
u=g/(2*nu)*(a*a-x.*x);

% x=1:0.01:nx;
x=x+nx/2;
plot(x,u,'r');
hold on
plot(nx_dat,dat,'kx');
xlabel('coordinate of Y');
ylabel('Velocity');
title(strcat('Channel width ',int2str(nx)));
legend('Anlytical Solution Half-way','Simulation');
hold off

% ana=g/(2*nu)*(a*a*(2*a)-1/3*2*2^3);
% 
% sim=0;
% for i=1:nx
%     sim=sim+dat(i);
% end
% error=abs(sim-ana)/ana
% error2=abs(sim-ana2)/ana2

