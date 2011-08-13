%Lattice Boltzmann Code:
%ostringstream name2;
%	name2<<"LBM_velocity_"<<m<<".out";
%	ofstream out2(name2.str().c_str());
%	for (int j=0;j<NY0;j++)
%		out2<<u[2][j][1][0]<<endl;

Pic_Num=20000;
dat=load(strcat('LBM_velocity_',int2str(Pic_Num),'.out'));

[nx,ny]=size(dat);
s_v=0.04;
nu=1.0/(3*s_v)-1.0/6.0;
nu1=0.2;
nu2=0.02
aa=8.5;
nu=0.02;
a=(nx-1)/2;
g=1.0e-6;
x=1:0.01:nx;
a2=(nx-2)/2;
[lal,lbl]=size(x);
ui=zeros(1,lbl);
ind=1;
for xl=1:0.01:nx
    xl=xl-(nx+1)/2;
    if (abs(xl)<=aa)
        ui(ind)=g/(2*nu1)*(a2*a2-aa*aa)+g/(2*nu2)*(aa*aa-xl*xl);
        ind=ind+1;
    else
        ui(ind)=g/(2*nu1)*(a2*a2-xl*xl);
        ind=ind+1;
    end
end


x=1:0.01:nx;
plot(x,ui,'r');
hold on
plot(dat,'kx');
xlabel('coordinate of Y');
ylabel('Velocity');
title(strcat('Channel width ',int2str(nx)));
legend('Anlytical Solution Half-way','Simulation');
hold off

