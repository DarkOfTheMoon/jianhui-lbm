nx_ori=20;
x=[];
for i=-0.5:1:nx_ori-2+0.5;
    x=[x;i];
end

[nx,ny]=size(x);
x(1)=0;
x(nx)=nx_ori-2;

G=1e-6;
m=0.001;
n=0.5;
L=nx_ori-2;
y=[];
for i=1:nx
    %y=[y;G/(2*m)*((L/2)^2-(L/2-x(i))^2)];                    %NF
    %y=[y;(G/(2*m))^(1/n)*(n/(n+1))*((L/2)^((n+1)/n)-(abs(L/2-x(i)))^((n+1)/n))];  %NNF
    %y=[y;(G/(2*m))^(1/n)*(n/(n+1))*((L/2)^((n+1)/n)-(abs(L/2-x(i)))^((n+1)/n))];  %NNF
    y=[y;(2*n+1)/(n+1)*(1-((abs(L/2-x(i)))/(L/2))^(1+1/n))]; %%NNF AVERAGE
end

%y=y/max(y);
plot (x,y);


%==========================================
hold on
dat=load('LBM_velocity_12000.out');
[nx,ny]=size(dat);

%====MAX NORMALISE==========
%dat=dat/max(dat);
%===========================

%=====AVE NORMALISE=========
sums=sum(dat);
sums=sums/(nx-2);
dat=dat/sums;
%===========================


xx=[];
for i=-0.5:1:nx-1-0.5;
    xx=[xx;i];
end
plot (xx,dat,'rx');
hold off
% %============================================

%===============error cal====================
error=0;
sums=0;
for i=2:nx-1
    error=error+abs(y(i)-dat(i));
    sums=sums+abs(y(i))
end

error=error/sums
%============================================
