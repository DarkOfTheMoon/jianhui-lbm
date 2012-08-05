nx=30;
x=[];
for i=0:1:nx-2;
    x=[x;i];
end

[nx,ny]=size(x);

G=1e-5;
m=0.001;
n=0.25;
L=nx-1;
y=[];
for i=1:nx
    %y=[y;G/(2*m)*((L/2)^2-(L/2-x(i))^2)];                    %NF
    %y=[y;(G/(2*m))^(1/n)*(n/(n+1))*((L/2)^((n+1)/n)-(abs(L/2-x(i)))^((n+1)/n))];  %NNF
    y=[y;(G/(2*m))^(1/n)*(n/(n+1))*((L/2)^((n+1)/n)-(abs(L/2-x(i)))^((n+1)/n))];  %NNF
end

%y=y/max(y);
plot (x,y);


%==========================================
hold on
dat=load('LBM_velocity_6000.out');
%dat=dat/max(dat);

[nx,ny]=size(dat);
xx=[];
for i=-0.5:1:nx-1-0.5;
    xx=[xx;i];
end
plot (xx,dat,'rx');
hold off
%============================================

%===============error cal====================
% error=0;
% sum=0;
% for i=2:nx-1
%     error=error+abs(y(i)-dat(i));
%     sum=sum+abs(y(i))
% end
% 
% error=error/sum
%============================================
