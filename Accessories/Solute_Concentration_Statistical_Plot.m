dat_general=[];str=[];xs=[];Disp=[];
n=10;
v=5e-3;
dt=1;


for i=1:n
Pic_Num=1000*i;
    dat=load(strcat('Statistical_data_concentration_X_',int2str(Pic_Num),'.sta'));
    [nx,ny]=size(dat);

xx=0:nx-1;
    for j=1:nx
        if xx(j)-v*dt*Pic_Num<0
            xx(j)=(nx-1)+(xx(j)-v*dt*Pic_Num);
        else
            xx(j)=xx(j)-v*dt*Pic_Num;
        end
    end
    
for j=1:nx-1
    for k=1:nx-j
        if (xx(k)>xx(k+1))
            t=xx(k);xx(k)=xx(k+1);xx(k+1)=t;
            t=dat(k);dat(k)=dat(k+1);dat(k+1)=t;
        end
    end
end

    
xs=[xs,xx'];    

str=[str;Pic_Num];
dat_general=[dat_general,dat];

dat=dat/(sum(dat));
xx=0:nx-1;
EX=xx*dat;
DX=xx.^2*dat-EX^2;
Disp=[Disp,DX*0.5/Pic_Num];

end

x=0:nx-1;
str=num2str(str);

plot(xs,dat_general);
legend(str);


 
 



