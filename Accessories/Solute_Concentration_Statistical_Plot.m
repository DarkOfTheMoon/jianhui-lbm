dat_general=[];str=[];xs=[];Disp=[];EX_general=[];DX_general=[];
dat_ori=[];
n=7;
v=1.262644e-4;
dt=1;


for i=1:n
Pic_Num=30000*i;
    %if (i==0)
    %    Pic_Num=5000;
    %end
    
    dat=load(strcat('Statistical_data_concentration_X_',int2str(Pic_Num),'.sta'));
    [nx,ny]=size(dat);
dat1=dat;
dat_ori=[dat_ori,dat];

td=v*dt*Pic_Num;
while td>nx-1
    td=td-(nx-1);
end

xx=0:nx-1;
    for j=1:nx
        if xx(j)-td<0
                
            xx(j)=(nx-1)+(xx(j)-td);
            
        else
            xx(j)=xx(j)-td;
        end
    end

    
%=================================    
%    for j=1:nx
%     if xx(j)>200
%         xx(j)=xx(j)-nx+1;
%     end
%    end
%================================       
       
       
for j=1:nx-1
    for k=1:nx-j
        if (xx(k)>xx(k+1))
            t=xx(k);xx(k)=xx(k+1);xx(k+1)=t;
            t=dat(k);dat(k)=dat(k+1);dat(k+1)=t;
        end
    end
end

    
xs=[xs,xx'];    

str=[str;Pic_Num*dt];
dat_general=[dat_general,dat];

dat=dat/(sum(dat));
xx=1:nx;
EX=xx*dat;EX_general=[EX_general,EX];
DX=xx.^2*dat-EX^2;DX_general=[DX_general,DX];
Disp=[Disp,DX*0.5/(dt*Pic_Num)];

end

x=0:nx-1;
str=num2str(str);

plot(xs,dat_general);
%legend(str);


 
 



