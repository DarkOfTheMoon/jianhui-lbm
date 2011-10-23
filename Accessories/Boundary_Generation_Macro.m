dat_general=[];str=[];xs=[];Disp=[];EX_general=[];DX_general=[];
dat_ori=[];
n=460;
v=2.4e-03;
dt=1;

D_ana=0.0013+(6^2*v^2)/210/0.0013;

for i=1:n
Pic_Num=5000*i;
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
        if xx(j)+td>nx-1
                
            xx(j)=(xx(j)+td)-(nx-1);
            
        else
            xx(j)=xx(j)+td;
        end
    end
  
   for j=1:nx
      
    
        %if (xx(k)>xx(k+1))
        %    t=xx(k);xx(k)=xx(k+1);xx(k+1)=t;
        %    t=dat(k);dat(k)=dat(k+1);dat(k+1)=t;
        %end
        dat(j)=dat1(round(xx(j)));
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

%plot(dat_general);
%legend(str);


 
 



