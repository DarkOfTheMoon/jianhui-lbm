dat_general=[];str=[];
n=5;

for i=1:n
Pic_Num=200*i;
    dat=load(strcat('Statistical_data_concentration_X_',int2str(Pic_Num),'.sta'));
    [nx,ny]=size(dat);
dat_general=[dat_general,dat];

str=[str;Pic_Num];

%bar(rand(3,5));    
%plot(dat);
%hold on
end

x=1:nx;
str=num2str(str);

plot(x,dat_general);
legend(str);


 
 



