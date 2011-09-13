

Pic_Num=500;

dat=load(strcat('LBM_psi_',int2str(Pic_Num),'.out'));

sum=trapz(dat(:,1),dat(:,2));

da=[];
for i=1:50
    ls=i*100;
    dat=load(strcat('LBM_psi_',int2str(ls),'.out'));
    sum=trapz(dat(:,1),dat(:,2));
    da=[da,sum];
    
    
end
