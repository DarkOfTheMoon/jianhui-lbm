

clear;
lx=150;
ly=120;
lz=3;
la=143/2;
lb=117/2;

ci=(lx+1)/2;
cj=(ly+1)/2;
A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
            if (j-cj)*(j-cj)/(lb*lb)+(i-ci)*(i-ci)/(la*la)<=1
                A(i,j,k)=1;
            end
            %if (j==29)
            %    A(i,j,k)=0;
            %end
            
        end
    end
end




 
fid = fopen('fibrous.dat','wt');

for k=1:lz
    for j=1:ly
        for i=1:lx
        fprintf(fid,'%1d\n',A(i,j,k));
        end
    end
end



x=0:0.05:lx;
[sx1,sx2]=size(x);
h=zeros(1,sx2);tmp=0;
for i=1:sx2
    if (abs(x(i)-ci)<=la)
    h(i)=ly/2-(sqrt(1-(x(i)-ci)*(x(i)-ci)/(la*la))*lb);h(i)=3/(h(i)^3);tmp=tmp+h(1,i)*(lx/sx2);
    else
        h(i)=ly/2;h(i)=3/(h(i)^3);tmp=tmp+h(1,i)*(lx/sx2);%h(i)=0;tmp+=h(i);
    end
   
end

x=x-lx/2;
I=trapz(x,h);
K=2*lx/ly/I/la/lb
por=1-pi/4*(la*lb*4/lx/ly)

