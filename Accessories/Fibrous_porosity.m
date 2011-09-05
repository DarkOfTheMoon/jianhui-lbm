lx=200;
ly=200;

la=(lx/2-2);
lb=(ly/2-1);



ci=(lx+1)/2;
cj=(ly+1)/2;

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

la=la*2;
lb=lb*2;