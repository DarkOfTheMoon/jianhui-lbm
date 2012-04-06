rel1=load('Relative_Permeability_Component1.txt');
rel2=load('Relative_Permeability_Component2.txt');
len=50;

[nx,ny]=size(rel1);

rel=zeros(nx,3);


for i=len:nx
    sumt1=0;
    sumt2=0;
    sumy1=0;
    sumy2=0;
    for j=0:len-1

    sumt1=sumt1+len-j;
    sumt2=sumt2+len-j;
    sumy1=sumy1+rel1(i-j,8);
    sumy2=sumy2+rel2(i-j,5);
    end
    
    sumt1=sumt1/len;
    sumt2=sumt2/len;
    sumy1=sumy1/len;
    sumy2=sumy2/len;
    
    x1_1=0;
    x1_2=0;
    x0_1=0;
    x0_2=0;
    t2=0;
    
    for j=0:len-1
        x1_1=x1_1+(rel1(i-j,8)-sumy1)*(len-j-sumt1);
        x1_2=x1_2+(rel2(i-j,5)-sumy2)*(len-j-sumt2);
        t2=t2+(len-j-sumt1)*(len-j-sumt1);
    end
    x1_1=x1_1/t2;
    x1_2=x1_2/t2;
    
    x0_1=sumy1-x1_1*sumt1;
    x0_2=sumy2-x1_2*sumt2;
    
    rel(i,2)=x0_1+len*x1_1;
    rel(i,3)=x0_2+len*x1_2;
    rel(i,1)=sumy2;
    
end


for i=1:len-1
    rel(i,2)=rel1(i,8);
    rel(i,3)=rel2(i,5);
   % rel(i,1)=rel1(i,2);
end

% for i=1:nx
%     rel(i,1)=rel1(i,2);
% end




figure
plot(rel(:,2),'r');
hold on
plot(rel1(:,8),'b');
hold off

figure
plot(rel(:,3),'r');
hold on
plot(rel2(:,5),'b');
hold off
