a1=load ('histo_test.dat');
a2=load ('histo_test2.dat');
a3=load ('histo_test3.dat');
a4=load ('histo_test4.dat');



[sa,sb]=size(a1);
x=1:sa;


% figure
% plot (a1);

% figure
% plot (a2);

% figure
% plotyy (x,a3,x,a1);

figure
plotyy(x,a2,x,a3);

figure
plotyy (x,a1-a4,x,a1);

% [nx,ny]=size(a1);
% test1=zeros(nx,1);
% test2=zeros(nx,1);
% test3=zeros(nx,1);
% 
% sum1=200000;
% for i=1:nx
%     sum1=200000;
%      for j=0:10
%          if ((i+j>=1) && (i+j<=nx) && (a1(i+j)<sum1))
%              sum1=a1(i+j);
%          end
%      end
%      test1(i)=sum1;
%      
%      
%      sum1=200000;
%      for j=-10:0    
%          if ((i+j>=1) && (i+j<=nx) && (a1(i+j)<sum1))
%              sum1=a1(i+j);
%          end
%      end
%      test2(i)=sum1;
%      
% end
% 
% for i=1:nx
%     if (i+10)<=nx
%      test3(i)=test1(i+10)-test1(i);
%     end
% end
% 
% figure
% plotyy(x,a1,x,test1);
% 
% figure
% plotyy(x,a1,x,test2);
% 
% figure
% plotyy(x,test1,x,test3);