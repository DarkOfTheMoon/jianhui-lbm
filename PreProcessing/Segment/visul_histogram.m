a1=load ('histo_test.dat');
a2=load ('histo_test2.dat');
a3=load ('histo_test3.dat');
a4=load ('histo_test4.dat');



[sa,sb]=size(a1);
x=1:sa;


figure
plot (a1);

figure
plot (a2);
figure
plot (a3);

% figure
% plotyy(x,a2,x,a3);

figure
plot (a4);
