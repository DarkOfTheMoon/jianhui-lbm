A=load('2_General_disp_concentration_X_2000000.sta');
n_0=350;
u_ave=3.81794e-06;
t=2600000;

u_0_b=u_ave*t;






[nx,ny]=size(A);
new_x=-(nx-n_0+1):n_0-2;
na=[];

for (i=n_0:nx)
    na=[na,A(i)];
end

for (i=1:n_0-1)
    na=[na,A(i)];
end

new_x_rescal=new_x/u_0_b;
na_rescal=na*u_0_b;


% plot (new_x,na);
% figure

plot (new_x_rescal,na_rescal);
xlabel('\zeta/ \langle \zeta \rangle_0');
ylabel('P(\zeta)x\langle \zeta \rangle_0')

