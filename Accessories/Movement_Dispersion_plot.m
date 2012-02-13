A=load('2_General_disp_concentration_X_2600000.sta');
n_0=350;

[nx,ny]=size(A);
new_x=-(nx-n_0+1):n_0-2;
na=[];

for (i=n_0:nx)
    na=[na,A(i)];
end

for (i=1:n_0-1)
    na=[na,A(i)];
end




plot (new_x,na);
