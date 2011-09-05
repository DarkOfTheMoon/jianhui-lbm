por=0.12;
L=32;
n_tube=1;%1:1,2:2,3:4;
N_t=[1,2,4];
c_area=L*L*por;
area_per=c_area/N_t(n_tube);
radius=sqrt(area_per/pi);
radius=9.0;

center=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5;0.25,0.75,0.25,0.75,0.75,0.25,0.75,0.25;0.25,0.25,0.25,0.75,0.75,0.25,0.75,0.75];

Lx=12;





A=ones(Lx,L,L);

for k=1:L
    for j=1:L
        for i=1:Lx
             for s=1:4
             jc=center(n_tube,(s-1)*2+1)*L;
             kc=center(n_tube,(s-1)*2+2)*L;
             if (sqrt((j-jc)*(j-jc)+(k-kc)*(k-kc))<=radius)
                 A(i,j,k)=0;
             end
             
             end
             
        end
    end
end

 
fid = fopen('Tube.dat','wt');

for k=1:L
    for j=1:L
        for i=1:Lx
        fprintf(fid,'%1d\n',A(i,j,k));
        end
    end
end


perm=pi*radius*radius*radius*radius/8/L/L*N_t(n_tube)
radius

