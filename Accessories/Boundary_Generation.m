

clear;
lx=161;
ly=60;
lz=60;

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
            if (sqrt((i-lx*0.15)*(i-lx*0.15)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=7)
                A(i,j,k)=1;
            end
            
        end
    end
end




 
fid = fopen('BC.dat','wt');

for k=1:lz
    for j=1:ly
        for i=1:lx
        fprintf(fid,'%1d\n',A(i,j,k));
        end
    end
end


