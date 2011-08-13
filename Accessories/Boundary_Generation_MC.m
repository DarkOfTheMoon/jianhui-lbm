

clear;
lx=20;
ly=40;
lz=3;

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
           %if (((k==1) || (k==lz)) && (abs(i-lx/2)<=10) && (abs(j-ly/2)<=10)) 
           %if ((k==12) && (abs(i-lx/2)<=9) && (abs(j-ly/2)<=9))
           %if (sqrt((i-lx*0.3)*(i-lx*0.3)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=6)
           if ((j==1) || (j==ly))
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

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
            %if (sqrt((i-lx*0.15)*(i-lx*0.15)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=7)
             %if ((abs(i-lx/2)<=5) && (abs(j-ly/2)<=5) && (abs(k-lz/2)<=5))
             if (j<=10) || (j>ly-10)
                A(i,j,k)=1;
            else
                A(i,j,k)=-1;
            end
            
        end
    end
end

fid = fopen('phase.dat','wt');

for k=1:lz
    for j=1:ly
        for i=1:lx
        fprintf(fid,'%1d\n',A(i,j,k));
        end
    end
end
