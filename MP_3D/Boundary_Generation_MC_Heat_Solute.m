

clear;
lx=21;
ly=80;
lz=3;

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
           %if (((k==1) || (k==lz)) && (abs(i-lx/2)<=10) && (abs(j-ly/2)<=10)) 
           %if ((k==12) && (abs(i-lx/2)<=9) && (abs(j-ly/2)<=9))
            %if (sqrt((i-lx*0.3)*(i-lx*0.3)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=6)
           %if ((j==13) && (i>=5) && (i<=lx-5))% && (k>=5) && (k<=lz-5))
           %if ((j==1) || (j==ly))
           %if ((j==1) || (j==ly) || (k==1) || (k==lz))
           if ((i==1) || (i==lx)) && (j>=ly/2-5) && (j<=ly-5)
           %if ((j==1) || (j==ly) || (k==1) || (k==lz)) && (i>=lx/2-15) && (i<=lx-15)
               %if ((j==1) || (j==ly) || (k==1) || (k==lz))
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
            %if (sqrt((i-lx*0.5)*(i-lx*0.5)+(j-ly*0.4+3)*(j-ly*0.4+3))<=8)
             if ((abs(i-lx/2)<=6) && (abs(j-19)<=6))
             %if ((abs(i-lx/2)<=7) && (abs(j-ly/2)<=7) && (abs(k-lz/2)<=7))
             %if (i<=lx/2-10)
             %if ((sqrt((i-40)*(i-40)+(j-8)*(j-8))<=7) ||(sqrt((i-56)*(i-56)+(j-8)*(j-8))<=7)) %100,60,3
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


%==================================HEAT/SOLUTE==============================================

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
            %if (sqrt((i-lx*0.5)*(i-lx*0.5)+(j-ly*0.5)*(j-ly*0.5))<=13)
             %if ((abs(i-lx/2)<=8) && (abs(j-ly/2)<=8))
             %if ((abs(i-lx/2)<=7) && (abs(j-ly/2)<=7) && (abs(k-lz/2)<=7))
             %if (i<=lx/2-10)
             %if (j>=ly/2-3) && (j<=ly/2+3)
             %if ((sqrt((i-40)*(i-40)+(j-8)*(j-8))<=7) ||(sqrt((i-56)*(i-56)+(j-8)*(j-8))<=7)) %100,60,3
                A(i,j,k)=20.0+(rand-0.5)*7;
%             else
%                 A(i,j,k)=0;
%             end
            
        end
    end
end

fid = fopen('phase2.dat','wt');

for k=1:lz
    for j=1:ly
        for i=1:lx
        fprintf(fid,'%1d\n',A(i,j,k));
        end
    end
end

