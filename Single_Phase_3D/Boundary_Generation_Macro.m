

clear;
lx=30;
ly=60;
lz=3;

A=zeros(lx,ly,lz);

for k=1:lz
    for j=1:ly
        for i=1:lx
           %if (((k==1) || (k==lz)) && (abs(i-lx/2)<=10) && (abs(j-ly/2)<=10)) 
           %if ((k==12) && (abs(i-lx/2)<=9) && (abs(j-ly/2)<=9))
            %if (sqrt((i-lx*0.3)*(i-lx*0.3)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=6)
           %if ((j==13) && (i>=5) && (i<=lx-5) && (k>=5) && (k<=lz-5))
           if ((j==1) || (j==ly))
           %if (j==1)
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
B=zeros(lx,ly,lz);

% for k=1:lz
%     for j=1:ly
%         for i=1:lx
%             %if (sqrt((i-lx*0.15)*(i-lx*0.15)+(j-ly*0.5)*(j-ly*0.5)+(k-lz*0.5)*(k-lz*0.5))<=7)
%              %if ((abs(i-lx/2)<=6) && (abs(j-ly/2)<=6))
%              %if ((abs(i-lx/2)<=7) && (abs(j-ly/2)<=7) && (abs(k-lz/2)<=7))
%              %if (j>=5) && (j<=18)
%              %if ((sqrt((i-20)*(i-20)+(j-8)*(j-8))<=7) ||(sqrt((i-6)*(i-6)+(j-8)*(j-8))<=7)) %100,60,3
%                 A(i,j,k)=0.8;B(i,j,k)=0.7;
%             else
%                 A(i,j,k)=0.2;B(i,j,k)=0.1;
%             end
%             
%         end
%     end
% end


%=======================A=Permeability,B=Porosity================

for k=1:lz
    for j=1:ly
        for i=1:lx
           
                A(i,j,k)=0.21;B(i,j,k)=0.1;
         
            
        end
    end
end





fid = fopen('perm.dat','wt');

for k=1:lz
    for j=1:ly
        for i=1:lx
        fprintf(fid,'%1d %1d\n',A(i,j,k),B(i,j,k));
        end
    end
end
