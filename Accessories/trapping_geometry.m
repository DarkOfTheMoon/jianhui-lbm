
ny=40;
nx=97;
a=zeros(nx,ny,ny);
min=6;
half=((ny/2)-1-min)/2-1;
l=ny/2-1;
for i=1:nx
    for j=1:ny
        for k=1:ny
            
%============TRIANGLE=============================== 
        l=min+half*(sin(0.12*i-pi/2*0.8))+half;     
        if ((j-ny/2)*(-sqrt(3))+l<=k-ny/2)
                a(i,j,k)=1;
        end
        if (sqrt(3)*(j-ny/2)+l<=k-ny/2)
                a(i,j,k)=1;
       end
        if (k-ny/2<=-l)
               a(i,j,k)=1;
       end
%================================================== 

%=========SQURE====================================
%        l=min+half*(sin(0.15*i))+half; 
%        if (-(j-ny/2)+l<=k-ny/2)
%            a(i,j,k)=1;
%        end
%        if (-(j-ny/2)-l>=k-ny/2)
%            a(i,j,k)=1;
%        end
%        if ((j-ny/2)+l<=k-ny/2)
%            a(i,j,k)=1;
%        end
%        if ((j-ny/2)-l>=k-ny/2)
%            a(i,j,k)=1;
%        end
%==================================================
    
        end
    end
    
end

fid = fopen('Capillary_Trapping.vtk','wt');
% fprintf(fid,'# vtk DataFile Version 2.0\n');
% fprintf(fid,'J.Yang Lattice Boltzmann Simulation 3D Capillary Trapping\n');
% fprintf(fid,'ASCII\n');
% fprintf(fid,'DATASET STRUCTURED_POINTS\n');
% fprintf(fid,'DIMENSIONS         %i         %i         %i\n',nx,ny,ny);
% fprintf(fid,'ORIGIN 0 0 0\n');
% fprintf(fid,'SPACING 1 1 1\n');
% fprintf(fid,'POINT_DATA     %i\n',nx*ny*ny);
% fprintf(fid,'SCALARS sample_scalars float\n');
% fprintf(fid,'LOOKUP_TABLE default\n');
% 
for k=1:ny
    for j=1:ny
        for i=1:nx
        fprintf(fid,'%1d\n',a(i,j,k));
        end
    end
end



